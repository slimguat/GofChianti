#!/usr/bin/env python3
"""Maintainer tool: convert IDL GOFNT ``.dat`` output into the GofChianti dataset.

THIS IS NOT PART OF THE PUBLIC PACKAGE.  End users never run it.  It is used by
the maintainers who generate the CHIANTI data with IDL/SSW to:

1. parse each ``*_gofnt_v-*.dat`` file (with robustness fixes for two IDL output
   bugs, see below),
2. divide out the elemental abundance so the stored ``G(T, n_e)`` is *bare*,
3. write a compressed ``.npz`` per line,
4. copy the CHIANTI ``.abund`` files into the dataset,
5. build ``catalog.parquet`` (also copied into the package for offline use) and
   ``manifest.json`` (with SHA256 hashes),
6. optionally publish everything to a GitHub release via ``gh``.

IDL output quirks handled here
------------------------------
* **Asterisk wavelength** — older IDL used ``format='(F7.3)'`` for the header
  wavelength, which overflows to ``*******`` for 4-integer-digit wavelengths
  (e.g. 1025.723).  The true wavelength is recovered from the *filename*.
* **Abundance line-jump** — IDL's 80-column auto-wrap pushed the long abundance
  path onto the next line without a ``#``.  We re-join it and keep only the
  basename.

Usage
-----
    python maintainers/convert_dat_to_npz.py \
        --dat-dir ../gofnt \
        --abund-src /usr/local/ssw/packages/chianti/dbase/abundance \
        --out-dir ./dataset
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

# Make the package importable when running from a source checkout.
_SRC = Path(__file__).resolve().parent.parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from gofchianti.abundance import element_to_z, parse_abund_file, z_to_symbol  # noqa: E402
from gofchianti.core import ContributionFunction  # noqa: E402

# Matches e.g. "fe_12_195.119_gofnt_v-11.0.2.dat" -> ion, wl_str, version
_FNAME_RE = re.compile(
    r"^(?P<ion>[A-Za-z]+_\d+)_(?P<wl>\d+(?:\.\d+)?)_gofnt_v-(?P<ver>.+)\.dat$"
)


# --------------------------------------------------------------------------- #
# Parsing (.dat)
# --------------------------------------------------------------------------- #
def _filename_meta(path: Path) -> Tuple[Optional[str], Optional[float], Optional[str]]:
    m = _FNAME_RE.match(path.name)
    if not m:
        return None, None, None
    return m.group("ion"), float(m.group("wl")), m.group("ver")


def parse_gofnt_dat(path: str | Path) -> ContributionFunction:
    """Parse a GOFNT ``.dat`` file into a :class:`ContributionFunction`.

    The returned object reflects the file *as written*: ``abundance_applied`` is
    set from the header's "Abundance multiplication" flag (abundance is removed
    later by :func:`strip_abundance`).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    fname_ion, fname_wl, fname_ver = _filename_meta(path)

    raw_lines = [ln.rstrip("\n") for ln in path.open("r")]

    # Split into header region (before the first "De:" block) and data lines.
    first_de = next(
        (i for i, ln in enumerate(raw_lines) if ln.strip().startswith("De:")),
        len(raw_lines),
    )
    header_region = raw_lines[:first_de]
    data_lines = [ln.strip() for ln in raw_lines[first_de:] if ln.strip()]

    # Parse header key/values, re-joining the wrapped abundance path.
    header_dict: Dict[str, str] = {}
    first_header_line = ""
    last_key: Optional[str] = None
    for ln in header_region:
        s = ln.strip()
        if not s:
            continue
        if s.startswith("#"):
            text = s.lstrip("#").strip()
            if not first_header_line:
                first_header_line = text
            if ":" in text:
                key, val = (x.strip() for x in text.split(":", 1))
                header_dict[key.lower()] = val
                last_key = key.lower()
            else:
                last_key = None
        else:
            # Non-# line inside the header region == wrapped continuation,
            # in practice the abundance file path.
            if last_key == "abundance file":
                header_dict["abundance file"] = (header_dict.get("abundance file", "") + " " + s).strip()
            # otherwise ignore stray header noise

    # --- ion & wavelength, with asterisk recovery from the filename ---------
    m = re.search(r"line:\s*(\S+)\s+at\s+([\d.]+)", first_header_line, re.IGNORECASE)
    header_ion = m.group(1) if m else None
    header_wl = float(m.group(2)) if m else None

    ion = header_ion or fname_ion
    wavelength = header_wl if header_wl is not None else fname_wl
    if wavelength is None:
        raise ValueError(f"Could not determine wavelength for {path.name}")
    if ion is None:
        raise ValueError(f"Could not determine ion for {path.name}")
    if fname_ion and ion.lower() != fname_ion.lower():
        print(f"  [warn] header ion {ion!r} != filename ion {fname_ion!r} in {path.name}")

    chianti_version = header_dict.get("chianti version", "") or (fname_ver or "")

    def _to_int(key: str, default: int = 0) -> int:
        try:
            return int(float(header_dict.get(key, default)))
        except (ValueError, TypeError):
            return default

    def _to_float(key: str, default: float = float("nan")) -> float:
        val = header_dict.get(key)
        if val is None:
            return default
        mnum = re.search(r"([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)", val)
        return float(mnum.group(1)) if mnum else default

    description = header_dict.get("description", "")
    ch_lower = _to_int("lower level")
    ch_upper = _to_int("upper level")
    f_value = _to_float("f")
    A_value = _to_float("a")
    Nmin = _to_float("minn", 0.0)
    Nmax = _to_float("maxn", 0.0)

    abund_raw = header_dict.get("abundance file", "").strip()
    source_abundance_file = Path(abund_raw).name if abund_raw else ""

    try:
        abundance_applied = bool(int(header_dict.get("abundance multiplication", "0")))
    except ValueError:
        abundance_applied = False

    units_raw = header_dict.get("rows and units", "")
    um = re.search(r"G\(T,n\)\s*\(([^)]+)\)", units_raw)
    units = um.group(1) if um else "erg cm^3 s^-1 sr^-1"

    generation_date = header_dict.get("generated on", "")

    # --- numeric De:/T/G blocks --------------------------------------------
    densities: List[float] = []
    all_temps: List[np.ndarray] = []
    gofnt_rows: List[np.ndarray] = []
    i = 0
    while i < len(data_lines):
        line = data_lines[i]
        if line.startswith("De:"):
            parts = line.split()
            densities.append(float(parts[1]))
            if i + 2 >= len(data_lines):
                raise ValueError(f"Incomplete De: block in {path.name}")
            temps = np.array([float(x) for x in data_lines[i + 1].split()], dtype=float)
            gvals = np.array([float(x) for x in data_lines[i + 2].split()], dtype=float)
            all_temps.append(temps)
            gofnt_rows.append(gvals)
            i += 3
        else:
            i += 1

    if not densities:
        raise ValueError(f"No De: blocks found in {path.name}")

    # --- harmonise onto a single temperature grid --------------------------
    lengths = [len(t) for t in all_temps]
    if len(set(lengths)) == 1:
        ref_temp = all_temps[0]
        rows = []
        for tarr, garr in zip(all_temps, gofnt_rows):
            if np.array_equal(tarr, ref_temp):
                rows.append(garr)
            else:
                rows.append(interp1d(tarr, garr, bounds_error=False, fill_value=0.0)(ref_temp))
        gofnt_rows = rows
    else:
        ref_idx = lengths.index(max(Counter(lengths).keys()))
        ref_temp = all_temps[ref_idx]
        rows = []
        for tarr, garr in zip(all_temps, gofnt_rows):
            if tarr.shape == ref_temp.shape and np.array_equal(tarr, ref_temp):
                rows.append(garr)
            else:
                rows.append(interp1d(tarr, garr, bounds_error=False, fill_value=0.0)(ref_temp))
        gofnt_rows = rows

    gofnt_matrix = np.vstack(gofnt_rows)
    densities_arr = np.array(densities, dtype=float)
    temperatures_arr = np.array(ref_temp, dtype=float)

    return ContributionFunction(
        ion=ion,
        wavelength=float(wavelength),
        chianti_version=chianti_version,
        description=description,
        ch_lower_level=ch_lower,
        ch_upper_level=ch_upper,
        f_value=f_value,
        A_value=A_value,
        Nmin=Nmin,
        Nmax=Nmax,
        source_abundance_file=source_abundance_file,
        abundance_applied=abundance_applied,
        units=units,
        generation_date=generation_date,
        densities=densities_arr,
        temperature=temperatures_arr,
        gofnt_matrix=gofnt_matrix,
        filename=str(path),
    )


# --------------------------------------------------------------------------- #
# Abundance removal
# --------------------------------------------------------------------------- #
def strip_abundance(cf: ContributionFunction, abund_src: Path) -> ContributionFunction:
    """Divide out the elemental abundance so ``gofnt_matrix`` becomes bare.

    Uses the abundance file named in the source metadata, looked up under
    ``abund_src``.  If abundance was not applied, the object is returned
    unchanged.
    """
    if not cf.abundance_applied:
        return cf
    if not cf.source_abundance_file:
        raise ValueError(
            f"{cf.ion} {cf.wavelength}: abundance was applied but no abundance "
            "file is recorded; cannot divide it out."
        )
    abund_path = abund_src / cf.source_abundance_file
    if not abund_path.exists():
        raise FileNotFoundError(f"Abundance file not found for division: {abund_path}")
    linear = parse_abund_file(abund_path)
    z = element_to_z(cf.ion)
    if z not in linear:
        raise KeyError(f"Z={z} ({z_to_symbol(z)}) not in {abund_path.name}")
    factor = linear[z]
    cf.gofnt_matrix = cf.gofnt_matrix / factor
    cf.abundance_applied = False
    return cf


# --------------------------------------------------------------------------- #
# Hash / asset helpers
# --------------------------------------------------------------------------- #
def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _npz_asset_name(path: Path) -> str:
    ion, wl, ver = _filename_meta(path)
    wl_str = _FNAME_RE.match(path.name).group("wl")  # preserve exact decimals
    return f"{ion}_{wl_str}_v-{ver}.npz"


# --------------------------------------------------------------------------- #
# Dataset build
# --------------------------------------------------------------------------- #
def build_dataset(
    dat_dir: Path,
    abund_src: Path,
    out_dir: Path,
    package_data_dir: Optional[Path] = None,
) -> Path:
    """Convert every ``.dat`` and assemble the full dataset under ``out_dir``."""
    dat_dir = Path(dat_dir)
    abund_src = Path(abund_src)
    out_dir = Path(out_dir)
    data_out = out_dir / "data"
    abund_out = out_dir / "abundance"
    data_out.mkdir(parents=True, exist_ok=True)
    abund_out.mkdir(parents=True, exist_ok=True)

    catalog_rows: List[dict] = []
    manifest_files: List[dict] = []

    dat_files = sorted(dat_dir.glob("*_gofnt_v-*.dat"))
    print(f"Converting {len(dat_files)} .dat files from {dat_dir} ...")
    for dat in dat_files:
        cf = parse_gofnt_dat(dat)
        cf = strip_abundance(cf, abund_src)
        asset = _npz_asset_name(dat)
        npz_path = data_out / asset
        cf.to_npz(npz_path)
        sha = _sha256(npz_path)
        catalog_rows.append(
            {
                "ion": cf.ion,
                "wavelength": cf.wavelength,
                "chianti_version": cf.chianti_version,
                "description": cf.description,
                "ch_lower_level": cf.ch_lower_level,
                "ch_upper_level": cf.ch_upper_level,
                "f_value": cf.f_value,
                "A_value": cf.A_value,
                "Nmin": cf.Nmin,
                "Nmax": cf.Nmax,
                "units": cf.units,
                "filename": asset,
                "sha256": sha,
            }
        )
        manifest_files.append({"name": asset, "kind": "data", "sha256": sha, "size": npz_path.stat().st_size})
        print(f"  ok {asset}")

    # Copy abundance files (top-level CHIANTI .abund only).
    for ab in sorted(abund_src.glob("*.abund")):
        dest = abund_out / ab.name
        shutil.copy2(ab, dest)
        manifest_files.append(
            {"name": ab.name, "kind": "abundance", "sha256": _sha256(dest), "size": dest.stat().st_size}
        )
    print(f"Copied {len(list(abund_out.glob('*.abund')))} abundance files.")

    # Catalogue parquet.
    catalog_df = pd.DataFrame(catalog_rows)
    catalog_path = out_dir / "catalog.parquet"
    catalog_df.to_parquet(catalog_path, index=False)
    manifest_files.append(
        {"name": "catalog.parquet", "kind": "catalog", "sha256": _sha256(catalog_path), "size": catalog_path.stat().st_size}
    )
    if package_data_dir is not None:
        package_data_dir = Path(package_data_dir)
        package_data_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(catalog_path, package_data_dir / "catalog.parquet")
        print(f"Bundled catalogue into {package_data_dir / 'catalog.parquet'}")

    # Manifest.
    versions = sorted({r["chianti_version"] for r in catalog_rows})
    manifest = {
        "dataset_version": versions[0] if len(versions) == 1 else versions,
        "generated": datetime.now(timezone.utc).isoformat(),
        "n_lines": len(catalog_rows),
        "files": manifest_files,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"Wrote catalogue ({len(catalog_rows)} lines) and manifest to {out_dir}")
    return out_dir


def upload_release(out_dir: Path, repo: str, tag: str, *, title: Optional[str] = None) -> None:
    """Publish all dataset assets to a GitHub release via ``gh`` (maintainer)."""
    out_dir = Path(out_dir)
    assets = (
        list((out_dir / "data").glob("*.npz"))
        + list((out_dir / "abundance").glob("*.abund"))
        + [out_dir / "catalog.parquet", out_dir / "manifest.json"]
    )
    subprocess.run(
        ["gh", "release", "create", tag, "--repo", repo, "--title", title or tag, "--notes", "GofChianti dataset"],
        check=True,
    )
    subprocess.run(
        ["gh", "release", "upload", tag, "--repo", repo, *map(str, assets), "--clobber"],
        check=True,
    )


def main(argv: Optional[List[str]] = None) -> None:
    p = argparse.ArgumentParser(description="Build the GofChianti dataset from IDL .dat files.")
    p.add_argument("--dat-dir", required=True, type=Path, help="Directory of *_gofnt_v-*.dat files.")
    p.add_argument(
        "--abund-src",
        type=Path,
        default=Path("/usr/local/ssw/packages/chianti/dbase/abundance"),
        help="Directory of CHIANTI .abund files.",
    )
    p.add_argument("--out-dir", required=True, type=Path, help="Output dataset directory.")
    p.add_argument(
        "--package-data-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "src" / "gofchianti" / "data",
        help="Where to copy the bundled catalogue (set empty to skip).",
    )
    args = p.parse_args(argv)
    pkg_data = args.package_data_dir if str(args.package_data_dir) else None
    build_dataset(args.dat_dir, args.abund_src, args.out_dir, package_data_dir=pkg_data)


if __name__ == "__main__":
    main()

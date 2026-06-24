#!/usr/bin/env python3
"""End-user dummy-install validation for GofChianti.

Exercises the full caching system and the interpolation API, then renders a
single figure of G(T) for every available line at log10(n_e)=9.

Phases
------
1. clear cache + download_all  (populate from the dataset source)
2. download_all again          (idempotent: cache hit, nothing re-copied)
3. clear cache + download_all  ("empty and download again")
4. offline-from-cache          (no dataset dir, unreachable URL -> serve from cache)
5. interpolate every line @ n_e=1e9 -> figure (+ peak-value CSV)

Source of the dataset
---------------------
``--source local`` (default) copies from the maintainer's local ``dataset/``
directory, mirroring exactly what a *public* GitHub release would serve.
``--source release`` hits the real release URL (only works if the repo/release
is public; a private repo returns 404 for anonymous downloads).
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

REPO = Path("/home/smzergua/workshop/2026/IDL_gofnt_tool/GofChianti")
DEFAULT_DATASET = REPO / "dataset"


def _section(title: str) -> None:
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", choices=["local", "release"], default="local",
                    help="Where download_all pulls from (default: local dataset dir).")
    ap.add_argument("--dataset-dir", default=str(DEFAULT_DATASET))
    ap.add_argument("--density", type=float, default=1e9, help="n_e in cm^-3.")
    ap.add_argument(
        "--outdir", default=str(Path(__file__).resolve().parent / "outputs"))
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cache = Path(__file__).resolve().parent / "_cache"

    # Clean, isolated cache so the "empty + download again" steps are meaningful.
    os.environ["GOFCHIANTI_CACHE"] = str(cache)
    if args.source == "local":
        os.environ["GOFCHIANTI_DATASET_DIR"] = str(
            Path(args.dataset_dir).resolve())
    else:
        os.environ.pop("GOFCHIANTI_DATASET_DIR", None)

    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import gofchianti as gc

    _section("Environment")
    print("gofchianti :", gc.version())
    print("source     :", args.source)
    print("cache dir  :", gc.get_cache_dir())
    print("dataset dir:", os.environ.get(
        "GOFCHIANTI_DATASET_DIR", "(none / release)"))

    # ---------------------------------------------------------------- Phase 1
    _section("Phase 1 — clear cache + download_all")
    gc.clear_cache()
    paths = gc.download_all(verbose=False)
    print(f"downloaded/cached {len(paths)} files")
    n_npz = len(list((cache / "data").glob("*.npz")))
    n_abund = len(list((cache / "abundance").glob("*.abund")))
    print(f"cache now holds: {n_npz} .npz, {n_abund} .abund")
    assert n_npz == 56, f"expected 56 npz, got {n_npz}"

    # ---------------------------------------------------------------- Phase 2
    _section("Phase 2 — download_all again (idempotent cache hit)")
    mtimes = {p.name: p.stat().st_mtime_ns for p in (
        cache / "data").glob("*.npz")}
    paths2 = gc.download_all(verbose=False)
    changed = [n for p in (cache / "data").glob("*.npz")
               if (n := p.name) in mtimes and p.stat().st_mtime_ns != mtimes[n]]
    print(
        f"re-fetched {len(paths2)} files; files re-written this pass: {len(changed)}")

    # ---------------------------------------------------------------- Phase 3
    _section('Phase 3 — clear cache + download_all ("empty and download again")')
    gc.clear_cache()
    assert not (
        cache / "data").exists() or not any((cache / "data").glob("*.npz"))
    print("cache emptied OK")
    paths3 = gc.download_all(verbose=False)
    n_npz3 = len(list((cache / "data").glob("*.npz")))
    print(f"re-downloaded {len(paths3)} files; cache holds {n_npz3} .npz")
    assert n_npz3 == 56

    # ---------------------------------------------------------------- Phase 4
    _section("Phase 4 — offline-from-cache (no dataset dir, unreachable URL)")
    gc.set_dataset_dir(None)
    # guaranteed unreachable
    gc.set_base_url("http://127.0.0.1:9/offline-test/")
    import gofchianti.fetch as _fetch
    _fetch._manifest_cache = None  # force re-resolution
    df = gc.available_lines()
    ok = 0
    for _, r in df.iterrows():
        cf = gc.get_line(str(r["ion"]), float(r["wavelength"]))
        _ = cf.get_gofnt(
            args.density, cf.temperature[len(cf.temperature) // 2])
        ok += 1
    print(
        f"served {ok}/{len(df)} lines purely from cache (network was unreachable)")
    assert ok == len(df)

    # ---------------------------------------------------------------- Phase 5
    _section(
        f"Phase 5 — interpolate all lines @ n_e={args.density:.1e} + plot")
    elements = sorted(df["ion"].astype(str).str.split(
        "_").str[0].str.title().unique())
    cmap = plt.get_cmap("tab20")
    ecolour = {el: cmap(i % cmap.N) for i, el in enumerate(elements)}

    fig, ax = plt.subplots(figsize=(11, 7))
    rows = []
    global_peak = 0.0
    for _, r in df.iterrows():
        ion, wl = str(r["ion"]), float(r["wavelength"])
        cf = gc.get_line(ion, wl)
        T = np.asarray(cf.temperature, dtype=float)
        g = cf.get_gofnt(args.density, T)
        g = np.where(g > 0, g, np.nan)
        el = ion.split("_")[0].title()
        ax.plot(T, g, color=ecolour[el], alpha=0.75, lw=1.2)
        if np.isfinite(g).any():
            k = int(np.nanargmax(g))
            rows.append((ion, wl, float(T[k]), float(g[k])))
            global_peak = max(global_peak, float(g[k]))

    ax.set_xscale("log")
    ax.set_yscale("log")
    # Focus on the physically meaningful band: the steep low-T tails span
    # ~200 decades, so clip to 8 decades below the strongest line's peak.
    if global_peak > 0:
        ax.set_ylim(global_peak / 1e8, global_peak * 3)
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(
        r"$G(n_e, T)$  [bare, no abundance]  (erg cm$^{3}$ s$^{-1}$ sr$^{-1}$)")
    ax.set_title(f"GofChianti — all {len(df)} lines at "
                 rf"$\log_{{10}} n_e = {np.log10(args.density):.1f}$")
    handles = [plt.Line2D([], [], color=ecolour[el], label=el)
               for el in elements]
    ax.legend(handles=handles, title="Element", ncol=2, fontsize="small",
              loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()
    fig_path = outdir / "gofnt_logne9_all_lines.png"
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    print("figure ->", fig_path)

    csv_path = outdir / "gofnt_logne9_peaks.csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ion", "wavelength_A", "peak_T_K", "peak_G"])
        w.writerows(sorted(rows))
    print("peak table ->", csv_path)

    _section("ALL PHASES PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""The line catalogue: which spectral lines are available and where.

The catalogue is a small table (one row per precomputed line) shipped inside
the package as ``data/catalog.parquet`` so :func:`available_lines` works fully
offline.  When online it can be refreshed from the dataset release.
"""

from __future__ import annotations

from importlib import resources
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from . import fetch

# Columns guaranteed to be present in the catalogue.
CATALOG_COLUMNS = [
    "ion",
    "wavelength",
    "chianti_version",
    "description",
    "ch_lower_level",
    "ch_upper_level",
    "f_value",
    "A_value",
    "Nmin",
    "Nmax",
    "units",
    "filename",
    "sha256",
]

_cache: Optional[pd.DataFrame] = None


def _bundled_catalog_path() -> Optional[Path]:
    try:
        res = resources.files("gofchianti.data").joinpath(fetch.CATALOG_NAME)
        if res.is_file():
            with resources.as_file(res) as p:
                return Path(p)
    except (ModuleNotFoundError, FileNotFoundError):
        pass
    return None


def load_catalog(refresh: bool = False) -> pd.DataFrame:
    """Load the line catalogue as a :class:`pandas.DataFrame`.

    Parameters
    ----------
    refresh : bool
        Re-read from disk/remote instead of returning the in-memory copy.
    """
    global _cache
    if _cache is not None and not refresh:
        return _cache

    path: Optional[Path] = None
    # Prefer a freshly fetched copy when we can reach the dataset; otherwise the
    # bundled file keeps us working offline.
    if refresh:
        try:
            path = fetch.fetch_catalog()
        except Exception:
            path = None
    if path is None:
        path = _bundled_catalog_path()
    if path is None:
        try:
            path = fetch.fetch_catalog()
        except Exception as exc:
            raise FileNotFoundError(
                "No line catalogue available (neither bundled nor downloadable)."
            ) from exc

    df = pd.read_parquet(path)
    _cache = df
    return df


def available_lines(version: Optional[str] = None, refresh: bool = False) -> pd.DataFrame:
    """Return the table of available lines, optionally filtered by version.

    Parameters
    ----------
    version : str, optional
        Keep only rows for this CHIANTI version (e.g. ``"11.0.2"``).
    refresh : bool
        Force a reload (and remote refresh when possible).
    """
    df = load_catalog(refresh=refresh)
    if version is not None:
        df = df[df["chianti_version"].astype(str) == str(version)]
    return df.reset_index(drop=True).copy()


def resolve_line(
    ion: str,
    wavelength: float,
    version: Optional[str] = None,
    tol: float = 0.5,
) -> pd.Series:
    """Find the catalogue row for ``ion`` nearest to ``wavelength``.

    Parameters
    ----------
    ion : str
        Ion tag such as ``"Fe_12"`` (case-insensitive).
    wavelength : float
        Target wavelength in Ångström.
    version : str, optional
        Restrict to a CHIANTI version.
    tol : float
        Maximum allowed |Δλ| in Ångström for a match.
    """
    df = load_catalog()
    mask = df["ion"].astype(str).str.lower() == str(ion).lower()
    if version is not None:
        mask &= df["chianti_version"].astype(str) == str(version)
    sub = df[mask]
    if sub.empty:
        raise KeyError(f"No catalogue entries for ion {ion!r}"
                       + (f" (version {version})" if version else ""))
    dwl = (sub["wavelength"].astype(float) - float(wavelength)).abs()
    idx = dwl.idxmin()
    if dwl.loc[idx] > tol:
        nearest = sub.loc[idx, "wavelength"]
        raise KeyError(
            f"No {ion} line within {tol} Å of {wavelength} Å "
            f"(nearest is {nearest} Å)."
        )
    return sub.loc[idx]

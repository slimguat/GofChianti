"""GofChianti — lightweight Python access to precomputed CHIANTI G(n_e, T).

Typical use::

    import astropy.units as u
    import gofchianti as gc

    gc.available_lines()                    # what is available (DataFrame)
    cf = gc.get_line("Fe_12", 195.119)      # download + cache + load
    g = cf.get_gofnt(density=1e9 * u.cm**-3, # bare G at (n_e, T)
                     temperature=1.5e6 * u.K)

    # any two of (density, temperature, pressure) may be given
    g = cf.get_gofnt(pressure=1e15 * u.cm**-3 * u.K,  # reduced pressure n_e*T
                     temperature=1.5e6 * u.K)
    g = cf.get_gofnt(pressure=0.02 * u.Pa,            # thermal pressure
                     density=1e9 * u.cm**-3)

    # abundance-scaled values
    cf = gc.get_line("Fe_12", 195.119,
                     abundance="sun_photospheric_2021_asplund")
    g = cf.get_gofnt(density=1e9 * u.cm**-3,  # now multiplied by Fe/H
                     temperature=1.5e6 * u.K)
"""

from __future__ import annotations

from typing import Optional, Union
from pathlib import Path
import tomllib
from importlib.metadata import version as _get_installed_version, PackageNotFoundError

from .abundance import Abundance, available_abundances
from .catalog import available_lines
from .core import ContributionFunction
from .fetch import (
    clear_cache,
    download_all,
    get_cache_dir,
    set_base_url,
    set_cache_dir,
    set_dataset_dir,
)
from .utils import _vprint


def version(verbose=0) -> str:
    """Return the package version.

    Prefer the installed distribution (importlib.metadata). Fall back to
    reading the repository's pyproject.toml via `tomllib`.
    """
    try:
        return _get_installed_version("gofchianti")
    except Exception:
        _vprint(verbose, 2, "Package not found, falling back to pyproject.toml")
        pass

    try:
        py = Path(__file__).resolve().parents[2] / "pyproject.toml"
        if py.exists():
            with py.open("rb") as fh:
                data = tomllib.load(fh)
            ver = data.get("project", {}).get("version")
            if ver:
                return ver
            ver = data.get("tool", {}).get("poetry", {}).get("version")
            if ver:
                return ver
    except Exception:
        pass
    return "0.0.0"


def __getattr__(name: str):
    if name == "__version__":
        return version()
    raise AttributeError(name)


__all__ = [
    "get_line",
    "available_lines",
    "available_abundances",
    "Abundance",
    "ContributionFunction",
    "download_all",
    "set_cache_dir",
    "set_dataset_dir",
    "set_base_url",
    "get_cache_dir",
    "clear_cache",
]


def get_line(
    ion: str,
    wavelength: float,
    version: Optional[str] = None,
    tol: float = 0.5,
    abundance: Union[Abundance, str, Path, dict, None] = None,
):
    """Load the contribution function for a spectral line.

    The corresponding ``.npz`` file is fetched from the cache (or downloaded if
    missing) and returned as a :class:`~gofchianti.core.ContributionFunction`.

    Parameters
    ----------
    ion : str
        Ion tag, e.g. ``"Fe_12"`` (case-insensitive).
    wavelength : float
        Wavelength in Ångström; the nearest line within ``tol`` is selected.
    version : str, optional
        CHIANTI version to restrict to (e.g. ``"11.0.2"``).
    tol : float
        Maximum |Δλ| in Ångström allowed when matching the line.
    abundance : Abundance | str | Path | dict, optional
        If given, attach this abundance set and return abundance-scaled values
        from :meth:`~gofchianti.core.ContributionFunction.get_gofnt`.  By
        default the contribution function is returned *bare* (no abundance).

    Returns
    -------
    ContributionFunction
    """
    from . import catalog, fetch

    row = catalog.resolve_line(ion, wavelength, version=version, tol=tol)
    local = fetch.fetch_data_file(str(row["filename"]))
    cf = ContributionFunction.from_npz(local)
    if abundance is not None:
        cf.set_abundance(abundance, apply=True)
    return cf

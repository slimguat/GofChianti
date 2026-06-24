"""GofChianti — lightweight Python access to precomputed CHIANTI G(T, n_e).

Typical use::

    import gofchianti as gc

    gc.available_lines()                    # what is available (DataFrame)
    cf = gc.get_line("Fe_12", 195.119)      # download + cache + load
    g = cf.get_gofnt(1e9, 1.5e6)            # bare G at (n_e, T)

    # abundance-scaled values
    cf = gc.get_line("Fe_12", 195.119,
                     abundance="sun_photospheric_2021_asplund")
    g = cf.get_gofnt(1e9, 1.5e6)            # now multiplied by Fe/H
"""

from __future__ import annotations

from typing import Optional, Union
from pathlib import Path

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

__version__ = "0.1.0"

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

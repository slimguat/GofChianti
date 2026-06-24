"""End-to-end and correctness tests for GofChianti."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import gofchianti as gc
from gofchianti.abundance import parse_abund_file
from gofchianti.core import ContributionFunction

from conftest import ABUND_SRC, DATASET_DIR

ASPLUND = "sun_photospheric_2021_asplund"


# --------------------------------------------------------------------------- #
# IDL output bug fixes
# --------------------------------------------------------------------------- #
def test_asterisk_wavelength_recovered_from_filename(dataset_available):
    """h_1 1025.723 had a '*******' header; λ must be recovered from the name."""
    cf = ContributionFunction.from_npz(DATASET_DIR / "data" / "h_1_1025.723_v-11.0.2.npz")
    assert cf.ion.lower() == "h_1"
    assert cf.wavelength == pytest.approx(1025.723, abs=1e-3)


def test_abundance_file_is_basename_not_path(dataset_available):
    """The wrapped abundance path must be reduced to a basename."""
    cf = ContributionFunction.from_npz(DATASET_DIR / "data" / "fe_12_195.119_v-11.0.2.npz")
    assert cf.source_abundance_file == "sun_photospheric_2021_asplund.abund"
    assert "/" not in cf.source_abundance_file


# --------------------------------------------------------------------------- #
# Stored data is bare; abundance round-trips exactly
# --------------------------------------------------------------------------- #
def test_stored_gofnt_is_bare(dataset_available):
    cf = ContributionFunction.from_npz(DATASET_DIR / "data" / "fe_12_195.119_v-11.0.2.npz")
    assert cf.abundance_applied is False


def test_abundance_roundtrip_matches_original_idl(dataset_available):
    """bare_G × (Fe/H) must reproduce the original abundance-multiplied .dat."""
    from maintainers.convert_dat_to_npz import parse_gofnt_dat  # noqa: WPS433

    dat = Path(__file__).resolve().parent.parent.parent / "gofnt" / "fe_12_195.119_gofnt_v-11.0.2.dat"
    raw = parse_gofnt_dat(dat)  # original, abundance applied
    assert raw.abundance_applied is True

    bare = ContributionFunction.from_npz(DATASET_DIR / "data" / "fe_12_195.119_v-11.0.2.npz")
    factor = parse_abund_file(ABUND_SRC / (ASPLUND + ".abund"))[26]  # Fe

    np.testing.assert_allclose(bare.gofnt_matrix * factor, raw.gofnt_matrix, rtol=1e-10, atol=0.0)


# --------------------------------------------------------------------------- #
# Public API: catalogue + get_line + interpolation
# --------------------------------------------------------------------------- #
def test_available_lines_dataframe():
    df = gc.available_lines()
    assert len(df) > 0
    for col in ("ion", "wavelength", "chianti_version", "f_value", "A_value", "filename"):
        assert col in df.columns
    assert (gc.available_lines(version="11.0.2")["chianti_version"] == "11.0.2").all()


def test_get_line_bare(dataset_available):
    cf = gc.get_line("Fe_12", 195.119)
    assert isinstance(cf, ContributionFunction)
    assert cf.get_abundance is False
    g = cf.get_gofnt(1e9, 1.5e6)
    assert np.ndim(g) == 0 or g.shape == (1,)
    assert np.isfinite(g).all() and (np.asarray(g) > 0).all()


def test_get_line_with_abundance_scales_by_factor(dataset_available):
    cf_bare = gc.get_line("Fe_12", 195.119)
    cf_ab = gc.get_line("Fe_12", 195.119, abundance=ASPLUND)
    factor = parse_abund_file(ABUND_SRC / (ASPLUND + ".abund"))[26]

    n, t = 1e9, 1.5e6
    g_bare = np.asarray(cf_bare.get_gofnt(n, t))
    g_ab = np.asarray(cf_ab.get_gofnt(n, t))
    np.testing.assert_allclose(g_ab, g_bare * factor, rtol=1e-10)


def test_interpolation_recovers_grid_values(dataset_available):
    cf = ContributionFunction.from_npz(DATASET_DIR / "data" / "fe_12_195.119_v-11.0.2.npz")
    di, ti = 5, 10
    n = cf.densities[di]
    t = cf.temperature[ti]
    g = np.asarray(cf.get_gofnt(n, t)).ravel()[0]
    assert g == pytest.approx(cf.gofnt_matrix[di, ti], rel=1e-6)


def test_set_abundance_toggle(dataset_available):
    cf = gc.get_line("Fe_12", 195.119)
    n, t = 1e9, 1.5e6
    bare = np.asarray(cf.get_gofnt(n, t))
    factor = parse_abund_file(ABUND_SRC / (ASPLUND + ".abund"))[26]  # Fe/H < 1
    cf.set_abundance(ASPLUND)
    assert cf.get_abundance is True
    scaled = np.asarray(cf.get_gofnt(n, t))
    np.testing.assert_allclose(scaled, bare * factor, rtol=1e-10)
    assert not np.allclose(scaled, bare, rtol=1e-6, atol=0.0)
    cf.set_abundance(None)
    assert cf.get_abundance is False
    np.testing.assert_allclose(np.asarray(cf.get_gofnt(n, t)), bare, rtol=1e-12)


# --------------------------------------------------------------------------- #
# Abundance helper
# --------------------------------------------------------------------------- #
def test_inline_abundance_no_loading():
    ab = gc.Abundance.from_dex({"Fe": 8.0, "H": 12.0})
    assert ab.is_inline
    assert ab.linear("Fe") == pytest.approx(10.0 ** (8.0 - 12.0))
    assert ab["H"] == pytest.approx(1.0)


def test_available_abundances_lists_sets(dataset_available):
    names = gc.available_abundances()
    assert ASPLUND in names


def test_download_all_populates_cache(dataset_available):
    paths = gc.download_all(verbose=False)
    assert len(paths) > 0
    assert all(Path(p).exists() for p in paths)

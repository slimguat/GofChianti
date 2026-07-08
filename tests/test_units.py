"""Tests for the astropy unit-aware ``get_gofnt`` API and ``gofchianti.units``.

The helper-level tests run without the dataset; the interpolation-consistency
tests use the locally built dataset via the ``dataset_available`` fixture.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np
import pytest
from astropy.constants import k_B

import gofchianti as gc
from gofchianti import units as gu
from gofchianti.core import ContributionFunction

from conftest import DATASET_DIR

FE12 = "fe_12_195.119_v-11.0.2.npz"
GOFNT_UNIT = u.erg * u.cm**3 / (u.s * u.sr)


# --------------------------------------------------------------------------- #
# Unit-string parsing
# --------------------------------------------------------------------------- #
def test_parse_gofnt_units_caret_string():
    unit = gu.parse_gofnt_units("erg cm^3 s^-1 sr^-1")
    assert unit.is_equivalent(GOFNT_UNIT)


def test_parse_gofnt_units_empty_is_dimensionless():
    assert gu.parse_gofnt_units("") == u.dimensionless_unscaled
    assert gu.parse_gofnt_units(None) == u.dimensionless_unscaled


# --------------------------------------------------------------------------- #
# Quantity assertions
# --------------------------------------------------------------------------- #
def test_require_quantity_rejects_bare_number():
    with pytest.raises(TypeError):
        gu.require_quantity(1e9, "density")


def test_to_number_density_rejects_wrong_physical_type():
    with pytest.raises(u.UnitsError):
        gu.to_number_density(1.5e6 * u.K)


def test_to_temperature_accepts_celsius_equivalency():
    t = gu.to_temperature(0 * u.deg_C)
    assert t.unit == u.K
    assert t.value == pytest.approx(273.15)


def test_to_number_density_converts_m3_to_cm3():
    ne = gu.to_number_density(1e15 * u.m**-3)
    assert ne.unit == u.cm**-3
    assert ne.value == pytest.approx(1e9)


# --------------------------------------------------------------------------- #
# Pressure -> reduced pressure
# --------------------------------------------------------------------------- #
def test_reduced_pressure_passthrough():
    reduced = gu.pressure_to_reduced(1e15 * u.cm**-3 * u.K)
    assert reduced.unit.is_equivalent(u.cm**-3 * u.K)
    assert reduced.to_value(u.cm**-3 * u.K) == pytest.approx(1e15)


def test_thermal_pressure_divided_by_kb():
    # P = n_e k_B T  ->  reduced = P / k_B = n_e T
    ne, t = 1e9 * u.cm**-3, 1.5e6 * u.K
    p = (ne * k_B * t).to(u.Pa)
    reduced = gu.pressure_to_reduced(p)
    assert reduced.to_value(u.cm**-3 * u.K) == pytest.approx(1.5e15, rel=1e-10)


def test_pressure_wrong_units_raises():
    with pytest.raises(u.UnitsError):
        gu.pressure_to_reduced(5 * u.cm)


# --------------------------------------------------------------------------- #
# resolve_ne_temperature
# --------------------------------------------------------------------------- #
def test_resolve_requires_exactly_two():
    with pytest.raises(ValueError):
        gu.resolve_ne_temperature(density=1e9 * u.cm**-3)
    with pytest.raises(ValueError):
        gu.resolve_ne_temperature(
            density=1e9 * u.cm**-3,
            temperature=1e6 * u.K,
            pressure=1e15 * u.cm**-3 * u.K,
        )


def test_resolve_pressure_temperature_gives_density():
    ne, t = gu.resolve_ne_temperature(
        pressure=1.5e15 * u.cm**-3 * u.K, temperature=1.5e6 * u.K
    )
    assert ne.to_value(u.cm**-3) == pytest.approx(1e9)
    assert t.to_value(u.K) == pytest.approx(1.5e6)


def test_resolve_pressure_density_gives_temperature():
    ne, t = gu.resolve_ne_temperature(
        pressure=1.5e15 * u.cm**-3 * u.K, density=1e9 * u.cm**-3
    )
    assert ne.to_value(u.cm**-3) == pytest.approx(1e9)
    assert t.to_value(u.K) == pytest.approx(1.5e6)


# --------------------------------------------------------------------------- #
# get_gofnt: unit-aware behaviour + input-combination consistency
# --------------------------------------------------------------------------- #
def _cf():
    return ContributionFunction.from_npz(DATASET_DIR / "data" / FE12)


def test_get_gofnt_returns_quantity(dataset_available):
    cf = _cf()
    g = cf.get_gofnt(density=1e9 * u.cm**-3, temperature=1.5e6 * u.K)
    assert isinstance(g, u.Quantity)
    assert g.unit.is_equivalent(GOFNT_UNIT)
    assert (g.value > 0).all()


def test_get_gofnt_rejects_bare_numbers(dataset_available):
    cf = _cf()
    with pytest.raises(TypeError):
        cf.get_gofnt(1e9, 1.5e6)


def test_get_gofnt_rejects_wrong_units(dataset_available):
    cf = _cf()
    with pytest.raises(u.UnitsError):
        cf.get_gofnt(density=1e9 * u.K, temperature=1.5e6 * u.K)


def test_get_gofnt_reduced_pressure_matches_density(dataset_available):
    cf = _cf()
    ne, t = 1e9 * u.cm**-3, 1.5e6 * u.K
    g_nt = cf.get_gofnt(density=ne, temperature=t)
    g_pt = cf.get_gofnt(pressure=ne * t, temperature=t)
    np.testing.assert_allclose(g_pt.value, g_nt.value, rtol=1e-12)


def test_get_gofnt_thermal_pressure_matches_density(dataset_available):
    cf = _cf()
    ne, t = 1e9 * u.cm**-3, 1.5e6 * u.K
    p = (ne * k_B * t).to(u.Pa)
    g_nt = cf.get_gofnt(density=ne, temperature=t)
    g_pt = cf.get_gofnt(pressure=p, temperature=t)
    np.testing.assert_allclose(g_pt.value, g_nt.value, rtol=1e-10)


def test_get_gofnt_pressure_density_matches_density_temperature(dataset_available):
    cf = _cf()
    ne, t = 1e9 * u.cm**-3, 1.5e6 * u.K
    g_nt = cf.get_gofnt(density=ne, temperature=t)
    g_pn = cf.get_gofnt(pressure=ne * t, density=ne)
    np.testing.assert_allclose(g_pn.value, g_nt.value, rtol=1e-12)


def test_get_gofnt_array_inputs(dataset_available):
    cf = _cf()
    ne = np.array([1e8, 1e9, 1e10]) * u.cm**-3
    t = np.array([1e6, 1.5e6, 2e6]) * u.K
    g = cf.get_gofnt(density=ne, temperature=t)
    assert g.shape == (3,)
    assert (g.value > 0).all()


def test_get_line_pressure_api(dataset_available):
    cf = gc.get_line("Fe_12", 195.119)
    g = cf.get_gofnt(pressure=1.5e15 * u.cm**-3 * u.K, temperature=1.5e6 * u.K)
    assert isinstance(g, u.Quantity)
    assert g.unit.is_equivalent(GOFNT_UNIT)

"""Astropy-based unit handling for the contribution function API.

The public :meth:`~gofchianti.core.ContributionFunction.get_gofnt` accepts the
plasma state as :class:`astropy.units.Quantity` objects and returns ``G`` as a
``Quantity``.  A user may specify *any two* of

* electron density ``n_e``   (convertible to ``cm^-3``),
* temperature ``T``          (convertible to ``K``),
* pressure ``P``             (see below),

and this module resolves them to the ``(n_e, T)`` grid the stored data lives on.

Two flavours of "pressure" are accepted and told apart purely by their units:

* **True (thermal) pressure** — anything convertible to :data:`astropy.units.Pa`
  (i.e. energy per volume, ``P = n_e k_B T``).  Electron density is recovered as
  ``n_e = P / (k_B T)``.
* **Reduced pressure** ``P_e = n_e T`` — anything convertible to ``cm^-3 K``
  (no ``k_B``; this is the quantity CHIANTI's IDL ``gofnt`` calls "pressure").
  Electron density is recovered as ``n_e = P_e / T``.

The two are never dimensionally ambiguous: a true pressure has dimensions of
energy/volume while a reduced pressure has dimensions of number-density ×
temperature, and the ratio between them is exactly ``k_B``.
"""

from __future__ import annotations

from typing import Optional, Tuple

import astropy.units as u
from astropy.constants import k_B

__all__ = [
    "NE_UNIT",
    "T_UNIT",
    "REDUCED_PRESSURE_UNIT",
    "require_quantity",
    "to_number_density",
    "to_temperature",
    "pressure_to_reduced",
    "resolve_ne_temperature",
    "parse_gofnt_units",
]

#: Canonical working units the stored grid is expressed in.
NE_UNIT = u.cm**-3
T_UNIT = u.K
#: Reduced pressure ``P_e = n_e * T`` (CHIANTI-style, no ``k_B``).
REDUCED_PRESSURE_UNIT = u.cm**-3 * u.K


def require_quantity(value: object, name: str) -> u.Quantity:
    """Return ``value`` as a :class:`~astropy.units.Quantity` or raise.

    A bare number (``float``/``int``/ndarray) is rejected: the whole point of
    the unit-aware API is that the caller states the physical units explicitly.
    """
    if not isinstance(value, u.Quantity):
        raise TypeError(
            f"{name} must be an astropy Quantity with explicit units "
            f"(e.g. 1e9 * u.cm**-3), got {type(value).__name__}: {value!r}."
        )
    return value


def to_number_density(value: object, name: str = "density") -> u.Quantity:
    """Validate and convert ``value`` to an electron density in ``cm^-3``."""
    q = require_quantity(value, name)
    if not q.unit.is_equivalent(NE_UNIT):
        raise u.UnitsError(
            f"{name} must be a number density convertible to cm^-3, "
            f"got units {q.unit!s} (physical type: {q.unit.physical_type})."
        )
    return q.to(NE_UNIT)


def to_temperature(value: object, name: str = "temperature") -> u.Quantity:
    """Validate and convert ``value`` to a temperature in ``K``."""
    q = require_quantity(value, name)
    try:
        return q.to(T_UNIT, equivalencies=u.temperature())
    except u.UnitConversionError as exc:
        raise u.UnitsError(
            f"{name} must be a temperature convertible to K, "
            f"got units {q.unit!s} (physical type: {q.unit.physical_type})."
        ) from exc


def pressure_to_reduced(value: object, name: str = "pressure") -> u.Quantity:
    """Return the reduced pressure ``n_e * T`` (in ``cm^-3 K``) for ``value``.

    Accepts either a true thermal pressure (convertible to ``Pa``; divided by
    ``k_B``) or an already-reduced pressure (convertible to ``cm^-3 K``).
    """
    q = require_quantity(value, name)
    if q.unit.is_equivalent(u.Pa):
        return (q / k_B).to(REDUCED_PRESSURE_UNIT)
    if q.unit.is_equivalent(REDUCED_PRESSURE_UNIT):
        return q.to(REDUCED_PRESSURE_UNIT)
    raise u.UnitsError(
        f"{name} must be either a thermal pressure (convertible to Pa, i.e. "
        f"P = n_e k_B T) or a reduced pressure n_e*T (convertible to cm^-3 K); "
        f"got units {q.unit!s} (physical type: {q.unit.physical_type})."
    )


def resolve_ne_temperature(
    density: object = None,
    temperature: object = None,
    pressure: object = None,
) -> Tuple[u.Quantity, u.Quantity]:
    """Resolve any two of ``(density, temperature, pressure)`` to ``(n_e, T)``.

    Returns a ``(n_e, T)`` tuple of :class:`~astropy.units.Quantity` in
    ``cm^-3`` and ``K`` respectively.  Exactly two inputs must be supplied.
    """
    provided = [x is not None for x in (density, temperature, pressure)]
    if sum(provided) != 2:
        raise ValueError(
            "Provide exactly two of (density, temperature, pressure). "
            f"Got density={density is not None}, "
            f"temperature={temperature is not None}, "
            f"pressure={pressure is not None}."
        )

    ne: Optional[u.Quantity] = (
        to_number_density(density) if density is not None else None
    )
    temp: Optional[u.Quantity] = (
        to_temperature(temperature) if temperature is not None else None
    )

    if pressure is not None:
        reduced = pressure_to_reduced(pressure)
        if temp is not None:
            ne = (reduced / temp).to(NE_UNIT)
        else:  # density is the other provided quantity
            assert ne is not None
            temp = (reduced / ne).to(T_UNIT)

    assert ne is not None and temp is not None
    return ne, temp


def parse_gofnt_units(units_str: Optional[str]) -> u.UnitBase:
    """Parse a stored G(T,n) units string into an astropy unit.

    The dataset stores units like ``"erg cm^3 s^-1 sr^-1"``.  Astropy's generic
    parser does not accept the ``^`` caret, so it is stripped first
    (``cm^3`` -> ``cm3``, ``s^-1`` -> ``s-1``).  Falls back to a dimensionless
    unit if the string is empty or cannot be parsed.
    """
    if not units_str or not units_str.strip():
        return u.dimensionless_unscaled
    cleaned = units_str.replace("^", "").strip()
    try:
        return u.Unit(cleaned)
    except Exception:
        try:
            return u.Unit(cleaned, parse_strict="silent")
        except Exception:
            return u.dimensionless_unscaled

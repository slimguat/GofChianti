"""The :class:`ContributionFunction` object and ``.npz`` (de)serialisation.

A :class:`ContributionFunction` stores the *bare* contribution function
``G(T, n_e)`` (no elemental abundance applied) on a density/temperature grid,
together with the relevant CHIANTI metadata.  ``G`` is multiplied by an
elemental abundance only on request (see :meth:`ContributionFunction.get_gofnt`).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import RegularGridInterpolator

from .abundance import Abundance, coerce_abundance, element_to_z, z_to_symbol

# Metadata keys serialised alongside the arrays in the ``.npz`` file.
_META_KEYS = (
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
    "source_abundance_file",
    "abundance_applied",
    "units",
    "generation_date",
)


@dataclass
class ContributionFunction:
    """Bare CHIANTI contribution function ``G(T, n_e)`` with metadata.

    The stored :attr:`gofnt_matrix` does **not** include any elemental
    abundance.  Use :meth:`set_abundance` / :meth:`get_gofnt` to obtain
    abundance-scaled values.
    """

    ion: str
    wavelength: float
    chianti_version: str
    description: str
    ch_lower_level: int
    ch_upper_level: int
    f_value: float
    A_value: float
    Nmin: float
    Nmax: float
    source_abundance_file: str
    abundance_applied: bool
    units: str
    generation_date: str
    densities: NDArray
    temperature: NDArray
    gofnt_matrix: NDArray
    filename: Optional[str] = None

    # Runtime-only state (not serialised).
    get_abundance: bool = field(default=False, compare=False)
    abundance: Optional[Abundance] = field(default=None, repr=False, compare=False)
    interpolator: Optional[RegularGridInterpolator] = field(
        default=None, repr=False, compare=False
    )

    # ------------------------------------------------------------------ #
    # Identity helpers
    # ------------------------------------------------------------------ #
    @property
    def element(self) -> str:
        """Element symbol of the ion (e.g. ``"Fe"`` for ``"fe_12"``)."""
        return z_to_symbol(self.Z)

    @property
    def Z(self) -> int:
        """Atomic number of the ion."""
        return element_to_z(self.ion)

    def __repr__(self) -> str:
        YELLOW, RESET = "\033[93m", "\033[0m"
        ab = "+abund" if self.get_abundance else "bare"
        return (
            f"<ContributionFunction {YELLOW}{self.ion} {self.wavelength:.3f} Å{RESET} | "
            f"{len(self.densities)}n × {len(self.temperature)}T | {ab}>"
        )

    def describe(self) -> str:
        """Human-readable multi-line summary."""
        lines = [
            "Contribution Function Details",
            f"Ion: {self.ion}  (Z={self.Z}, {self.element})",
            f"Wavelength: {self.wavelength:.3f} Å",
            f"CHIANTI version: {self.chianti_version}",
            f"Description: {self.description}",
            f"Transition levels: {self.ch_lower_level} → {self.ch_upper_level}",
            f"Oscillator strength (f): {self.f_value:.3e}",
            f"Einstein A value (s^-1): {self.A_value:.3e}",
            f"Density range: {self.Nmin:.3e} – {self.Nmax:.3e} cm^-3",
            f"Stored G includes abundance: {self.abundance_applied}",
            f"Source abundance file: {self.source_abundance_file}",
            f"Units: {self.units}",
            f"Generated on: {self.generation_date}",
            f"Densities: {len(self.densities)} pts "
            f"({self.densities[0]:.2e} → {self.densities[-1]:.2e})",
            f"Temperatures: {len(self.temperature)} pts "
            f"({self.temperature[0]:.2e} → {self.temperature[-1]:.2e})",
            f"G matrix shape: {self.gofnt_matrix.shape}",
        ]
        if self.filename:
            lines.append(f"Source file: {self.filename}")
        return "\n".join(lines)

    # ------------------------------------------------------------------ #
    # Abundance handling
    # ------------------------------------------------------------------ #
    def set_abundance(
        self,
        abundance: Union[Abundance, str, Path, dict, None],
        *,
        apply: bool = True,
    ) -> "ContributionFunction":
        """Attach an abundance set and toggle whether it is applied.

        ``abundance`` may be an :class:`~gofchianti.abundance.Abundance`, a
        reference name, a path to a ``.abund`` file, a ``{Z: linear}`` mapping,
        or ``None`` (which disables abundance scaling).
        """
        self.abundance = coerce_abundance(abundance)
        self.get_abundance = bool(apply and self.abundance is not None)
        return self

    def _abundance_factor(self) -> float:
        if self.abundance is None:
            raise ValueError(
                "Abundance scaling requested but no abundance set has been "
                "attached. Call set_abundance(...) or pass abundance= to "
                "get_line()."
            )
        return float(self.abundance.linear(self.ion))

    # ------------------------------------------------------------------ #
    # Interpolation
    # ------------------------------------------------------------------ #
    def get_gofnt(
        self,
        density: Union[float, NDArray],
        temperature: Union[float, NDArray],
        include_abundance: Optional[bool] = None,
        recompute_interpolation: bool = False,
    ) -> NDArray:
        """Interpolate ``G(T, n_e)`` at the requested density/temperature.

        Parameters
        ----------
        density, temperature : float or ndarray
            Query points (cm^-3 and K).  Either may be scalar; otherwise they
            must broadcast to a common shape.
        include_abundance : bool, optional
            Whether to multiply by the attached elemental abundance.  Defaults
            to the instance flag :attr:`get_abundance`.
        recompute_interpolation : bool
            Force rebuilding of the cached interpolator.
        """
        if include_abundance is None:
            include_abundance = self.get_abundance

        density = np.atleast_1d(np.asarray(density, dtype=float))
        temperature = np.atleast_1d(np.asarray(temperature, dtype=float))

        if (
            density.shape != temperature.shape
            and density.size > 1
            and temperature.size > 1
        ):
            raise ValueError(
                "Density and temperature arrays must share a shape or one must "
                f"be scalar-like. Got {density.shape} and {temperature.shape}."
            )
        if density.shape != temperature.shape:
            density, temperature = np.broadcast_arrays(density, temperature)

        original_shape = density.shape
        log_dens_q = np.log10(density.ravel())
        log_temp_q = np.log10(temperature.ravel())

        if self.interpolator is None or recompute_interpolation:
            log_dens_grid = np.log10(np.asarray(self.densities, dtype=float))
            log_temp_grid = np.log10(np.asarray(self.temperature, dtype=float))
            with np.errstate(divide="ignore"):
                log_g_grid = np.log10(np.asarray(self.gofnt_matrix, dtype=float))
            expected = (log_dens_grid.size, log_temp_grid.size)
            if log_g_grid.shape != expected:
                raise ValueError(
                    f"gofnt_matrix shape {log_g_grid.shape} does not match "
                    f"(n_density, n_temperature)={expected}."
                )
            self.interpolator = RegularGridInterpolator(
                (log_dens_grid, log_temp_grid),
                log_g_grid,
                bounds_error=False,
                fill_value=-np.inf,  # -> 0 after the 10** below
            )

        points = np.column_stack([log_dens_q, log_temp_q])
        g = np.power(10.0, self.interpolator(points)).reshape(original_shape)

        if include_abundance:
            g = g * self._abundance_factor()
        return g

    # ------------------------------------------------------------------ #
    # Serialisation
    # ------------------------------------------------------------------ #
    def to_npz(self, path: str | Path) -> Path:
        """Write this contribution function to a compressed ``.npz`` file."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        meta = {k: getattr(self, k) for k in _META_KEYS}
        np.savez_compressed(
            path,
            densities=np.asarray(self.densities, dtype=float),
            temperature=np.asarray(self.temperature, dtype=float),
            gofnt=np.asarray(self.gofnt_matrix, dtype=float),
            meta=np.asarray(json.dumps(meta)),
        )
        # numpy appends .npz if missing
        return path if path.suffix == ".npz" else path.with_suffix(".npz")

    @classmethod
    def from_npz(cls, path: str | Path) -> "ContributionFunction":
        """Load a contribution function from a ``.npz`` file."""
        path = Path(path)
        with np.load(path, allow_pickle=False) as data:
            meta = json.loads(str(data["meta"]))
            densities = np.asarray(data["densities"], dtype=float)
            temperature = np.asarray(data["temperature"], dtype=float)
            gofnt = np.asarray(data["gofnt"], dtype=float)
        return cls(
            densities=densities,
            temperature=temperature,
            gofnt_matrix=gofnt,
            filename=str(path),
            **meta,
        )

    def copy(self) -> "ContributionFunction":
        """Return a copy sharing the (immutable) cached interpolator."""
        new = ContributionFunction(
            ion=self.ion,
            wavelength=self.wavelength,
            chianti_version=self.chianti_version,
            description=self.description,
            ch_lower_level=self.ch_lower_level,
            ch_upper_level=self.ch_upper_level,
            f_value=self.f_value,
            A_value=self.A_value,
            Nmin=self.Nmin,
            Nmax=self.Nmax,
            source_abundance_file=self.source_abundance_file,
            abundance_applied=self.abundance_applied,
            units=self.units,
            generation_date=self.generation_date,
            densities=np.copy(self.densities),
            temperature=np.copy(self.temperature),
            gofnt_matrix=np.copy(self.gofnt_matrix),
            filename=self.filename,
        )
        new.get_abundance = self.get_abundance
        new.abundance = self.abundance
        new.interpolator = self.interpolator
        return new

"""CHIANTI abundance files and the :class:`Abundance` helper.

``G(T, n_e)`` is stored *without* any elemental abundance applied.  When a user
asks for abundance-scaled values, the contribution function is multiplied by the
linear abundance ``N_X / N_H`` of the relevant element, Available abundances are 
read from a available CHIANTI ``.abund`` file.

CHIANTI ``.abund`` format
-------------------------
Each data row is ``Z  A(X)  Symbol`` where ``A(X) = 12 + log10(N_X / N_H)`` (so
hydrogen is ``12.00``).  Lines starting with ``%`` are comments and the data
block is terminated by a line containing ``-1``.  The linear abundance is
therefore ``10 ** (A(X) - 12)``.

A :class:`Abundance` can be one of:

* a *reference* identified by ``name`` (e.g. ``"sun_photospheric_2021_asplund"``)
  whose ``.abund`` file is fetched from the cache/remote on first use;
* a *local file* given by ``path``;
* a *self-contained* set given by ``values`` (linear ``N_X / N_H`` keyed by
  atomic number ``Z``) which needs no loading at all.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Union

# Atomic symbols indexed by Z (1..30 covers every element in the CHIANTI
# database).  Index 0 is a placeholder so that ``_SYMBOLS[Z]`` is the symbol.
_SYMBOLS: List[str] = [
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
]
_SYMBOL_TO_Z: Dict[str, int] = {
    sym.lower(): z for z, sym in enumerate(_SYMBOLS) if sym}

ElementLike = Union[int, str]


def element_to_z(element: ElementLike) -> int:
    """Resolve an element symbol, ion tag, or atomic number to ``Z``.

    Accepts ``26``, ``"Fe"``, ``"fe"`` or an ion tag such as ``"fe_12"``.
    """
    if isinstance(element, (int,)) and not isinstance(element, bool):
        return int(element)
    if isinstance(element, str):
        token = element.strip().split("_")[0].lower()
        if token in _SYMBOL_TO_Z:
            return _SYMBOL_TO_Z[token]
        if token.isdigit():
            return int(token)
    raise KeyError(f"Unknown element: {element!r}")


def z_to_symbol(z: int) -> str:
    """Return the atomic symbol for atomic number ``z``."""
    if 1 <= z < len(_SYMBOLS):
        return _SYMBOLS[z]
    raise KeyError(f"No symbol for Z={z}")


def parse_abund_file(path: str | Path) -> Dict[int, float]:
    """Parse a CHIANTI ``.abund`` file into ``{Z: linear N_X/N_H}``."""
    path = Path(path)
    values: Dict[int, float] = {}
    with path.open("r") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("%"):
                continue
            parts = line.split()
            try:
                z = int(parts[0])
            except (ValueError, IndexError):
                continue
            if z < 0:  # terminator (-1)
                break
            if len(parts) < 2:
                continue
            a_x = float(parts[1])  # 12 + log10(N_X/N_H)
            values[z] = 10.0 ** (a_x - 12.0)
    if not values:
        raise ValueError(f"No abundance values parsed from {path}")
    return values


@dataclass
class Abundance:
    """An elemental abundance set used to scale ``G(T, n_e)``.

    Parameters
    ----------
    name : str, optional
        Reference name of a CHIANTI abundance file *without* extension, e.g.
        ``"sun_photospheric_2021_asplund"``.  Resolved via the cache/remote.
    path : str or pathlib.Path, optional
        A local ``.abund`` file to load.
    values : dict, optional
        Self-contained mapping ``{Z: linear N_X/N_H}`` (or use
        :meth:`from_dex`).  When provided, no file is ever loaded.
    """

    name: Optional[str] = None
    path: Optional[Union[str, Path]] = None
    values: Optional[Dict[int, float]] = None
    _resolved: bool = field(default=False, repr=False)

    def __post_init__(self) -> None:
        if self.values is not None:
            self.values = {int(z): float(v) for z, v in self.values.items()}
            self._resolved = True
        if self.path is not None:
            self.path = Path(self.path)

    # -- constructors ------------------------------------------------------- #
    @classmethod
    def from_file(cls, path: str | Path) -> "Abundance":
        """Build an abundance set from a local ``.abund`` file."""
        return cls(path=Path(path))

    @classmethod
    def from_name(cls, name: str) -> "Abundance":
        """Build a reference abundance resolved from the dataset by name."""
        return cls(name=name)

    @classmethod
    def from_values(cls, values: Dict[int, float]) -> "Abundance":
        """Build a self-contained abundance from linear ``{Z: N_X/N_H}``."""
        return cls(values=dict(values))

    @classmethod
    def from_dex(cls, values: Dict[ElementLike, float]) -> "Abundance":
        """Build a self-contained abundance from ``A(X) = 12 + log10(N_X/N_H)``.

        Keys may be atomic numbers or element symbols.
        """
        linear = {element_to_z(k): 10.0 ** (v - 12.0)
                  for k, v in values.items()}
        return cls(values=linear)

    # -- resolution --------------------------------------------------------- #
    @property
    def is_inline(self) -> bool:
        """True if the abundance values are self-contained (no loading needed)."""
        return self.values is not None and self.name is None and self.path is None

    def resolve(self) -> "Abundance":
        """Ensure ``values`` are populated, loading/fetching if required."""
        if self._resolved and self.values is not None:
            return self
        if self.path is not None:
            self.values = parse_abund_file(self.path)
        elif self.name is not None:
            from . import fetch

            filename = self.name if self.name.endswith(
                ".abund") else f"{self.name}.abund"
            local = fetch.fetch_abundance_file(filename)
            self.values = parse_abund_file(local)
        else:
            raise ValueError(
                "Abundance has no name, path, or inline values to resolve.")
        self._resolved = True
        return self

    # -- access ------------------------------------------------------------- #
    def linear(self, element: ElementLike) -> float:
        """Return the linear abundance ``N_X / N_H`` for ``element``."""
        self.resolve()
        z = element_to_z(element)
        assert self.values is not None
        if z not in self.values:
            raise KeyError(
                f"Element Z={z} ({z_to_symbol(z) if z < len(_SYMBOLS) else '?'}) "
                f"not present in abundance set {self.label!r}."
            )
        return self.values[z]

    def __getitem__(self, element: ElementLike) -> float:
        return self.linear(element)

    @property
    def label(self) -> str:
        if self.name:
            return self.name
        if self.path:
            return Path(self.path).stem
        return "inline"

    def __repr__(self) -> str:
        kind = "inline" if self.is_inline else (
            "file" if self.path else "name")
        return f"<Abundance {self.label!r} ({kind})>"


def coerce_abundance(spec: Union["Abundance", str, Path, Dict[int, float], None]) -> Optional["Abundance"]:
    """Turn a user-supplied abundance specification into an :class:`Abundance`.

    Accepts an :class:`Abundance`, a reference name, a path to a ``.abund``
    file, a ``{Z: linear}`` mapping, or ``None``.
    """
    if spec is None:
        return None
    if isinstance(spec, Abundance):
        return spec
    if isinstance(spec, dict):
        return Abundance.from_values(spec)
    if isinstance(spec, (str, Path)):
        p = Path(spec)
        if p.suffix == ".abund" or p.exists():
            return Abundance.from_file(p)
        return Abundance.from_name(str(spec))
    raise TypeError(f"Cannot interpret abundance specification: {spec!r}")


def available_abundances() -> List[str]:
    """List abundance set names available in the dataset manifest."""
    from . import fetch

    manifest = fetch.load_manifest()
    names = [
        Path(name).stem
        for name, entry in manifest.items()
        if entry.get("kind") == "abundance"
    ]
    return sorted(names)

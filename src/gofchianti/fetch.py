"""Cache directory management and on-demand dataset downloading.

The GofChianti dataset (precomputed ``G(n_e, T)`` ``.npz`` files, CHIANTI
``.abund`` abundance files and the line catalogue) is hosted as flat release
assets on GitHub.  Files are downloaded the
first time they are needed and cached on disk so the tool works offline
afterwards.

Cache location resolution order
-------------------------------
1. An explicit directory set via :func:`set_cache_dir`.
2. The ``GOFCHIANTI_CACHE`` environment variable.
3. The OS-standard per-user cache directory (via :mod:`platformdirs`),
   e.g. ``~/.cache/gofchianti`` on Linux.

For development / fully offline use a local dataset directory may be provided
through the ``GOFCHIANTI_DATASET_DIR`` environment variable (or
:func:`set_dataset_dir`); files are then copied from there instead of being
downloaded.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import time
from pathlib import Path
from typing import Dict, List, Optional

import platformdirs

__all__ = [
    "get_cache_dir",
    "set_cache_dir",
    "set_dataset_dir",
    "set_base_url",
    "clear_cache",
    "fetch_file",
    "fetch_data_file",
    "fetch_abundance_file",
    "fetch_catalog",
    "load_manifest",
    "download_all",
]

# Default remote location (GitHub release "latest" assets).  Overridable with
# the GOFCHIANTI_BASE_URL environment variable or :func:`set_base_url` so the
#TODO: In the future the dataset can be served from a lab URL without code changes. So set th needed tests for this chang by then.
_DEFAULT_BASE_URL = "https://github.com/slimguat/GofChianti/releases/latest/download/"

MANIFEST_NAME = "manifest.json"
CATALOG_NAME = "catalog.parquet"

# Subdirectory used inside the cache to keep each kind of file tidy.  The remote
# assets themselves are flat (release assets cannot contain ``/``).
_KIND_SUBDIR = {
    "data": "data",
    "abundance": "abundance",
    "catalog": "",
    "manifest": "",
}

# Module-level overrides (set via the public setters).
_cache_dir_override: Optional[Path] = None
_dataset_dir_override: Optional[Path] = None
_base_url_override: Optional[str] = None

# In-process manifest cache (logical name -> entry dict) with a short TTL so we
# honour "check the remote for updates" without issuing a request per file.
_manifest_cache: Optional[Dict[str, dict]] = None
_manifest_fetched_at: float = 0.0
_MANIFEST_TTL_SECONDS = 60.0


# --------------------------------------------------------------------------- #
# Configuration setters
# --------------------------------------------------------------------------- #
def set_cache_dir(path: str | os.PathLike) -> Path:
    """Override the cache directory for the current session."""
    global _cache_dir_override
    _cache_dir_override = Path(path).expanduser().resolve()
    _cache_dir_override.mkdir(parents=True, exist_ok=True)
    return _cache_dir_override


def set_dataset_dir(path: str | os.PathLike | None) -> Optional[Path]:
    """Point at a local dataset directory used instead of downloading.

    Passing ``None`` clears the override.
    """
    global _dataset_dir_override
    if path is None:
        _dataset_dir_override = None
    else:
        _dataset_dir_override = Path(path).expanduser().resolve()
    return _dataset_dir_override


def set_base_url(url: str | None) -> Optional[str]:
    """Override the remote base URL used for downloads."""
    global _base_url_override
    _base_url_override = url.rstrip("/") + "/" if url else None
    return _base_url_override


def get_cache_dir() -> Path:
    """Return (and create) the active cache directory."""
    if _cache_dir_override is not None:
        cache = _cache_dir_override
    else:
        env = os.environ.get("GOFCHIANTI_CACHE")
        if env:
            cache = Path(env).expanduser().resolve()
        else:
            cache = Path(platformdirs.user_cache_dir("gofchianti"))
    cache.mkdir(parents=True, exist_ok=True)
    return cache


def _base_url() -> str:
    if _base_url_override is not None:
        return _base_url_override
    env = os.environ.get("GOFCHIANTI_BASE_URL")
    if env:
        return env.rstrip("/") + "/"
    return _DEFAULT_BASE_URL


def _dataset_dir() -> Optional[Path]:
    if _dataset_dir_override is not None:
        return _dataset_dir_override
    env = os.environ.get("GOFCHIANTI_DATASET_DIR")
    if env:
        return Path(env).expanduser().resolve()
    return None


# --------------------------------------------------------------------------- #
# Hashing helpers
# --------------------------------------------------------------------------- #
def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _cache_path(name: str, kind: str) -> Path:
    sub = _KIND_SUBDIR.get(kind, "")
    base = get_cache_dir()
    return (base / sub / name) if sub else (base / name)


# --------------------------------------------------------------------------- #
# Manifest handling
# --------------------------------------------------------------------------- #
def _local_source(name: str) -> Optional[Path]:
    """Return a path to ``name`` inside the local dataset dir, if available.

    The local dataset directory is laid out with ``data/`` and ``abundance/``
    subfolders plus the catalogue/manifest at the top level, mirroring the
    converter output.
    """
    ddir = _dataset_dir()
    if ddir is None:
        return None
    for candidate in (ddir / name, ddir / "data" / name, ddir / "abundance" / name):
        if candidate.exists():
            return candidate
    return None


def load_manifest(force: bool = False) -> Dict[str, dict]:
    """Load the dataset manifest mapping ``name -> {kind, sha256, size}``.

    Tries, in order: a fresh in-process copy, the local dataset directory, the
    remote URL, then any previously cached manifest.  Returns an empty mapping
    if nothing can be found (callers then fall back to whatever is cached).
    """
    global _manifest_cache, _manifest_fetched_at
    now = time.time()
    if not force and _manifest_cache is not None and (now - _manifest_fetched_at) < _MANIFEST_TTL_SECONDS:
        return _manifest_cache

    raw: Optional[str] = None

    # 1) Local dataset directory (offline / dev).
    local = _local_source(MANIFEST_NAME)
    if local is not None:
        raw = local.read_text()

    # 2) Remote.
    if raw is None:
        try:
            import requests  # pooch dependency; always available

            resp = requests.get(_base_url() + MANIFEST_NAME, timeout=15)
            if resp.ok:
                raw = resp.text
                # Persist a copy so we can keep working offline later.
                cached = get_cache_dir() / MANIFEST_NAME
                cached.write_text(raw)
        except Exception:
            raw = None

    # 3) Previously cached manifest.
    if raw is None:
        cached = get_cache_dir() / MANIFEST_NAME
        if cached.exists():
            raw = cached.read_text()

    if raw is None:
        _manifest_cache = {}
        _manifest_fetched_at = now
        return {}

    data = json.loads(raw)
    files = data.get("files", [])
    mapping = {entry["name"]: entry for entry in files}
    _manifest_cache = mapping
    _manifest_fetched_at = now
    return mapping


# --------------------------------------------------------------------------- #
# Fetching
# --------------------------------------------------------------------------- #
def fetch_file(name: str, kind: str, *, expected_sha256: Optional[str] = None) -> Path:
    """Return a local path to ``name``, downloading/copying it if necessary.

    Parameters
    ----------
    name : str
        Flat asset filename, e.g. ``fe_12_195.119_v-11.0.2.npz``.
    kind : str
        One of ``"data"``, ``"abundance"``, ``"catalog"``, ``"manifest"``.
    expected_sha256 : str, optional
        If given, overrides the hash looked up from the manifest.
    """
    dest = _cache_path(name, kind)
    dest.parent.mkdir(parents=True, exist_ok=True)

    # Resolve the expected hash from the manifest when not supplied.
    if expected_sha256 is None:
        entry = load_manifest().get(name)
        if entry is not None:
            expected_sha256 = entry.get("sha256")

    # If a cached copy already matches the expected hash, use it as-is.
    if dest.exists() and expected_sha256 and _sha256(dest) == expected_sha256:
        return dest

    # 1) Local dataset directory (copy in).
    local = _local_source(name)
    if local is not None:
        if not (dest.exists() and local.samefile(dest)):
            shutil.copy2(local, dest)
        if expected_sha256 and _sha256(dest) != expected_sha256:
            raise ValueError(
                f"Hash mismatch for local dataset file {name!r} "
                f"(expected {expected_sha256})."
            )
        return dest

    # 2) Remote download via pooch (handles hashing + atomic writes).
    url = _base_url() + name
    try:
        import pooch

        known_hash = f"sha256:{expected_sha256}" if expected_sha256 else None
        fetched = pooch.retrieve(
            url=url,
            known_hash=known_hash,
            fname=name,
            path=str(dest.parent),
            progressbar=False,
        )
        return Path(fetched)
    except Exception as exc:
        # 3) Last resort: an existing cached copy, even if we could not verify it.
        if dest.exists():
            return dest
        raise FileNotFoundError(
            f"Could not obtain {name!r} from local dataset dir, remote "
            f"({url}), or cache. Last error: {exc}"
        ) from exc


def fetch_data_file(filename: str) -> Path:
    """Fetch a precomputed ``G(n_e, T)`` ``.npz`` file."""
    return fetch_file(filename, "data")


def fetch_abundance_file(filename: str) -> Path:
    """Fetch a CHIANTI ``.abund`` abundance file."""
    return fetch_file(filename, "abundance")


def fetch_catalog() -> Path:
    """Fetch the line catalogue parquet file."""
    return fetch_file(CATALOG_NAME, "catalog")


# --------------------------------------------------------------------------- #
# Bulk / maintenance
# --------------------------------------------------------------------------- #
def download_all(kinds: Optional[List[str]] = None, *, verbose: bool = True) -> List[Path]:
    """Download the entire dataset so the tool works fully offline.

    Parameters
    ----------
    kinds : list of str, optional
        Restrict to specific kinds (``"data"``, ``"abundance"``, ``"catalog"``).
        Defaults to everything in the manifest.
    """
    manifest = load_manifest(force=True)
    if not manifest:
        raise RuntimeError(
            "No manifest available; cannot enumerate the dataset to download. "
            "Set a base URL or local dataset directory first."
        )
    paths: List[Path] = []
    for name, entry in manifest.items():
        kind = entry.get("kind", "data")
        if kinds is not None and kind not in kinds:
            continue
        if verbose:
            print(f"[gofchianti] fetching {name} ...")
        paths.append(fetch_file(
            name, kind, expected_sha256=entry.get("sha256")))
    if verbose:
        print(
            f"[gofchianti] dataset ready in {get_cache_dir()} ({len(paths)} files).")
    return paths


def clear_cache(name: Optional[str] = None, kind: str = "data") -> None:
    """Delete cached files.

    With no arguments the whole cache directory is removed.  Pass ``name`` (and
    optionally ``kind``) to remove a single file.
    """
    if name is None:
        cache = get_cache_dir()
        if cache.exists():
            shutil.rmtree(cache)
        return
    path = _cache_path(name, kind)
    if path.exists():
        path.unlink()

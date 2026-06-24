"""Shared test fixtures.

Tests run fully offline against the locally built dataset (``GofChianti/dataset``)
by pointing ``GOFCHIANTI_DATASET_DIR`` at it and isolating the cache in a temp
directory.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
DATASET_DIR = REPO_ROOT / "dataset"
ABUND_SRC = Path("/usr/local/ssw/packages/chianti/dbase/abundance")
# IDL-generated ``.dat`` inputs (maintainer-only; not shipped in the package).
GOFNT_DAT_DIR = REPO_ROOT.parent / "gofnt"

# Make the maintainer tools importable (they live at the repo root, not in src).
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Use the local dataset for all fetches.
os.environ["GOFCHIANTI_DATASET_DIR"] = str(DATASET_DIR)


@pytest.fixture(autouse=True)
def isolated_cache(tmp_path, monkeypatch):
    """Give every test a clean, throwaway cache directory."""
    monkeypatch.setenv("GOFCHIANTI_CACHE", str(tmp_path / "cache"))
    import gofchianti.fetch as fetch

    # Reset module-level state that may persist between tests.
    fetch._cache_dir_override = None
    fetch._manifest_cache = None
    fetch._manifest_fetched_at = 0.0
    yield


@pytest.fixture(scope="session")
def dataset_available():
    if not (DATASET_DIR / "manifest.json").exists():
        pytest.skip(
            "dataset not built; run maintainers/convert_dat_to_npz.py first")
    return DATASET_DIR


@pytest.fixture(scope="session")
def dat_dir():
    """Directory of IDL-generated ``*_gofnt_v-*.dat`` inputs (maintainer tests)."""
    if not GOFNT_DAT_DIR.exists() or not any(GOFNT_DAT_DIR.glob("*_gofnt_v-*.dat")):
        pytest.skip("IDL .dat inputs not present")
    return GOFNT_DAT_DIR


@pytest.fixture(scope="session")
def abund_available():
    """CHIANTI abundance source directory (skips when SSW is unavailable)."""
    if not ABUND_SRC.exists():
        pytest.skip("CHIANTI abundance source not present")
    return ABUND_SRC

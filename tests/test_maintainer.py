"""Tests for the maintainer converter's release/build helpers.

These run fully offline: the dry-run upload path never touches the network, and
the build test writes into a temporary directory.
"""

from __future__ import annotations

import shutil
import subprocess

import pytest

from maintainers.convert_dat_to_npz import (
    _release_assets,
    build_dataset,
    upload_release,
)


def test_release_assets_order_and_contents(dataset_available):
    assets = _release_assets(dataset_available)
    names = [a.name for a in assets]
    # catalogue + manifest are always last, in that order.
    assert names[-2:] == ["catalog.parquet", "manifest.json"]
    # 56 per-line .npz files in the canonical dataset.
    assert sum(n.endswith(".npz") for n in names) == 56
    # at least one abundance file is included.
    assert sum(n.endswith(".abund") for n in names) >= 1


def test_upload_missing_assets_raises(tmp_path):
    """Validation happens before any remote/gh interaction."""
    with pytest.raises(FileNotFoundError):
        upload_release(tmp_path, "owner/repo", "tag", dry_run=True)


def test_upload_dry_run_does_not_invoke_subprocess(dataset_available, monkeypatch):
    if shutil.which("gh") is None:
        pytest.skip("GitHub CLI 'gh' not available")

    def fail(*args, **kwargs):  # pragma: no cover - must never run
        raise AssertionError("subprocess.run must not be called during dry-run")

    monkeypatch.setattr(subprocess, "run", fail)
    assets = upload_release(dataset_available, "owner/repo", "tag", dry_run=True)
    # 56 npz + >=1 abund + catalog + manifest.
    assert len(assets) >= 59
    assert [a.name for a in assets][-2:] == ["catalog.parquet", "manifest.json"]


def test_build_dataset_into_tmp_no_bundle(tmp_path, dat_dir, abund_available):
    out = build_dataset(dat_dir, abund_available, tmp_path / "ds",
                        package_data_dir=None, verbose=0)
    assert (out / "manifest.json").exists()
    assert (out / "catalog.parquet").exists()
    assert len(list((out / "data").glob("*.npz"))) == 56
    # bundling disabled → nothing copied into the package data dir.
    assert not (out / "src").exists()

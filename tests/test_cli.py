"""Tests for the ``gofchianti-build-dataset`` console entry point."""

from __future__ import annotations

import gofchianti.cli as cli


def test_main_is_callable():
    assert callable(cli.main)


def test_main_returns_2_on_converter_failure(monkeypatch):
    """Exceptions from the converter are caught and surfaced as exit code 2."""
    import maintainers.convert_dat_to_npz as conv

    def boom(argv=None):
        raise RuntimeError("simulated failure")

    monkeypatch.setattr(conv, "main", boom)
    assert cli.main([]) == 2


def test_main_builds_dataset_and_skips_bundle(tmp_path, dat_dir, abund_available, monkeypatch):
    """End-to-end build through the CLI, with bundling disabled via empty string.

    Also guards the ``--package-data-dir ""`` sentinel fix: no stray
    ``catalog.parquet`` should be written into the current directory.
    """
    monkeypatch.chdir(tmp_path)
    out_dir = tmp_path / "ds"
    rc = cli.main([
        "--dat-dir", str(dat_dir),
        "--abund-src", str(abund_available),
        "--out-dir", str(out_dir),
        "--package-data-dir", "",
    ])
    assert rc == 0
    assert (out_dir / "manifest.json").exists()
    assert not (tmp_path / "catalog.parquet").exists()

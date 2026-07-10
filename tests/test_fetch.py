"""Tests for the download base-URL configuration.

The default download location is the IAS SPICE data server. The URL is
overridable at runtime (``set_base_url``) or via the ``GOFCHIANTI_BASE_URL``
environment variable, always normalised to a single trailing slash.
"""

from __future__ import annotations

import gofchianti as gc
import gofchianti.fetch as fetch


def test_default_base_url_points_to_ias():
    assert "spice.osups.universite-paris-saclay.fr" in fetch._DEFAULT_BASE_URL
    assert "contribution_functions" in fetch._DEFAULT_BASE_URL
    assert fetch._DEFAULT_BASE_URL.endswith("/")


def test_base_url_resolves_to_default(monkeypatch):
    monkeypatch.delenv("GOFCHIANTI_BASE_URL", raising=False)
    fetch.set_base_url(None)
    try:
        assert fetch._base_url() == fetch._DEFAULT_BASE_URL
        assert fetch._base_url().endswith("/")
    finally:
        fetch.set_base_url(None)


def test_set_base_url_overrides_and_normalises(monkeypatch):
    monkeypatch.delenv("GOFCHIANTI_BASE_URL", raising=False)
    try:
        gc.set_base_url("https://example.test/gofchianti")
        assert fetch._base_url() == "https://example.test/gofchianti/"
    finally:
        fetch.set_base_url(None)


def test_env_var_base_url(monkeypatch):
    fetch.set_base_url(None)
    monkeypatch.setenv("GOFCHIANTI_BASE_URL", "https://env.test/data")
    try:
        assert fetch._base_url() == "https://env.test/data/"
    finally:
        monkeypatch.delenv("GOFCHIANTI_BASE_URL", raising=False)
        fetch.set_base_url(None)

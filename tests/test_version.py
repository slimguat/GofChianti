"""Tests for dynamic version resolution."""

from __future__ import annotations

import gofchianti as gc


def test_version_callable_returns_nonempty_string():
    v = gc.version()
    assert isinstance(v, str)
    assert v.strip()


def test_dunder_version_matches_version_call():
    assert gc.__version__ == gc.version()


def test_version_looks_like_a_version():
    # Expect at least a leading numeric component (e.g. "0.1.0").
    assert gc.version()[0].isdigit()

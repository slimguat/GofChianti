"""Tests for the shared ``_vprint`` verbosity helper."""

from __future__ import annotations

from gofchianti.utils import _vprint


def test_vprint_suppressed_when_below_level(capsys):
    """Nothing is printed when ``verbose`` is below the message ``level``."""
    _vprint(0, 1, "should-not-appear")
    assert capsys.readouterr().out == ""


def test_vprint_emitted_when_at_or_above_level(capsys):
    """The message and its label are printed once the threshold is met."""
    _vprint(1, 1, "hello")
    out = capsys.readouterr().out
    assert "hello" in out
    assert "[Verbose]" in out


def test_vprint_info_label(capsys):
    _vprint(0, 0, "info-message")
    out = capsys.readouterr().out
    assert "[Info]" in out
    assert "info-message" in out


def test_vprint_warning_always_shown_at_default_verbosity(capsys):
    """Warnings (level -1) appear even at the lowest verbosity (0)."""
    _vprint(0, -1, "warn-message")
    out = capsys.readouterr().out
    assert "[Warning]" in out
    assert "warn-message" in out


def test_vprint_non_integer_verbose_treated_as_zero(capsys):
    """A non-int ``verbose`` degrades gracefully to 0 (info shows, verbose hides)."""
    _vprint("not-an-int", 0, "info")
    _vprint("not-an-int", 1, "verbose")
    out = capsys.readouterr().out
    assert "info" in out
    assert "verbose" not in out


def test_vprint_unknown_level_uses_generic_label(capsys):
    _vprint(9, 9, "deep")
    out = capsys.readouterr().out
    assert "[Level_9]" in out
    assert "deep" in out

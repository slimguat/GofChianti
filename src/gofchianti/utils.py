from __future__ import annotations


def _vprint(verbose: int, level: int, *args, **kwargs) -> None:
    """
    Print a message at a given verbosity level, prefixed by a label.

    Parameters
    ----------
    verbose : int
        Current verbosity setting.
    level : int
        The threshold level for this message:
          -1 → "Warning"
           0 → "Info"
           1 → "Verbose"
           2 → "Debug"
           3 → "Debug_Plot"
           4 → "Debug_Plot_Save"
    *args, **kwargs
        Passed to built-in print() after the prefix.
    """
    try:
        v = int(verbose)
    except Exception:
        v = 0
    if v < level:
        return

    labels = {
        -1: ("Warning", "\033[91m"),
        0: ("Info", "\033[0m"),
        1: ("Verbose", "\033[92m"),
        2: ("Debug", "\033[90m"),
        3: ("Debug_Plot", "\033[90m"),
        4: ("Debug_Plot_Save", "\033[90m"),
    }
    prefix, color = labels.get(level, (f"Level_{level}", "\033[0m"))
    reset = "\033[0m"
    print(f"{color}[{prefix}]{reset}", *args, **kwargs)

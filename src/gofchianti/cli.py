from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional, Sequence


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Wrapper entry point for building the GofChianti dataset.

    This attempts to import the maintainer converter from the repository and
    invoke its `main` function. It is intended to be used from a source
    checkout (or editable install).
    """
    # Ensure the repository root is on sys.path so the top-level `maintainers`
    # package can be imported when running from the checkout.
    repo_root = Path(__file__).resolve().parents[2]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    try:
        from maintainers.convert_dat_to_npz import main as converter_main  # type: ignore

        # Converter expects argv: Optional[List[str]]
        converter_main(list(argv) if argv is not None else None)
        return 0
    except Exception as exc:
        print("Could not run maintainer converter:", exc, file=sys.stderr)
        print("Ensure you run this command from the repository checkout where the maintainer scripts live.", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

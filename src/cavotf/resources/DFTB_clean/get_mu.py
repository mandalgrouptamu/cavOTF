"""Placeholder μ calculation script.

This default script simply writes a zero-valued ``dmu.dat`` so the workflow
can proceed. Replace its contents with your actual μ calculation.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np


def main() -> None:
    outfile = Path("dmu.dat")
    if not outfile.exists():
        outfile.write_text("0.0\n")
    print(f"Wrote placeholder mu to {outfile.resolve()}")


if __name__ == "__main__":
    main()

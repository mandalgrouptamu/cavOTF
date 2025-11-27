from __future__ import annotations

import argparse
from pathlib import Path

from .config import CavOTFConfig
from .workflow import CavOTFWorkflow


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Run the cavOTF workflow from a single configuration file.")
    parser.add_argument(
        "--config",
        default="input.txt",
        type=Path,
        help="Path to the configuration file (default: input.txt)",
    )

    args = parser.parse_args(argv)
    config = CavOTFConfig.from_file(args.config)
    workflow = CavOTFWorkflow(config)
    workflow.run_all()


if __name__ == "__main__":
    main()

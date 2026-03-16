# =============================================================================
#  Project:     cavOTF.py
#  File:        cli.py
#  Author:      Amir H. Amini <amiramini@tamu.edu>
#  Modified by: Sachith Wickramasinghe
#  Last update: 03/16/2026
#
#  Description:
#      Command-line interface for cavOTF.py simulations, including argument
#      parsing and workflow execution.
# =============================================================================

from __future__ import annotations

import time
import sys

import argparse
import logging
from pathlib import Path
from typing import Optional

from .config import ConfigError, load_config
from .workflow import run_workflow, validate_workflow


LOG_LEVELS = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
}


def _print_banner() -> None:
        banner = r"""
▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒
╔═══════════════════════════════════════════════════════════╗
║     ██████╗ █████╗ ██╗   ██╗ ██████╗ ████████╗███████╗    ║
║    ██╔════╝██╔══██╗██║   ██║██╔═══██╗╚══██╔══╝██╔════╝    ║
║    ██║     ███████║██║   ██║██║   ██║   ██║   █████╗      ║
║    ██║     ██╔══██║╚██╗ ██╔╝██║   ██║   ██║   ██╔══╝      ║
║    ╚██████╗██║  ██║ ╚████╔╝ ╚██████╔╝   ██║   ██║         ║
║     ╚═════╝╚═╝  ╚═╝  ╚═══╝   ╚═════╝    ╚═╝   ╚═╝         ║
╠═══════════════════════════════════════════════════════════╣
║      Mandal Group • Texas A&M University • 2025-2026      ║
╚═══════════════════════════════════════════════════════════╝
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 If you use this software (in whole or in part), please cite:
 DOI: 10.1039/D5CP04879F.
▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒
"""
        for line in banner.splitlines():
            print(line)
            sys.stdout.flush()
            time.sleep(0.05)   # adjust speed (seconds per line)



def _configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="[%(levelname)s] %(message)s")


def _parse_args(args: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="CAVOTF cavity MD workflow runner")
    parser.add_argument("--config", required=True, type=Path, help="Path to input.txt configuration file")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")

    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Execute the workflow")
    run_parser.set_defaults(command="run")

    validate_parser = subparsers.add_parser("validate", help="Validate configuration and print planned actions")
    validate_parser.set_defaults(command="validate")
    
    extend_parser = subparsers.add_parser("extend", help="Continue an existing workflow from saved midpoint data")
    extend_parser.set_defaults(command="extend")

    return parser.parse_args(args)


def main(argv: Optional[list[str]] = None) -> None:
    _print_banner()
    args = _parse_args(argv)
    _configure_logging(args.verbose)

    try:
        config = load_config(args.config)
    except ConfigError as exc:  # noqa: BLE001
        logging.error("Configuration error: %s", exc)
        raise SystemExit(1) from exc

    if args.command == "validate":
        validate_workflow(config)
    elif args.command == "run":
        run_workflow(config)
    elif args.command == "extend":
        run_workflow(config, extend=True)
    else:  # pragma: no cover - argparse ensures coverage
        raise SystemExit(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()

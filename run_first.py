from __future__ import annotations

from pathlib import Path

from cavotf.config import CavOTFConfig
from cavotf.workflow import CavOTFWorkflow


def main(config_path: Path | str = "input.txt") -> None:
    config = CavOTFConfig.from_file(config_path)
    workflow = CavOTFWorkflow(config)
    workflow.stage_clients()


if __name__ == "__main__":
    main()

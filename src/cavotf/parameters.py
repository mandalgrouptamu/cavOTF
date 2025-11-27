from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType
from typing import Tuple


def load_param(clean_template_dir: Path) -> Tuple[object, ModuleType]:
    func_path = clean_template_dir / "funcLM.py"
    if not func_path.exists():
        raise FileNotFoundError(f"Unable to locate funcLM.py in {clean_template_dir}")

    spec = importlib.util.spec_from_file_location("cavotf.funcLM", func_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load funcLM.py from {func_path}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    if not hasattr(module, "param"):
        raise AttributeError("funcLM.py does not define a 'param' callable")

    return module.param(), module

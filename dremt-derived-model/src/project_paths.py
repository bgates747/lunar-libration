# SPDX-License-Identifier: MIT
"""Portable paths anchored at the lunar-libration clone root, never at cwd."""
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
MODEL_ROOT = PROJECT_ROOT / "dremt-derived-model"
DATA_DIR = MODEL_ROOT / "data"
CONFIG_DIR = MODEL_ROOT / "config"
OUTPUT_DIR = MODEL_ROOT / "outputs" / "validation"
ANIMATION_OUTPUT_DIR = MODEL_ROOT / "outputs" / "animation"


def project_path(value: str | Path) -> Path:
    """Resolve a repository-relative input/output; accept explicit absolute paths."""
    path = Path(value)
    return (PROJECT_ROOT / path).resolve() if not path.is_absolute() else path.resolve()

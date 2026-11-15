"""Configuration and path access.

Thresholds live in config/config.yaml, family definitions in
config/families.yaml. Nothing in this package hard-codes a cutoff: the thesis
has to quote these numbers, and several of them are ablation parameters, so
there must be exactly one place to change them.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = REPO_ROOT / "config"


@dataclass(frozen=True)
class Config:
    """Parsed configuration, with paths resolved against the repository root."""

    raw: dict[str, Any]
    families: dict[str, Any]

    def __getitem__(self, key: str) -> Any:
        return self.raw[key]

    def get(self, key: str, default: Any = None) -> Any:
        return self.raw.get(key, default)

    def path(self, key: str) -> Path:
        """Resolve a `paths.<key>` entry against the repository root."""
        return REPO_ROOT / self.raw["paths"][key]

    @property
    def active_families(self) -> list[str]:
        """Families in scope for this round.

        Narrowed at the end of stage 1 if the data profile is thinner than
        expected -- that decision is made in config/families.yaml, not in code.
        """
        return list(self.families["active_families"])

    def family_def(self, family: str) -> dict[str, Any]:
        return self.families["families"][family]

    def family_label(self, family: str, language: str = "en") -> str:
        """Display name for tables and figures, in Turkish or English."""
        return self.family_def(family)[f"name_{language}"]

    @property
    def seeds(self) -> list[int]:
        return list(self.raw["project"]["random_seeds"])


@lru_cache(maxsize=1)
def load_config(
    config_path: Path | None = None,
    families_path: Path | None = None,
) -> Config:
    """Load and cache the configuration.

    Cached because it is read from nearly every module; pass explicit paths in
    tests to bypass the cache for a fixture.
    """
    config_path = config_path or CONFIG_DIR / "config.yaml"
    families_path = families_path or CONFIG_DIR / "families.yaml"

    with Path(config_path).open(encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    with Path(families_path).open(encoding="utf-8") as handle:
        families = yaml.safe_load(handle)

    unknown = set(families["active_families"]) - set(families["families"])
    if unknown:
        raise ValueError(
            f"active_families lists undefined families: {sorted(unknown)}"
        )

    return Config(raw=raw, families=families)

"""Portable repository path configuration.

Public notebooks use repository-relative defaults. Large local datasets can stay
outside the repository and be selected with environment variables, avoiding
machine-specific paths in committed code.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path


def _find_repository_root(start: Path | None = None) -> Path:
    current = (start or Path.cwd()).resolve()
    for candidate in (current, *current.parents):
        if (candidate / "pyproject.toml").is_file() and (
            candidate / "documentation_code"
        ).is_dir():
            return candidate
    raise FileNotFoundError(
        "Could not locate the repository root. Start inside the repository or "
        "set SCALE_AWARE_ST_REPO_ROOT."
    )


def _env_path(name: str, default: Path) -> Path:
    value = os.environ.get(name)
    return Path(value).expanduser().resolve() if value else default.resolve()


@dataclass(frozen=True)
class RepositoryConfig:
    """Resolved paths shared by notebooks and scripts."""

    root: Path
    data_dir: Path
    aging_gene_dir: Path
    additional_data_dir: Path
    external_data_dir: Path
    results_dir: Path
    resources_dir: Path

    @classmethod
    def from_env(cls, start: Path | None = None) -> "RepositoryConfig":
        root_value = os.environ.get("SCALE_AWARE_ST_REPO_ROOT")
        root = (
            Path(root_value).expanduser().resolve()
            if root_value
            else _find_repository_root(start)
        )
        data_dir = _env_path("SCALE_AWARE_ST_DATA_DIR", root / "data")
        return cls(
            root=root,
            data_dir=data_dir,
            aging_gene_dir=_env_path(
                "SCALE_AWARE_ST_AGING_GENE_DIR",
                root / "external_data" / "published_aging_gene_sets",
            ),
            additional_data_dir=_env_path(
                "SCALE_AWARE_ST_ADDITIONAL_DATA_DIR",
                root / "Supplementary tables",
            ),
            external_data_dir=_env_path(
                "SCALE_AWARE_ST_EXTERNAL_DATA_DIR", root / "external_data"
            ),
            results_dir=_env_path("SCALE_AWARE_ST_RESULTS_DIR", root / "results"),
            resources_dir=_env_path(
                "SCALE_AWARE_ST_RESOURCES_DIR", root / "resources"
            ),
        )

    def ensure_output_dirs(self) -> None:
        """Create only output directories; input locations are never created."""

        self.results_dir.mkdir(parents=True, exist_ok=True)

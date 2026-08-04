from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_conda_environment_does_not_install_relative_project_path():
    environment = (
        ROOT / "environment_exports" / "environment-core.yml"
    ).read_text(encoding="utf-8")
    readme = (ROOT / "README.md").read_text(encoding="utf-8")

    assert "- -e ." not in environment
    assert "python -m pip install -e ." in readme
    assert "directory containing\n`pyproject.toml`" in readme


def test_r_tutorial_preserves_long_cell_identifiers():
    tutorial = (
        ROOT
        / "documentation_code"
        / "aldex_main_crossstudy_repo"
        / "tutorial"
        / "tutorial_ALDEx3.Rmd"
    ).read_text(encoding="utf-8")

    assert 'colClasses = "character"' in tutorial
    assert "cell_ids <- metadata[[1]]" in tutorial
    assert "metadata[[1]] <- NULL" in tutorial
    assert "metadata <- metadata[, -1, drop = FALSE]" not in tutorial
    assert "metadata[] <- lapply(metadata, type.convert, as.is = TRUE)" in tutorial


def test_r_tutorial_uses_portable_model_and_output_helpers():
    tutorial_dir = (
        ROOT / "documentation_code" / "aldex_main_crossstudy_repo" / "tutorial"
    )
    tutorial = (tutorial_dir / "tutorial_ALDEx3.Rmd").read_text(encoding="utf-8")
    installer = (tutorial_dir / "install_tutorial_dependencies.R").read_text(
        encoding="utf-8"
    )

    assert "reformulas::nobars(r_formula)" in tutorial
    assert 'normalizePath(output_file, winslash = "/", mustWork = TRUE)' in tutorial
    assert "8c05ad40c41279dffa05dc808167ffcd53207740" in installer

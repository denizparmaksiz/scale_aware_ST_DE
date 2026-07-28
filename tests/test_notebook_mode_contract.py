import ast
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AFFECTED_NOTEBOOKS = [
    ROOT / "documentation_code" / name
    for name in (
        "03_main_merfish_aldex_prepost.ipynb",
        "04_main_MERFISH_bio_figs.ipynb",
        "05_primaryMERFISH_functional_annotations.ipynb",
        "06_Sun_main_binary_DE.ipynb",
        "07_sun_continuous_raw_validation.ipynb",
        "08_sun_supplementary_analyses.ipynb",
        "09_cosmx_figs_DE.ipynb",
    )
]


def notebook_code(path):
    notebook = json.loads(path.read_text(encoding="utf-8"))
    return [
        "".join(cell.get("source", []))
        for cell in notebook["cells"]
        if cell.get("cell_type") == "code"
    ]


def test_aldex_notebooks_have_valid_mode_preflights_and_python_syntax():
    for path in AFFECTED_NOTEBOOKS:
        cells = notebook_code(path)
        combined = "\n".join(cells)
        assert "missing_mode_inputs" in combined, path.name
        assert "RECOMPUTE_WORKFLOW.md" in combined, path.name
        for cell_number, source in enumerate(cells, start=1):
            ast.parse(source, filename=f"{path.name}:code-cell-{cell_number}")


def test_external_aldex_boundaries_reject_empty_workbook_collections():
    expected_messages = {
        "03_main_merfish_aldex_prepost.ipynb": "No per-stratum primary cell-type x region",
        "06_Sun_main_binary_DE.ipynb": "No individual Sun binary-age",
        "07_sun_continuous_raw_validation.ipynb": "No individual Sun continuous-age",
        "08_sun_supplementary_analyses.ipynb": "No individual Sun continuous subregion",
        "09_cosmx_figs_DE.ipynb": "No individual CosMx ALDEx",
    }
    for name, message in expected_messages.items():
        assert message in "\n".join(notebook_code(ROOT / "documentation_code" / name))

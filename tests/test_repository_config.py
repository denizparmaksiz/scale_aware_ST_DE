from pathlib import Path
from unittest.mock import patch

from scale_aware_st.config import RepositoryConfig


ROOT = Path(__file__).resolve().parents[1]


def test_default_data_layout_matches_extracted_zenodo_tree():
    with patch.dict(
        "os.environ",
        {
            "SCALE_AWARE_ST_REPO_ROOT": str(ROOT),
            "SCALE_AWARE_ST_DATA_DIR": "",
            "SCALE_AWARE_ST_ADDITIONAL_DATA_DIR": "",
        },
        clear=False,
    ):
        config = RepositoryConfig.from_env()

    assert config.data_dir == (ROOT / "data").resolve()
    assert config.additional_data_dir == (
        ROOT / "data" / "supplementary_files"
    ).resolve()


def test_notebook_04_ontology_resource_is_version_controlled():
    ontology = ROOT / "resources" / "allen_ccf" / "CCFv3OntologyStructure_u16.xlsx"
    assert ontology.is_file()

from pathlib import Path
import os
import unittest

import pandas as pd
from pandas.testing import assert_frame_equal

from scale_aware_st.aldex_results import (
    POOLED_ALDEX_SPECS,
    load_pooled_aldex_results,
)


ROOT = Path(__file__).resolve().parents[1]
ADDITIONAL_DIR = Path(
    os.environ.get(
        "SCALE_AWARE_ST_ADDITIONAL_DATA_DIR",
        ROOT / "data" / "supplementary_files",
    )
)


class TestPooledALDExResults(unittest.TestCase):
    def test_all_pooled_specs_load_and_preserve_unique_genes(self):
        if not ADDITIONAL_DIR.is_dir():
            self.skipTest(
                "Extract the Zenodo supplementary_files archive under data/ "
                "or set SCALE_AWARE_ST_ADDITIONAL_DATA_DIR."
            )
        for name, spec in POOLED_ALDEX_SPECS.items():
            with self.subTest(analysis=name):
                results = load_pooled_aldex_results(
                    ADDITIONAL_DIR / spec.additional_file,
                    sheet_name=spec.sheet_name,
                    stratum_column=spec.stratum_column,
                )
                self.assertTrue(results)
                for label, table in results.items():
                    self.assertTrue(label)
                    self.assertIn("gene", table.columns)
                    self.assertFalse(table["gene"].astype(str).duplicated().any())

    def test_representative_individual_matches_pooled_rows_exactly(self):
        individual_path = (
            ROOT
            / "data_for_deeper_test"
            / "representative_rawcount_metadata_aldexoutput"
            / "aldex_mem_cell_type_Isocortex_astro_results_07142025.xlsx"
        )
        if not individual_path.is_file() or not ADDITIONAL_DIR.is_dir():
            self.skipTest("Local deep-test and deposited supplementary inputs are absent.")
        individual = pd.read_excel(individual_path, engine="openpyxl")
        spec = POOLED_ALDEX_SPECS["primary_celltype_region"]
        pooled = load_pooled_aldex_results(
            ADDITIONAL_DIR / spec.additional_file,
            sheet_name=spec.sheet_name,
            stratum_column=spec.stratum_column,
        )["Isocortex_astro"]
        assert_frame_equal(individual, pooled, check_exact=True)


if __name__ == "__main__":
    unittest.main()

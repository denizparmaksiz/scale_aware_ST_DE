import unittest

import numpy as np
import pandas as pd

from scale_aware_st.aging_gene_sets import empirical_effect_cutoff


class TestAgingGeneSets(unittest.TestCase):
    def test_empirical_effect_cutoff_matches_original_definition(self):
        values = pd.Series([-5, -4, -3, -2, -1, 1, 2, 3, 4, 10, np.nan])
        expected = min(
            abs(np.nanquantile(values, 0.1)),
            abs(np.nanquantile(values, 0.9)),
        )
        self.assertEqual(empirical_effect_cutoff(values), expected)

    def test_all_missing_effects_fail_clearly(self):
        with self.assertRaisesRegex(ValueError, "all-missing"):
            empirical_effect_cutoff(pd.Series([np.nan, np.nan]))


if __name__ == "__main__":
    unittest.main()

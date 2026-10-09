"""A held-out sample must not influence panel selection."""
import unittest
import numpy as np
import pandas as pd
from src.analysis import select_cancer_panel

class ReplicationTest(unittest.TestCase):
    def test_heldout_expression_cannot_change_selected_panel(self):
        c = pd.DataFrame([[100, 20, 5], [120, 25, 4], [15, 40, 20], [20, 45, 20], [80, 20, 5]],
                         index=['ca1', 'ca2', 'n1', 'n2', 'heldout'], columns=['g1', 'g2', 'g3'])
        def norm(frame):
            return np.log1p(frame.div(frame.sum(axis=1), axis=0) * 1e6)
        before, effects_before = select_cancer_panel(c, norm(c), ['ca1', 'ca2'], ['n1', 'n2'])
        self.assertEqual(before, ['g1'])
        c.loc['heldout'] = [0, 10**9, 10**9]
        after, effects_after = select_cancer_panel(c, norm(c), ['ca1', 'ca2'], ['n1', 'n2'])
        self.assertEqual(before, after)
        pd.testing.assert_series_equal(effects_before, effects_after)

    def test_overlapping_groups_are_rejected(self):
        with self.assertRaises(ValueError):
            select_cancer_panel(pd.DataFrame(), pd.DataFrame(), ['same'], ['same'])

if __name__ == '__main__':
    unittest.main()

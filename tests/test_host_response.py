import unittest
import numpy as np
import pandas as pd
from src.host_response import contrast, viral_feature_matches


class HostResponseTest(unittest.TestCase):
    def test_omission_exposes_single_sample_driving_direction(self):
        values = pd.Series({'a': 0., 'b': 10., 'c': 4., 'd': 4.})
        result, folds = contrast(values, ['a', 'b'], ['c', 'd'])
        self.assertEqual(result['effect'], 1.)
        self.assertFalse(result['direction_retained'])
        self.assertEqual(result['omission_min'], -4.)
        self.assertEqual(result['omission_max'], 6.)
        self.assertEqual(len(folds), 4)

    def test_equal_sample_weights_and_reversing_groups(self):
        values = pd.Series({'a': 8., 'b': 12., 'c': 1., 'd': 3., 'e': 2.})
        forward, _ = contrast(values, ['a', 'b'], ['c', 'd', 'e'])
        reverse, _ = contrast(values, ['c', 'd', 'e'], ['a', 'b'])
        self.assertEqual(forward['effect'], 8.)
        self.assertTrue(forward['direction_retained'])
        self.assertEqual(reverse['effect'], -8.)
        self.assertEqual(reverse['omission_min'], -forward['omission_max'])

    def test_invalid_replicates_are_not_silently_used(self):
        values = pd.Series({'a': 1., 'b': 2., 'c': np.nan, 'd': 4.})
        for target, reference in [(['a'], ['b', 'd']), (['a', 'b'], ['a', 'd']),
                                  (['a', 'a'], ['b', 'd']), (['a', 'b'], ['c', 'd'])]:
            with self.assertRaises(ValueError):
                contrast(values, target, reference)

    def test_viral_label_scan_avoids_host_gene_name_collisions(self):
        features = pd.DataFrame({'symbol': ['HPV16_E6', 'HPV-18-E7', 'human_papillomavirus',
                                             'E6', 'UBE3A', 'HPVBP', 'TP53']})
        self.assertEqual(viral_feature_matches(features).tolist(), [True, True, True, False, False, False, False])


if __name__ == '__main__':
    unittest.main()

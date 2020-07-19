import unittest
import numpy as np

import rescore


class TestRescore(unittest.TestCase):
    def setUp(self):
        exp_spec = {
            'm/z array': np.array([100.01, 103.409, 103.410, 105, 113.32]),
            'intensity array': np.array([0.5, 0.2, 0.3, 1.0, 0.35])
        }
        pred_spec1 = {
            'm/z array': np.array([100.05, 102.0, 103.40, 104.95, 104.99, 105.05, 112.0]),
            'intensity array': np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]),
            'ion array': np.array(["y1+1", "y3-H2O+1", "y3", "b2-NH3+1", "b2-H2O+1", "b4", "y6"])
        }
        binsize = 0.1
        self.dp, self.n_b, self.n_y = rescore.calc_dot_product_and_ions(exp_spec, pred_spec1, binsize)
        # Ions that should match: y1+1, y3 (twice), b2-NH3+1, b2-H2O+1, b4 (all three to the same peak, so divide intensity by 3)
        # Thus the expected DP is:
        # (0.1*0.5 + 0.15*0.2+0.15*0.3 + (0.4*1.0+0.5*1.0+0.6*1.0)/3)/(||exp_intensities||*||pred_intensities||)
        self.dp_expected = 0.625 / (
                np.linalg.norm(exp_spec['intensity array']) * np.linalg.norm(pred_spec1['intensity array']))
        # Dot Product of a spectrum with itself should give 1 if the mass tolerance is small enough
        self.auto_dp, _, _ = rescore.calc_dot_product_and_ions(pred_spec1, pred_spec1, binsize=0.001)

    def test_ny(self):
        self.assertEqual(self.n_y, 2, "Number of y ions does not match! Expect: 2, Got {}".format(self.n_y))

    def test_nb(self):
        self.assertEqual(self.n_b, 2, "Number of b ions does not match! Expect: 2, Got {}".format(self.n_b))

    def test_dp(self):
        self.assertEqual(self.dp, self.dp_expected,
                         "Dot Product does not match!, Expect {:.4f}, Got {:.4f}".format(self.dp_expected, self.dp))
    def test_dp_self(self):
        self.assertEqual(self.auto_dp, 1,
                         "Dot Product of a spectrum with itself should give 1 if binsize is small enough!")


class TestUtils(unittest.TestCase):
    #TODO
    pass

if __name__ == '__main__':
    unittest.main()

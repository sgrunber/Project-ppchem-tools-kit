# test_error_calculation.py
import unittest

import sys
sys.path.insert(0, "./src")
from project_ppchem_tools_kit.error_calculation_interface import calculate_error_propagation

class TestErrorCalculation(unittest.TestCase):

    def test_calculate_error_propagation(self):
        derivatives = [1.0, 2.0, 3.0]
        uncertainties = [0.1, 0.2, 0.3]
        standard_dev, average = calculate_error_propagation(derivatives, uncertainties)
        
        self.assertAlmostEqual(standard_dev, 0.34641016151377546)
        self.assertAlmostEqual(average, 2.0)
    
    def test_calculate_error_propagation_zero_division(self):
        derivatives = [1.0, 0.0, 3.0]
        uncertainties = [0.1, 0.2, 0.3]
        with self.assertRaises(ZeroDivisionError):
            calculate_error_propagation(derivatives, uncertainties)

    def test_calculate_error_propagation_empty_lists(self):
        derivatives = []
        uncertainties = []
        with self.assertRaises(ZeroDivisionError):
            calculate_error_propagation(derivatives, uncertainties)

if __name__ == '__main__':
    unittest.main()

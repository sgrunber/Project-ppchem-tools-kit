import unittest
from unittest.mock import Mock, patch
import pandas as pd 

import sys
sys.path.insert(0, "./src")
from Chem_pack.linear_regression import linear_regression


class TestLinearRegression(unittest.TestCase):
    @patch('Chem_pack.linear_regression.pd.read_excel')
    @patch('Chem_pack.linear_regression.display_graph')
    def test_linear_regression(self, mock_display_graph, mock_read_excel):
        # Mock the return value of pd.read_excel
        mock_read_excel.return_value = pd.DataFrame({
            'X': [1, 2, 3, 4, 5],
            'Y': [2, 4, 6, 8, 10]
        })
        
        # Create a mock for entry_input
        mock_entry_input = Mock()
        mock_entry_input.get.return_value = 'mock_path.xlsx'
        
        # Call the function with the mock object
        linear_regression(mock_entry_input, "X", "Y", "Test")

        # Check if display_graph was called once
        mock_display_graph.assert_called_once()
if __name__ == '__main__':
    unittest.main()

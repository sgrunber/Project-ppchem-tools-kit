
import unittest
from unittest.mock import patch
import sys
from pathlib import Path

sys.path.append('/Users/meloenzinger/Documents/GitHub/Project-ppchem-tools-kit-bis')
import Project_ppchem
from Project_ppchem import linear_regression
from Project_ppchem import relative_to_assets, ASSETS_PATH
from Project_ppchem import name_to_smiles

'''def run_tests():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(RelativePathTest))
    runner = unittest.TextTestRunner()
    runner.run(suite)'''

#RELATIVE TO ASSETS
class TestRelativePath(unittest.TestCase):
    def test_relative_to_assets(self):
        # Relative path to be tested
        test_path = Project_ppchem.relative_to_assets
        
        # Expected result
        expected_path = ASSETS_PATH / Path(test_path)
        
        # Actual result from the function
        result_path = relative_to_assets(test_path)
        
        # Assert that the constructed path is correct
        self.assertEqual(result_path, expected_path, f"Expected path '{expected_path}', got '{result_path}' instead")

# If the test is run as main, run the tests
if __name__ == '__main__':
    unittest.main()
    
    '''suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestRelativePath))
    runner = unittest.TextTestRunner()
    runner.run(suite)

#LINEAR REGRESSION
class TestLinearRegression(unittest.TestCase):
    @patch('Project_ppchem.pd.read_excel')
    def test_file_not_found(self, mock_read_excel):
        mock_read_excel.side_effect = FileNotFoundError
        with self.assertRaises(FileNotFoundError):
            linear_regression("invalid_path.xlsx", "X", "Y", "Test Plot")

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestLinearRegression))
    runner = unittest.TextTestRunner()
    runner.run(suite)
'''


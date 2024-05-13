'''
import unittest
from unittest.mock import patch
import sys
from pathlib import Path

sys.path.append('/Users/meloenzinger/Documents/GitHub/Project-ppchem-tools-kit-bis')
import Project_ppchem
from Project_ppchem import linear_regression
from Project_ppchem import relative_to_assets, ASSETS_PATH
from Project_ppchem import name_to_smiles

def run_tests():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(RelativePathTest))
    runner = unittest.TextTestRunner()
    runner.run(suite)

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
    
    suite = unittest.TestSuite()
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

def setUp(self):
        # Setting up mocks for GUI elements
        self.mock_entry = MagicMock()
        self.mock_messagebox = MagicMock()
        self.mock_radiobutton = MagicMock(return_value="1")

        # Patching the Tkinter components
        self.patcher_entry = patch('tkinter.Entry', return_value=self.mock_entry)
        self.patcher_messagebox = patch('tkinter.messagebox', self.mock_messagebox)
        self.patcher_radiobutton = patch('tkinter.IntVar.get', self.mock_radiobutton)

        self.patcher_entry.start()
        self.patcher_messagebox.start()
        self.patcher_radiobutton.start()

    def tearDown(self):
        self.patcher_entry.stop()
        self.patcher_messagebox.stop()
        self.patcher_radiobutton.stop()

    def test_empty_input(self):
        """Test process_input with empty input."""
        self.mock_entry.get.return_value = ''
        process_input()
        self.mock_messagebox.showerror.assert_called_once_with("Error", "Please enter a molecule name, SMILES code, or file path.")



def setUp(self):
        # Mocking the Entry widget and other necessary GUI components
        self.mock_entry_input = MagicMock()
        self.mock_selected_radio = MagicMock()
        self.mock_messagebox = MagicMock()

        # Patching GUI components
        self.patch_entry = patch('tkinter.Entry', return_value=self.mock_entry_input)
        self.patch_radio = patch('tkinter.IntVar', return_value=self.mock_selected_radio)
        self.patch_msgbox = patch('tkinter.messagebox', self.mock_messagebox)

        self.mock_entry_input.get.return_value = ''

        self.entry = self.patch_entry.start()
        self.radio = self.patch_radio.start()
        self.msgbox = self.patch_msgbox.start()

    def tearDown(self):
        self.patch_entry.stop()
        self.patch_radio.stop()
        self.patch_msgbox.stop()

    def test_empty_input(self):
        """Test process_input with empty input."""
        # Set radio value to trigger a specific path, if needed
        self.mock_selected_radio.get.return_value = '1'
        
        process_input()
        self.mock_messagebox.showerror.assert_called_once_with("Error", "Please enter a molecule name, SMILES code, or file path.")  
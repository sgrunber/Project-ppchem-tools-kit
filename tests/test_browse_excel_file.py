import unittest
from unittest.mock import patch, MagicMock
from io import StringIO
import sys
import tkinter as tk  

sys.path.insert(0, "./src")
from project_ppchem_tools_kit.browse_excel_file import browse_excel_file

class TestBrowseExcelFile(unittest.TestCase):

    @patch('project_ppchem_tools_kit.browse_excel_file.filedialog.askopenfilename')
    @patch('project_ppchem_tools_kit.browse_excel_file.tk.Entry')
    def test_browse_excel_file(self, mock_entry, mock_askopenfilename):
        mock_askopenfilename.return_value = 'test_file.xlsx'
        
        mock_entry_input = MagicMock()
        mock_entry.return_value = mock_entry_input

        global entry_input
        entry_input = mock_entry_input

        browse_excel_file()

        mock_entry_input.delete.assert_called_once_with(0, tk.END)
        mock_entry_input.insert.assert_called_once_with(0, 'test_file.xlsx')
        mock_askopenfilename.assert_called_once_with(title="Select Excel File", filetypes=(("Excel files", "*.xlsx"), ("All files", "*.*")))

    @patch('project_ppchem_tools_kit.browse_excel_file.filedialog.askopenfilename')
    @patch('project_ppchem_tools_kit.browse_excel_file.tk.Entry')
    def test_browse_excel_file_no_selection(self, mock_entry, mock_askopenfilename):
        mock_askopenfilename.return_value = ''
        
        mock_entry_input = MagicMock()
        mock_entry.return_value = mock_entry_input

        global entry_input
        entry_input = mock_entry_input

        browse_excel_file()

        mock_entry_input.delete.assert_not_called()
        mock_entry_input.insert.assert_not_called()
        mock_askopenfilename.assert_called_once_with(title="Select Excel File", filetypes=(("Excel files", "*.xlsx"), ("All files", "*.*")))

    @patch('project_ppchem_tools_kit.browse_excel_file.filedialog.askopenfilename')
    @patch('project_ppchem_tools_kit.browse_excel_file.tk.Entry')
    def test_browse_excel_file_exception_handling(self, mock_entry, mock_askopenfilename):
        # Mock the return value of askopenfilename
        mock_askopenfilename.return_value = 'test_file.xlsx'
        
        # Create a mock entry widget and set it to raise an exception
        mock_entry_input = MagicMock()
        mock_entry_input.delete.side_effect = Exception("Test exception")
        mock_entry.return_value = mock_entry_input

        # Set the global variable entry_input to the mock entry
        global entry_input
        entry_input = mock_entry_input

        # Capture standard output
        captured_output = StringIO()
        sys.stdout = captured_output

        # Call the function to be tested
        browse_excel_file()

        # Restore standard output
        sys.stdout = sys.__stdout__

        # Check for the exception message in output
        self.assertIn("Error: Test exception", captured_output.getvalue())

if __name__ == '__main__':
    unittest.main()

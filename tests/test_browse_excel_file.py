import unittest
from unittest.mock import patch, MagicMock
import tkinter as tk

import sys
sys.path.insert(0, "./src")
from Chem_pack.browse_excel_file import browse_excel_file

class TestBrowseExcelFile(unittest.TestCase):

    def setUp(self):
        # Create a mock root window
        self.root = tk.Tk()

    @patch('tkinter.filedialog.askopenfilename', return_value='/Users/meloenzinger/Documents/GitHub/Project-ppchem-tools-kit-bis/data/test_graph_1.xlsx')
    @patch("Chem_pack.browse_excel_file.entry_input", MagicMock())
    def test_browse_excel_file(self, mock_askopenfilename):
        # Mock the entry_input widget
        entry_input = MagicMock()

        # Patch the global entry_input variable
        with patch.dict('__main__.__dict__', {'entry_input': entry_input}):
            # Call the function
            browse_excel_file()

            # Assertions
            mock_askopenfilename.assert_called_once_with(title="Select Excel File", filetypes=(("Excel files", "*.xlsx"), ("All files", "*.*")))
            entry_input.insert.assert_called_once_with(0, '/Users/meloenzinger/Documents/GitHub/Project-ppchem-tools-kit-bis/data')
            entry_input.clipboard_clear.assert_called_once()
            entry_input.clipboard_append.assert_called_once_with('/Users/meloenzinger/Documents/GitHub/Project-ppchem-tools-kit-bis/data')
            self.assertEqual(entry_input.update.call_count, 1)
            entry_input.delete.assert_called_once_with(0, tk.END)

    def tearDown(self):
        # Destroy the mock root window
        self.root.destroy()

if __name__ == '__main__':
    unittest.main()

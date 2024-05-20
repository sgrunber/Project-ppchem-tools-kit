import unittest
from unittest.mock import patch, MagicMock

import sys
sys.path.insert(0, "./src")
from project_ppchem_tools_kit.process_input import process_input


class TestProcessInput(unittest.TestCase):
    def setUp(self):
        # Patching the Entry widget
        self.patcher_entry = patch('tkinter.Entry', autospec=True)
        self.mock_entry = self.patcher_entry.start()
        self.mock_entry_instance = MagicMock()
        self.mock_entry.return_value = self.mock_entry_instance
        self.mock_entry_instance.get.return_value = ''

        # Mock the messagebox to prevent actual GUI dialog boxes during the test
        self.patcher_msgbox = patch('tkinter.messagebox', autospec=True)
        self.mock_messagebox = self.patcher_msgbox.start()

    def tearDown(self):
        self.patcher_entry.stop()
        self.patcher_msgbox.stop()

    def test_empty_input(self):
        with patch('project_ppchem_tools_kit.process_input.entry_input', self.mock_entry_instance):
            with patch('project_ppchem_tools_kit.process_input.tk.StringVar') as mock_stringvar:
                mock_radio_instance = MagicMock()
                mock_stringvar.return_value = mock_radio_instance
                mock_radio_instance.get.return_value = ''
                process_input()
                self.mock_messagebox.showerror.assert_called_once_with("Error", "Please enter a molecule name, SMILES code, or file path.")

    def test_valid_input(self):
        self.mock_entry_instance.get.return_value = 'C1=CC=CC=C1'  # Example of a valid SMILES code
        with patch('project_ppchem_tools_kit.process_input.entry_input', self.mock_entry_instance):
            with patch('project_ppchem_tools_kit.process_input.tk.StringVar') as mock_stringvar:
                mock_radio_instance = MagicMock()
                mock_stringvar.return_value = mock_radio_instance
                mock_radio_instance.get.return_value = '1'  # Assuming '1' is a valid radio button selection
                process_input()
                self.mock_messagebox.showerror.assert_not_called()  # Should not show an error for valid input

    def test_invalid_radio_selection(self):
        self.mock_entry_instance.get.return_value = 'C1=CC=CC=C1'  # Example of a valid SMILES code
        with patch('project_ppchem_tools_kit.process_input.entry_input', self.mock_entry_instance):
            with patch('project_ppchem_tools_kit.process_input.tk.StringVar') as mock_stringvar:
                mock_radio_instance = MagicMock()
                mock_stringvar.return_value = mock_radio_instance
                mock_radio_instance.get.return_value = ''  # Invalid radio button selection
                process_input()
                self.mock_messagebox.showerror.assert_called_once_with("Error", "Please select a valid option.")

if __name__ == '__main__':
    unittest.main()
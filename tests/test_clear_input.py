import unittest
from unittest.mock import patch, MagicMock

import unittest
from unittest.mock import patch, MagicMock

import sys
sys.path.insert(0, "./src")

from project_ppchem_tools_kit.clear_input import clear_input

class TestClearInput(unittest.TestCase):

    @patch('project_ppchem_tools_kit.clear_input.tk.END')
    def test_clear_input(self, mock_end):
        # Create a mock entry_input widget
        entry_input_mock = MagicMock()

        # Patch the tk.Entry widget
        with patch('project_ppchem_tools_kit.clear_input.entry_input', entry_input_mock):
            # Call the clear_input function
            clear_input()

            # Assert that the delete method was called on the entry_input widget
            entry_input_mock.delete.assert_called_once_with(0, mock_end)

if __name__ == '__main__':
    unittest.main()

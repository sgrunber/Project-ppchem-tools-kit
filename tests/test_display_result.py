import unittest
from unittest.mock import MagicMock
import tkinter as tk

import sys
sys.path.insert(0, "./src")
from project_ppchem_tools_kit.display_result import display_result

class TestDisplayResult(unittest.TestCase):

    def setUp(self):
        # Create a mock root window
        self.root = tk.Tk()

    def test_display_result(self):
        # Mock the Toplevel and Text classes
        mock_toplevel = MagicMock()
        mock_text = MagicMock()

        # Patch the Toplevel and Text classes
        with unittest.mock.patch('tkinter.Toplevel', return_value=mock_toplevel):
            with unittest.mock.patch('tkinter.Text', return_value=mock_text):
                # Call the function with mock parameters
                display_result("Test result")

                # Assertions
                mock_toplevel.title.assert_called_once_with("Result")
                mock_text.insert.assert_called_once_with("1.0", "Test result")
                mock_text.config.assert_called_once_with(state="disabled")
                mock_text.grid.assert_called_once_with(row=0, column=0, padx=0, pady=0)

    def tearDown(self):
        # Destroy the mock root window
        self.root.destroy()

if __name__ == '__main__':
    unittest.main()

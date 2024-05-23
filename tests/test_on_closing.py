import unittest
from unittest.mock import MagicMock, patch
from tkinter import messagebox

import sys

from project_ppchem_tools_kit.on_closing import on_closing

class TestOnClosing(unittest.TestCase):

    def setUp(self):
        # Create a mock window
        self.window = MagicMock()
        self.window.destroy = MagicMock()

    @patch('project_ppchem_tools_kit.on_closing.messagebox.askokcancel', return_value=True)
    def test_on_closing_quit(self, mock_askokcancel):
        # Call the on_closing function
        on_closing(self.window)

        # Assert that messagebox.askokcancel was called with the correct parameters
        mock_askokcancel.assert_called_once_with("Quit", "Are you sure you want to quit?")

        # Assert that the window's destroy method was called
        self.window.destroy.assert_called_once()

    @patch('project_ppchem_tools_kit.on_closing.messagebox.askokcancel', return_value=False)
    def test_on_closing_cancel(self, mock_askokcancel):
        # Call the on_closing function
        on_closing(self.window)

        # Assert that messagebox.askokcancel was called with the correct parameters
        mock_askokcancel.assert_called_once_with("Quit", "Are you sure you want to quit?")

        # Assert that the window's destroy method was not called
        self.window.destroy.assert_not_called()

if __name__ == '__main__':
    unittest.main()

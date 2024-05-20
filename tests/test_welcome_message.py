import unittest
from unittest.mock import patch, MagicMock

import sys
sys.path.insert(0, "./src")
from project_ppchem_tools_kit.welcome_message import welcome_message

class TestWelcomeMessage(unittest.TestCase):

    @patch('project_ppchem_tools_kit.welcome_message.tk.Toplevel', return_value=None)
    @patch('project_ppchem_tools_kit.welcome_message.tk.Text')
    def test_welcome_message(self, mock_text, mock_toplevel):
        mock_toplevel_instance = MagicMock()
        mock_toplevel.return_value = mock_toplevel_instance

        welcome_message.welcome_message()

        # Assert that Toplevel and Text were called with the correct parameters
        mock_toplevel.assert_called_once_with(welcome_message.window)
        mock_text.assert_called_once_with(mock_toplevel_instance, wrap="word", font=("Times New Roman", 25), fg="#BBE1FA", bg="#1B262C", height=7, width=45)

        # Assert that the title method was called on the Toplevel instance
        mock_toplevel_instance.title.assert_called_once_with("Welcome Message")

if __name__ == '__main__':
    unittest.main()
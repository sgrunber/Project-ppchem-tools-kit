import unittest
from unittest.mock import MagicMock, patch

import sys
sys.path.insert(0, "./src")
from project_ppchem_tools_kit.copy_text import copy_text
from io import StringIO

class TestCopyText(unittest.TestCase):
    @patch('project_ppchem_tools_kit.copy_text.window')
    def test_copy_text_with_valid_event_and_window(self, mock_window):
        # Mock window and event
        mock_window.winfo_exists.return_value = True
        
        mock_event = MagicMock()
        mock_event.widget.winfo_exists.return_value = True
        
        result = copy_text(mock_event)
        
        mock_event.widget.event_generate.assert_called_once_with("<<Copy>>")
        self.assertEqual(result, "break")

    @patch('project_ppchem_tools_kit.copy_text.window')
    def test_copy_text_with_no_window(self, mock_window):
        # Set window to None
        mock_window.winfo_exists.return_value = False
        
        mock_event = MagicMock()
        mock_event.widget.winfo_exists.return_value = True
        
        result = copy_text(mock_event)
        
        mock_event.widget.event_generate.assert_not_called()
        self.assertEqual(result, "break")

    def test_copy_text_with_no_event(self):
        result = copy_text(None)
        
        self.assertEqual(result, "break")

    @patch('project_ppchem_tools_kit.copy_text.window')
    def test_copy_text_with_no_widget_in_event(self, mock_window):
        # Mock window and event without widget
        mock_window.winfo_exists.return_value = True
        
        mock_event = MagicMock()
        mock_event.widget = None
        
        result = copy_text(mock_event)
        
        self.assertEqual(result, "break")

    @patch('project_ppchem_tools_kit.copy_text.window')
    def test_copy_text_with_event_widget_not_existing(self, mock_window):
        # Mock window and event with widget not existing
        mock_window.winfo_exists.return_value = True
        
        mock_event = MagicMock()
        mock_event.widget.winfo_exists.return_value = False
        
        result = copy_text(mock_event)
        
        mock_event.widget.event_generate.assert_not_called()
        self.assertEqual(result, "break")

    @patch('project_ppchem_tools_kit.copy_text.window')
    def test_copy_text_raises_exception(self, mock_window):
        # Mock window and event with an exception being raised
        mock_window.winfo_exists.side_effect = Exception("Test Exception")
        
        mock_event = MagicMock()

        # Redirect stdout to capture print statements
        captured_output = StringIO()
        sys.stdout = captured_output

        result = copy_text(mock_event)
        
        sys.stdout = sys.__stdout__  # Reset redirect.

        self.assertIn("Error while copying text: Test Exception", captured_output.getvalue())
        self.assertEqual(result, "break")

if __name__ == '__main__':
    unittest.main()
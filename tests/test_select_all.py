import unittest
from unittest.mock import MagicMock, patch

import sys

from project_ppchem_tools_kit.select_all import select_all

class TestSelectAll(unittest.TestCase):
    
    @patch('project_ppchem_tools_kit.select_all.window')
    def test_select_all_window_exists(self, mock_window):
        mock_window.winfo_exists.return_value = True
        
        event = MagicMock()
        event.widget = MagicMock()
        
        result = select_all(event)
        
        event.widget.tag_add.assert_called_once_with("sel", "1.0", "end")
        self.assertEqual(result, "break")

    @patch('project_ppchem_tools_kit.select_all.window')
    def test_select_all_window_not_exists(self, mock_window):
        mock_window.winfo_exists.return_value = False
        
        event = MagicMock()
        event.widget = MagicMock()
        
        result = select_all(event)
        
        event.widget.tag_add.assert_not_called()
        self.assertEqual(result, "break")

    def test_select_all_window_none(self):
        with patch('project_ppchem_tools_kit.select_all.window', None):
            event = MagicMock()
            event.widget = MagicMock()
            
            result = select_all(event)
            
            event.widget.tag_add.assert_not_called()
            self.assertEqual(result, "break")

if __name__ == '__main__':
    unittest.main()
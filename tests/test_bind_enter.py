import unittest
from unittest.mock import patch, MagicMock
from io import StringIO


import sys
sys.path.insert(0, "./src")
from Chem_pack.bind_enter import bind_enter
from Chem_pack.process_input import process_input


class TestBindEnter(unittest.TestCase):
    @patch('Chem_pack.bind_enter.window')
    @patch('Chem_pack.process_input')
    def test_bind_enter_with_return_key(self, mock_process_input, mock_window):
        # Mock event with keysym "Return"
        mock_event = MagicMock()
        mock_event.keysym = "Return"
        
        result = bind_enter(mock_event)
        
        mock_process_input.assert_called_once()
        self.assertIsNone(result)

    @patch('Chem_pack.bind_enter.window')
    @patch('Chem_pack.process_input')
    def test_bind_enter_with_non_return_key(self, mock_process_input, mock_window):
        # Mock event with keysym other than "Return"
        mock_event = MagicMock()
        mock_event.keysym = "Escape"
        
        result = bind_enter(mock_event)
        
        mock_process_input.assert_not_called()
        self.assertIsNone(result)

    @patch('Chem_pack.bind_enter.window')
    @patch('Chem_pack.process_input')
    def test_bind_enter_raises_exception(self, mock_process_input, mock_window):
        # Mock event with keysym "Return" and process_input raising an exception
        mock_event = MagicMock()
        mock_event.keysym = "Return"
        mock_process_input.side_effect = Exception("Test Exception")
        
        # Redirect stdout to capture print statements
        captured_output = StringIO()
        sys.stdout = captured_output

        result = bind_enter(mock_event)
        
        sys.stdout = sys.__stdout__  # Reset redirect.

        self.assertIn("Error while binding Enter key: Test Exception", captured_output.getvalue())
        self.assertIsNone(result)

if __name__ == '__main__':
    unittest.main()
import unittest
from unittest.mock import patch, MagicMock

import sys
sys.path.insert(0, "./src")
from Chem_pack.bind_enter import bind_enter
from Chem_pack.process_input import process_input

class TestBindEnter(unittest.TestCase):

    @patch('Chem_pack.process_input')
    def test_bind_enter_with_enter_key(self, mock_process_input):
        # Create a mock event object with keysym "Return" to simulate Enter key
        event = MagicMock()
        event.keysym = "Return"

        # Call the function with the mock event
        bind_enter(event)

        # Assert that process_input was called once
        mock_process_input.assert_called_once()

    @patch('builtins.print')
    def test_bind_enter_with_exception(self, mock_print):
        # Create a mock event object with keysym "Return" to simulate Enter key
        event = MagicMock()
        event.keysym = "Return"

        # Mock process_input to raise an exception
        with patch('Chem_pack.process_input', side_effect=Exception("Test exception")):
            # Call the function with the mock event
            bind_enter(event)

        # Assert that the error message is printed correctly
        mock_print.assert_called_once_with("Error while binding Enter key:", "Test exception")

if __name__ == '__main__':
    unittest.main()

import unittest
from unittest.mock import patch
import tkinter as tk

import sys
sys.path.insert(0, "./src")
from Chem_pack.clear_input import clear_input
from Chem_pack.on_radio_select import on_radio_select  # Import the function you want to test

class TestOnRadioSelect(unittest.TestCase):
    @patch('Chem_pack.clear_input')  # Patch the clear_input function
    def test_on_radio_select(self, mock_clear_input):
        # Initialize selected_radio as a tk.StringVar()
        selected_radio = tk.StringVar()

        value = "selected_value"
        on_radio_select(value)

        # Assert that selected_radio is updated with the correct value
        self.assertEqual(selected_radio.get(), value)

        # Assert that clear_input is called once
        mock_clear_input.assert_called_once()

if __name__ == '__main__':
    unittest.main()
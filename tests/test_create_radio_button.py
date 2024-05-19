import unittest
from unittest.mock import MagicMock
from tkinter import Canvas

import sys
sys.path.insert(0, "./src")
from project_ppchem_tools_kit.create_radio_button import create_radio_button

class TestCreateRadioButton(unittest.TestCase):

    def setUp(self):
        # Create a mock canvas
        self.canvas = MagicMock(spec=Canvas)
        # Create a mock selected_radio variable
        self.selected_radio = MagicMock()

    def test_create_radio_button(self):
        # Define the parameters for the radio button
        x = 100
        y = 100
        text = "Option 1"
        value = "option1"

        # Call the function to create the radio button
        create_radio_button(x, y, text, value)

        # Assert that the canvas's create_window method was called with the correct parameters
        self.canvas.create_window.assert_called_once()

        # Assert that the Radiobutton was created with the correct parameters
        self.canvas.create_window.return_value.Radiobutton.assert_called_once()

        # Assert that the command function is correctly set to on_radio_select
        self.canvas.create_window.return_value.Radiobutton.return_value.config.assert_called_once()

if __name__ == '__main__':
    unittest.main()
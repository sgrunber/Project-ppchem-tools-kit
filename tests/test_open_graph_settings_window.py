import unittest
from unittest.mock import patch, MagicMock
import tkinter as tk

import sys
sys.path.insert(0, "./src")
from Chem_pack.open_graph_settings_window import open_graph_settings_window

class TestOpenGraphSettingsWindow(unittest.TestCase):

    @patch('Chem_pack.set_axes_color')
    @patch('Chem_pack.set_grid_color')
    @patch('Chem_pack.set_label_color')
    @patch('Chem_pack.set_background_color')
    @patch('Chem_pack.tk.Button')
    def test_open_graph_settings_window(self, mock_button, mock_set_background_color, mock_set_label_color, mock_set_grid_color, mock_set_axes_color):
        # Create a mock figure and plot canvas
        fig = MagicMock()
        plot_canvas = MagicMock()

        # Create a mock Tkinter Toplevel window
        graph_settings_window = MagicMock(spec=tk.Toplevel)
        graph_settings_window.Button = mock_button
        tk.Toplevel = MagicMock(return_value=graph_settings_window)

        # Call the function
        open_graph_settings_window(fig, plot_canvas)

        # Assert that Tkinter Toplevel was created
        tk.Toplevel.assert_called_once()

        # Assert that buttons were created and their commands were bound correctly
        mock_button.assert_any_call(graph_settings_window, text="Set Line Color", command=mock_set_axes_color)
        mock_button.assert_any_call(graph_settings_window, text="Set Axis Color", command=mock_set_axes_color)
        mock_button.assert_any_call(graph_settings_window, text="Set grid", command=mock_set_grid_color)
        mock_button.assert_any_call(graph_settings_window, text="Set Label Color", command=mock_set_label_color)
        mock_button.assert_any_call(graph_settings_window, text="Set Background Color", command=mock_set_background_color)
        mock_button.assert_any_call(graph_settings_window, text="Close", command=graph_settings_window.destroy)

if __name__ == '__main__':
    unittest.main()

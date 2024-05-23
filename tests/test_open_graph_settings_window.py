import unittest
from unittest.mock import patch, MagicMock, call
import tkinter as tk

import sys

from project_ppchem_tools_kit.open_graph_settings_window import open_graph_settings_window


class TestOpenGraphSettingsWindow(unittest.TestCase):

    @patch('project_ppchem_tools_kit.open_graph_settings_window.tk.Button')
    @patch('project_ppchem_tools_kit.open_graph_settings_window.tk.Toplevel')
    @patch('project_ppchem_tools_kit.open_graph_settings_window.set_axes_color')
    @patch('project_ppchem_tools_kit.open_graph_settings_window.set_grid_color')
    @patch('project_ppchem_tools_kit.open_graph_settings_window.set_label_color')
    @patch('project_ppchem_tools_kit.open_graph_settings_window.set_background_color')
    def test_open_graph_settings_window(self, mock_set_background_color, mock_set_label_color, mock_set_grid_color, mock_set_axes_color, mock_toplevel, mock_button):
        # Mocking the Toplevel window
        mock_toplevel_instance = MagicMock()
        mock_toplevel.return_value = mock_toplevel_instance
        
        # Mocking the figure and plot_canvas
        mock_fig = MagicMock()
        mock_plot_canvas = MagicMock()

        # Call the function to create the settings window
        open_graph_settings_window(mock_fig, mock_plot_canvas)

        # Ensure Toplevel window was created
        mock_toplevel.assert_called_once()
        mock_toplevel_instance.title.assert_called_once_with("Graph Settings")

        # Ensure the correct number of buttons were created
        self.assertEqual(mock_button.call_count, 6)

        # Find all Button instances and their commands
        button_calls = mock_button.call_args_list
        button_commands = [call[1]['command'] for call in button_calls]

        # Simulate button clicks
        button_commands[0]()
        mock_set_axes_color.assert_called_with(mock_fig.axes[0], mock_plot_canvas)

        button_commands[1]()
        mock_set_axes_color.assert_called_with(mock_fig.axes[0], mock_plot_canvas)

        button_commands[2]()
        mock_set_grid_color.assert_called_with(mock_fig.axes[0], mock_plot_canvas)

        button_commands[3]()
        mock_set_label_color.assert_called_with(mock_fig.axes[0], mock_plot_canvas)

        button_commands[4]()
        mock_set_background_color.assert_called_with(mock_fig.axes[0], mock_plot_canvas)

        # Simulate the close button click
        button_commands[5]()
        mock_toplevel_instance.destroy.assert_called_once()

if __name__ == '__main__':
    unittest.main()


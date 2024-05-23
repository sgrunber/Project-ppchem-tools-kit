import sys

from project_ppchem_tools_kit.display_graph import display_graph 

# test_display_graph.py
import unittest
from unittest.mock import patch, MagicMock
import tkinter as tk
from matplotlib.figure import Figure
from project_ppchem_tools_kit.display_graph import display_graph, graph_window

class TestDisplayGraph(unittest.TestCase):

    @patch('project_ppchem_tools_kit.display_graph.tk.Toplevel')
    @patch('project_ppchem_tools_kit.display_graph.FigureCanvasTkAgg')
    @patch('project_ppchem_tools_kit.display_graph.NavigationToolbar2Tk')
    @patch('project_ppchem_tools_kit.display_graph.set_custom_labels_and_title')
    @patch('project_ppchem_tools_kit.display_graph.add_data_point')
    @patch('project_ppchem_tools_kit.display_graph.set_scale')
    @patch('project_ppchem_tools_kit.display_graph.toggle_grid')
    @patch('project_ppchem_tools_kit.display_graph.display_max_point_and_coords')
    @patch('project_ppchem_tools_kit.display_graph.open_graph_settings_window')
    @patch('project_ppchem_tools_kit.display_graph.messagebox.askokcancel', return_value=True)
    @patch('project_ppchem_tools_kit.display_graph.tk.Button')
    def test_display_graph(self, mock_Button, mock_askokcancel, mock_open_graph_settings_window, mock_display_max_point_and_coords,
                           mock_toggle_grid, mock_set_scale, mock_add_data_point, mock_set_custom_labels_and_title,
                           mock_NavigationToolbar2Tk, mock_FigureCanvasTkAgg, mock_Toplevel):
        # Mocking the Toplevel window
        mock_toplevel_instance = MagicMock()
        mock_Toplevel.return_value = mock_toplevel_instance

        # Mocking FigureCanvasTkAgg and NavigationToolbar2Tk
        mock_canvas_instance = MagicMock()
        mock_FigureCanvasTkAgg.return_value = mock_canvas_instance
        mock_toolbar_instance = MagicMock()
        mock_NavigationToolbar2Tk.return_value = mock_toolbar_instance

        # Creating a mock figure
        mock_fig = MagicMock(spec=Figure)

        # Ensure graph_window is reset before the test
        global graph_window
        graph_window = None

        # Call the function to display the graph
        display_graph(mock_fig)

        # Ensure Toplevel window was created
        mock_Toplevel.assert_called_once()
        mock_toplevel_instance.title.assert_called_once_with("Graph")
        mock_toplevel_instance.geometry.assert_called_once_with("1000x600")

        # Ensure FigureCanvasTkAgg was created and packed
        mock_FigureCanvasTkAgg.assert_called_once_with(mock_fig, master=mock_toplevel_instance)
        mock_canvas_instance.draw.assert_called_once()
        mock_canvas_instance.get_tk_widget().pack.assert_called_once_with(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Ensure NavigationToolbar2Tk was created and updated
        mock_NavigationToolbar2Tk.assert_called_once_with(mock_canvas_instance, mock_toplevel_instance)
        mock_toolbar_instance.update.assert_called_once()

        # Ensure buttons were created and packed
        self.assertEqual(mock_Button.call_count, 6)  # Ensure 6 buttons were created

        # Simulate button clicks to verify commands
        button_calls = mock_Button.call_args_list
        button_commands = [call[1]['command'] for call in button_calls]

        button_commands[0]()
        mock_set_custom_labels_and_title.assert_called_once_with(mock_fig.axes[0], mock_canvas_instance)

        button_commands[1]()
        mock_add_data_point.assert_called_once_with(mock_fig.axes[0], mock_canvas_instance)

        button_commands[2]()
        mock_set_scale.assert_called_once_with(mock_fig.axes[0], mock_canvas_instance)

        button_commands[3]()
        mock_toggle_grid.assert_called_once_with(mock_fig.axes[0], mock_canvas_instance)

        button_commands[4]()
        mock_display_max_point_and_coords.assert_called_once_with(mock_fig.axes[0], mock_canvas_instance)

        button_commands[5]()
        mock_open_graph_settings_window.assert_called_once_with(mock_fig, mock_canvas_instance)

        # Simulate closing the window
        on_closing = mock_toplevel_instance.protocol.call_args[0][1]
        with patch('project_ppchem_tools_kit.display_graph.graph_window', mock_toplevel_instance):
            on_closing()
            mock_toplevel_instance.withdraw.assert_called_once()
            mock_toplevel_instance.quit.assert_called_once()

if __name__ == '__main__':
    unittest.main()
import unittest
from unittest.mock import MagicMock, patch, call
import tkinter as tk
from tkinter import messagebox
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import sys
sys.path.insert(0, "./src")
from Chem_pack.display_graph import display_graph  # Ensure this path is correct

class TestDisplayGraph(unittest.TestCase):

    @patch('tkinter.Toplevel')
    @patch('tkinter.Button')
    @patch('matplotlib.backends.backend_tkagg.FigureCanvasTkAgg')
    @patch('matplotlib.backends.backend_tkagg.NavigationToolbar2Tk')
    @patch('tkinter.messagebox.askokcancel', return_value=True)
    def test_display_graph(self, mock_askokcancel, mock_toolbar, mock_canvas, mock_button, mock_toplevel):
        # Set up the mock objects
        mock_window = MagicMock()
        mock_toplevel.return_value = mock_window
        mock_canvas_instance = MagicMock()
        mock_canvas.return_value = mock_canvas_instance

        fig = Figure()
        ax = fig.add_subplot(111)
        ax.plot([1, 2, 3], [4, 5, 6])

        global graph_window
        graph_window = None  # Ensure global graph_window is None to start

        # Call the function
        display_graph(fig)  # Correctly calling the function

        # Check that a new Toplevel window is created
        mock_toplevel.assert_called_once()
        mock_window.title.assert_called_once_with("Graph")
        mock_window.geometry.assert_called_once_with("1000x600")

        # Check that the FigureCanvasTkAgg and NavigationToolbar2Tk are created and set up correctly
        mock_canvas.assert_called_once_with(fig, master=mock_window)
        mock_canvas_instance.draw.assert_called_once()
        mock_canvas_instance.get_tk_widget().pack.assert_called_once_with(side=tk.TOP, fill=tk.BOTH, expand=1)
        mock_toolbar.assert_called_once_with(mock_canvas_instance, mock_window)
        mock_toolbar.return_value.update.assert_called_once()

        # Check that all buttons are created and packed correctly
        self.assertEqual(mock_button.call_count, 6)
        button_calls = [
            call(mock_window, text="Custom Labels and Title", command=unittest.mock.ANY),
            call(mock_window, text="Add Data Point", command=unittest.mock.ANY),
            call(mock_window, text="Set Scale", command=unittest.mock.ANY),
            call(mock_window, text="Toggle Grid", command=unittest.mock.ANY),
            call(mock_window, text="Display Max Point", command=unittest.mock.ANY),
            call(mock_window, text="Graph Settings", command=unittest.mock.ANY)
        ]
        mock_button.assert_has_calls(button_calls, any_order=True)

        for btn in mock_button.return_value:
            btn.pack.assert_called_once_with(side=tk.LEFT, padx=5, pady=5 if 'Graph Settings' in btn.call_args[1]['text'] else 0)

        # Check the protocol for window closing
        mock_window.protocol.assert_called_once_with("WM_DELETE_WINDOW", unittest.mock.ANY)

        # Simulate closing the window
        close_callback = mock_window.protocol.call_args[0][1]
        close_callback()
        mock_window.withdraw.assert_called_once()
        mock_window.quit.assert_called_once()

if __name__ == '__main__':
    unittest.main()
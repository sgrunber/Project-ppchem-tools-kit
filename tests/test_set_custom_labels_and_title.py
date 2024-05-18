import unittest
from unittest.mock import MagicMock, patch

import sys
sys.path.insert(0, "./src")
from Chem_pack.set_custom_labels_and_title import set_custom_labels_and_title

class TestSetCustomLabelsAndTitle(unittest.TestCase):

    @patch('tkinter.simpledialog.askstring', side_effect=["X Label", "Y Label", "Graph Title"])
    @patch('matplotlib.pyplot.tight_layout')
    def test_set_custom_labels_and_title(self, mock_tight_layout, mock_askstring):
        # Create a mock axis and plot canvas
        ax = MagicMock()
        plot_canvas = MagicMock()

        # Call the function
        set_custom_labels_and_title(ax, plot_canvas)

        # Assert that simpledialog.askstring() was called three times with correct parameters
        mock_askstring.assert_any_call("Custom Labels and Title", "Enter X-axis Label:")
        mock_askstring.assert_any_call("Custom Labels and Title", "Enter Y-axis Label:")
        mock_askstring.assert_any_call("Custom Labels and Title", "Enter Graph Title:")

        # Assert that ax.set_xlabel(), ax.set_ylabel(), ax.set_title(), and plot_canvas.draw_idle() were called
        ax.set_xlabel.assert_called_once_with("X Label", fontsize=20)
        ax.set_ylabel.assert_called_once_with("Y Label", fontsize=20)
        ax.set_title.assert_called_once_with("Graph Title", fontsize=25, fontweight='bold')
        mock_tight_layout.assert_called_once()
        plot_canvas.draw_idle.assert_called_once()

if __name__ == '__main__':
    unittest.main()
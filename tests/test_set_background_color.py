import unittest
from unittest.mock import patch, MagicMock
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "./src")
from project_ppchem_tools_kit.set_background_color import set_background_color

class TestSetBackgroundColor(unittest.TestCase):

    @patch('Chem_pack.set_background_color.colorchooser.askcolor', return_value=((255, 0, 0), '#ff0000'))
    def test_set_background_color(self, mock_askcolor):
        # Create a mock figure and axis
        fig = MagicMock()
        ax = MagicMock()
        ax.figure = fig
        fig.get_axes = MagicMock(return_value=[ax])
        plot_canvas = MagicMock()
        plot_canvas.get_tk_widget.return_value.configure.return_value = None

        # Call the function
        set_background_color(ax, plot_canvas)

        # Assert that the color chooser was called
        mock_askcolor.assert_called_once()

        # Assert that background colors are set correctly
        fig.set_facecolor.assert_called_once_with('#ff0000')
        plot_canvas.get_tk_widget.return_value.configure.assert_called_once_with(bg='#ff0000')
        ax.set_facecolor.assert_called_once_with('#ff0000')

        # Assert that plot_canvas.draw_idle() was called
        plot_canvas.draw_idle.assert_called_once()

if __name__ == '__main__':
    unittest.main()

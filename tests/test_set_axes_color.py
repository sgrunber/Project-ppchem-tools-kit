import unittest
from unittest.mock import patch, MagicMock
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "./src")
from project_ppchem_tools_kit.set_axes_color import set_axes_color

class TestSetAxesColor(unittest.TestCase):

    @patch('Chem_pack.set_axes_color.colorchooser.askcolor', return_value=((255, 0, 0), '#ff0000'))
    def test_set_axes_color(self, mock_askcolor):
        # Create a mock axis and plot canvas
        ax = MagicMock()
        ax.spines = {'left': MagicMock(), 'right': MagicMock(), 'bottom': MagicMock(), 'top': MagicMock()}
        plot_canvas = MagicMock()

        # Call the function
        set_axes_color(ax, plot_canvas)

        # Assert that the color chooser was called
        mock_askcolor.assert_called_once()

        # Assert that spine colors are set correctly
        for spine in ax.spines.values():
            spine.set_color.assert_called_once_with('#ff0000')

        # Assert that tick colors are set correctly
        ax.tick_params.assert_any_call(axis='x', colors='#ff0000')
        ax.tick_params.assert_any_call(axis='y', colors='#ff0000')

        # Assert that plot_canvas.draw_idle() was called
        plot_canvas.draw_idle.assert_called_once()

if __name__ == '__main__':
    unittest.main()
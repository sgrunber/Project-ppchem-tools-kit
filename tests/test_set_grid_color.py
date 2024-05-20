import unittest
from unittest.mock import patch, MagicMock
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "./src")
from project_ppchem_tools_kit.set_grid_color import set_grid_color

class TestSetGridColor(unittest.TestCase):

    @patch('project_ppchem_tools_kit.set_grid_color.colorchooser.askcolor', return_value=((255, 0, 0), '#ff0000'))
    def test_set_grid_color(self, mock_askcolor):
        # Create a mock axis and plot canvas
        ax = MagicMock(spec=plt.Axes)
        plot_canvas = MagicMock()

        # Call the function
        set_grid_color(ax, plot_canvas)

        # Assert that the color chooser was called
        mock_askcolor.assert_called_once()

        # Assert that grid color is set correctly
        ax.grid.assert_called_once_with(color='#ff0000')

        # Assert that plot_canvas.draw_idle() was called
        plot_canvas.draw_idle.assert_called_once()

if __name__ == '__main__':
    unittest.main()

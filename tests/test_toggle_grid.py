import unittest
from unittest.mock import MagicMock

import sys
sys.path.insert(0, "./src")
from project_ppchem_tools_kit.toggle_grid import toggle_grid

class TestToggleGrid(unittest.TestCase):

    def test_toggle_grid(self):
        # Create a mock axis and plot canvas
        ax = MagicMock()
        plot_canvas = MagicMock()

        # Call the function
        toggle_grid(ax, plot_canvas)

        # Assert that grid lines visibility is set to True for the axis
        ax.grid.assert_called_once_with(True)

        # Assert that plot_canvas.draw_idle() was called
        plot_canvas.draw_idle.assert_called_once()

if __name__ == '__main__':
    unittest.main()

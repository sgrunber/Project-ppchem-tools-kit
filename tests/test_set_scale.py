import unittest
from unittest.mock import patch, MagicMock, call
from matplotlib.ticker import MultipleLocator

import sys
sys.path.insert(0, "./src")
from Chem_pack.set_scale import set_scale

class TestSetScale(unittest.TestCase):

    @patch('Chem_pack.set_scale.simpledialog.askfloat', side_effect=[2.0, 3.0])  # Mock user input for x and y spacings
    def test_set_scale(self, mock_askfloat):
        # Create a mock axis and plot canvas
        ax = MagicMock()
        ax.figure = MagicMock()
        plot_canvas = MagicMock()

        # Call the function
        set_scale(ax, plot_canvas)

        # Assert that simpledialog.askfloat() was called twice with correct parameters
        mock_askfloat.assert_any_call("X Spacing", "Enter spacing between x ticks:")
        mock_askfloat.assert_any_call("Y Spacing", "Enter spacing between y ticks:")

        # Capture the actual calls made to ax.xaxis.set_major_locator
        actual_calls = ax.xaxis.set_major_locator.call_args

        # Assert that major locators for x and y axes are set correctly
        expected_calls = call(MultipleLocator(2.0))
        actual_calls.assert_called_once_with(expected_calls)

        # Assert that plt.tight_layout() and plot_canvas.draw_idle() were called
        ax.figure.tight_layout.assert_called_once()
        plot_canvas.draw_idle.assert_called_once()

if __name__ == '__main__':
    unittest.main()
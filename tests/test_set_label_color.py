import unittest
from unittest.mock import patch, MagicMock
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "./src")
from Chem_pack.set_label_color import set_label_color


class TestSetLabelColor(unittest.TestCase):

    @patch('Chem_pack.set_label_color.colorchooser.askcolor', return_value=((255, 0, 0), '#ff0000'))
    def test_set_label_color(self, mock_askcolor):
        # Create a mock axis and plot canvas
        ax = MagicMock(spec=plt.Axes)
        ax.xaxis = MagicMock()
        ax.yaxis = MagicMock()
        ax.title = MagicMock()
        plot_canvas = MagicMock()

        # Call the function
        set_label_color(ax, plot_canvas)

        # Assert that the color chooser was called
        mock_askcolor.assert_called_once()

        # Assert that label colors are set correctly
        ax.xaxis.label.set_color.assert_called_once_with('#ff0000')
        ax.yaxis.label.set_color.assert_called_once_with('#ff0000')
        ax.title.set_color.assert_called_once_with('#ff0000')

        # Assert that plot_canvas.draw_idle() was called
        plot_canvas.draw_idle.assert_called_once()

if __name__ == '__main__':
    unittest.main()

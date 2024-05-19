import unittest
from unittest.mock import patch, MagicMock
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "./src")
from project_ppchem_tools_kit.set_line_color import set_line_color

class TestSetLineColor(unittest.TestCase):
    def setUp(self):
        # Create a mock Axes object
        self.fig, self.ax = plt.subplots()
        x = np.linspace(0, 10, 100)
        y = np.sin(x)
        self.ax.plot(x, y)

        # Create a mock plot_canvas object
        self.plot_canvas = MagicMock()

        # Patch the colorchooser and simpledialog from tkinter
        self.patcher_colorchooser = patch('tkinter.colorchooser.askcolor', autospec=True)
        self.mock_colorchooser = self.patcher_colorchooser.start()

        self.patcher_simpledialog = patch('tkinter.simpledialog.askstring', autospec=True)
        self.mock_simpledialog = self.patcher_simpledialog.start()

    def tearDown(self):
        self.patcher_colorchooser.stop()
        self.patcher_simpledialog.stop()

    def test_set_line_color_with_valid_inputs(self):
        # Mock the colorchooser to return a valid color
        self.mock_colorchooser.return_value = ((255, 0, 0), '#ff0000')
        # Mock the simpledialog to return a valid line style
        self.mock_simpledialog.return_value = '--'

        set_line_color(self.ax, self.plot_canvas)

        # Check if the line color is changed
        self.assertEqual(self.ax.lines[0].get_color(), '#ff0000')
        # Check if the line style is changed
        self.assertEqual(self.ax.lines[0].get_linestyle(), '--')
        # Check if the plot canvas is updated
        self.plot_canvas.draw_idle.assert_called_once()

    def test_set_line_color_cancel_color(self):
        # Mock the colorchooser to return None (simulate cancel action)
        self.mock_colorchooser.return_value = (None, None)
        # Mock the simpledialog to return a valid line style
        self.mock_simpledialog.return_value = '--'

        set_line_color(self.ax, self.plot_canvas)

        # Ensure that the line color is not changed
        self.assertNotEqual(self.ax.lines[0].get_color(), '#ff0000')
        # Ensure that the plot canvas is not updated
        self.plot_canvas.draw_idle.assert_not_called()

    def test_set_line_color_cancel_line_style(self):
        # Mock the colorchooser to return a valid color
        self.mock_colorchooser.return_value = ((255, 0, 0), '#ff0000')
        # Mock the simpledialog to return None (simulate cancel action)
        self.mock_simpledialog.return_value = None

        set_line_color(self.ax, self.plot_canvas)

        # Ensure that the line color is not changed
        self.assertNotEqual(self.ax.lines[0].get_color(), '#ff0000')
        # Ensure that the plot canvas is not updated
        self.plot_canvas.draw_idle.assert_not_called()

if __name__ == '__main__':
    unittest.main()

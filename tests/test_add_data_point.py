import unittest
from unittest.mock import patch, MagicMock
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "./src")
from project_ppchem_tools_kit.add_data_point import add_data_point

class TestAddDataPoint(unittest.TestCase):
    def setUp(self):
        # Create a mock Axes object
        self.fig, self.ax = plt.subplots()
        x = np.linspace(0, 10, 100)
        y = np.sin(x)
        self.ax.plot(x, y)

        # Create a mock plot_canvas object
        self.plot_canvas = MagicMock()

        # Patch the colorchooser from tkinter
        self.patcher_colorchooser = patch('tkinter.colorchooser.askcolor', autospec=True)
        self.mock_colorchooser = self.patcher_colorchooser.start()

    def tearDown(self):
        self.patcher_colorchooser.stop()

    def test_add_data_point_with_color(self):
        # Mock the colorchooser to return a valid color
        self.mock_colorchooser.return_value = ((255, 0, 0), '#ff0000')

        add_data_point(self.ax, self.plot_canvas)

        # Check if the scatter plot is added to the axes
        self.assertEqual(len(self.ax.collections), 1)

        # Check if the plot canvas is updated
        self.plot_canvas.draw_idle.assert_called_once()

    def test_add_data_point_cancel_color(self):
        # Mock the colorchooser to return None (simulate cancel action)
        self.mock_colorchooser.return_value = (None, None)

        add_data_point(self.ax, self.plot_canvas)

        # Ensure that no scatter plot is added to the axes
        self.assertEqual(len(self.ax.collections), 0)

        # Ensure that the plot canvas is not updated
        self.plot_canvas.draw_idle.assert_not_called()

if __name__ == '__main__':
    unittest.main()

import unittest
from unittest.mock import patch, MagicMock
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "./src")
from project_ppchem_tools_kit.display_max_point_and_coords import display_max_point_and_coords

class TestDisplayMaxPointAndCoords(unittest.TestCase):
    def setUp(self):
        # Create a mock Axes object
        self.fig, self.ax = plt.subplots()
        x = np.linspace(0, 10, 100)
        y = np.sin(x)
        self.ax.plot(x, y)

        # Create a mock plot_canvas object
        self.plot_canvas = MagicMock()

        # Patch the messagebox and colorchooser from tkinter
        self.patcher_messagebox = patch('tkinter.messagebox.showwarning', autospec=True)
        self.mock_messagebox = self.patcher_messagebox.start()

        self.patcher_colorchooser = patch('tkinter.colorchooser.askcolor', autospec=True)
        self.mock_colorchooser = self.patcher_colorchooser.start()

    def tearDown(self):
        self.patcher_messagebox.stop()
        self.patcher_colorchooser.stop()

    def test_display_max_point_in_range(self):
        # Mock the colorchooser to return a valid color
        self.mock_colorchooser.return_value = ((255, 0, 0), '#ff0000')

        display_max_point_and_coords(self.ax, self.plot_canvas)

        # Check if the point is plotted on the axes
        self.assertEqual(len(self.ax.collections), 1)

        # Check if the text annotation is added to the axes
        self.assertEqual(len(self.ax.texts), 1)
        
        # Check if the plot canvas is updated
        self.plot_canvas.draw_idle.assert_called_once()

    def test_display_max_point_out_of_range(self):
        # Modify the data to make the max point out of range
        self.ax.set_xlim(5, 10)
        self.ax.set_ylim(-0.5, 0.5)

        # Mock the colorchooser to return a valid color
        self.mock_colorchooser.return_value = ((255, 0, 0), '#ff0000')

        display_max_point_and_coords(self.ax, self.plot_canvas)

        # Check if the warning messagebox is shown
        self.mock_messagebox.assert_called_once_with("Warning", "Max point coordinates are out of range.")

    def test_choose_color_cancel(self):
        # Mock the colorchooser to return None (simulate cancel action)
        self.mock_colorchooser.return_value = (None, None)

        display_max_point_and_coords(self.ax, self.plot_canvas)

        # Ensure that no point is plotted on the axes
        self.assertEqual(len(self.ax.collections), 0)

        # Ensure that no text annotation is added to the axes
        self.assertEqual(len(self.ax.texts), 0)
        
        # Ensure that the plot canvas is not updated
        self.plot_canvas.draw_idle.assert_not_called()

if __name__ == '__main__':
    unittest.main()

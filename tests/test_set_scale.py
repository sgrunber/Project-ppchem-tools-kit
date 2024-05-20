import unittest
from unittest.mock import patch, MagicMock

import sys
sys.path.insert(0, "./src")
from project_ppchem_tools_kit.set_scale import set_scale


class TestSetScale(unittest.TestCase):

    @patch('project_ppchem_tools_kit.set_scale.simpledialog.askfloat')
    @patch('project_ppchem_tools_kit.set_scale.plt.MultipleLocator')
    def test_set_scale(self, mock_MultipleLocator, mock_askfloat):
        # Mocking the return values for askfloat
        mock_askfloat.side_effect = [1.0, 2.0]
        
        # Creating mock objects for ax and plot_canvas
        mock_ax = MagicMock()
        mock_plot_canvas = MagicMock()

        # Call the function
        set_scale(mock_ax, mock_plot_canvas)

        # Verify askfloat was called twice
        self.assertEqual(mock_askfloat.call_count, 2)
        
        # Verify MultipleLocator was called with the correct arguments
        mock_MultipleLocator.assert_any_call(1.0)
        mock_MultipleLocator.assert_any_call(2.0)
        
        # Verify ax.xaxis.set_major_locator and ax.yaxis.set_major_locator were called
        self.assertEqual(mock_ax.xaxis.set_major_locator.call_count, 1)
        self.assertEqual(mock_ax.yaxis.set_major_locator.call_count, 1)
        
        # Verify plot_canvas.draw_idle was called
        mock_plot_canvas.draw_idle.assert_called_once()

if __name__ == '__main__':
    unittest.main()


import unittest
from unittest.mock import MagicMock
from tkinter import Label
from datetime import datetime

import sys
sys.path.insert(0, "./src")
from Chem_pack.update_clock import update_clock

class TestUpdateClock(unittest.TestCase):

    def setUp(self):
        # Create a mock window with the 'after' method
        self.window = MagicMock()
        self.window.after = MagicMock()

        # Mock the clock and date labels
        self.clock_label = MagicMock()
        self.date_label = MagicMock()

    def test_update_clock(self):
        # Call the update_clock function with the mocked window
        update_clock(self.window, self.clock_label, self.date_label)

        # Get the current time and date
        expected_time = datetime.now().strftime('%H:%M:%S')
        expected_date = datetime.now().strftime('%Y-%m-%d')

        # Assert that the clock and date labels were updated with the expected values
        self.clock_label.config.assert_called_once_with(text=expected_time)
        self.date_label.config.assert_called_once_with(text=expected_date)

        # Assert that the 'after' method was called on the window
        self.window.after.assert_called_once_with(1000, update_clock, self.window, self.clock_label, self.date_label)

if __name__ == '__main__':
    unittest.main()
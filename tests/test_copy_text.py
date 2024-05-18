import unittest
from unittest.mock import MagicMock

import sys
sys.path.insert(0, "./src")
from Chem_pack.copy_text import copy_text

class TestCopyText(unittest.TestCase):
    def test_copy_text(self):
        global window
        event = MagicMock()
        event.widget = MagicMock()
        event.widget.winfo_exists.return_value = True
        window = MagicMock()
        window.winfo_exists.return_value = True

        # Call the copy_text function
        copy_text(event)

        # Assert that event_generate method is called with "<<Copy>>"
        event.widget.event_generate.assert_called_once_with("<<Copy>>")

if __name__ == "__main__":
    unittest.main()
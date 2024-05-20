import unittest

from pathlib import Path
import os
OUTPUT_PATH = Path(os.getcwd())
ASSETS_PATH = OUTPUT_PATH / Path("../assets/frame0")

import sys
sys.path.insert(0, "./src")
from project_ppchem_tools_kit.relative_to_assets import relative_to_assets

#1.RELATIVETOASSETS
class TestRelativePath(unittest.TestCase):
    def test_relative_to_assets(self):
        # Test the path computation is correct
        test_path = "subfolder/file.txt"
        expected = ASSETS_PATH / "subfolder/file.txt"
        result = relative_to_assets(test_path)
        self.assertEqual(result, expected)
if __name__ == '__main__':
    unittest.main()
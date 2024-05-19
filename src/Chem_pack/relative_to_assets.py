import sys
from pathlib import Path
import os
OUTPUT_PATH = Path(os.getcwd())
ASSETS_PATH = OUTPUT_PATH / Path("../assets/frame0")



def relative_to_assets(path: str) -> Path:
    return ASSETS_PATH / Path(path)
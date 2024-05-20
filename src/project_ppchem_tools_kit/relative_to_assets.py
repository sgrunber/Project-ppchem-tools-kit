from pathlib import Path
import os
OUTPUT_PATH = Path(os.getcwd())
ASSETS_PATH = OUTPUT_PATH / Path("../assets/frame0")



def relative_to_assets(path: str) -> Path:
    """Constructs a path to a file located in the assets directory by combining the provided relative path with the ASSETS_PATH.


    Args:
        path (str): The relative path to be appended.

    Returns:
        Path: The full path to the asset.
    """
    return ASSETS_PATH / Path(path)


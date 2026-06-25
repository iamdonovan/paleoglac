"""Load example glacier data"""
import os
from importlib.resources import files


# can't just use paleoglac, because paleoglac is inside src/paleoglac
_SAMPLE_DIRECTORY = os.path.abspath(files('paleoglac') / '..' / '..' / 'sample_data')

_SAMPLE_DATA = {
    'ref_dem': os.path.join(_SAMPLE_DIRECTORY, 'dem.tif'),
    'outlines': os.path.join(_SAMPLE_DIRECTORY, 'outlines.gpkg'),
}

available = list(_SAMPLE_DATA.keys())

def get_path(name: str) -> str:
    """
    Get path to test/sample dataset. List of available files can be found in 'examples.available'

    :param name: Name of sample data.
    :return:
    """
    if name in _SAMPLE_DATA.keys():
        return _SAMPLE_DATA[name]
    else:
        raise ValueError(f"Data name should be one of: {list(_SAMPLE_DATA.keys())}")

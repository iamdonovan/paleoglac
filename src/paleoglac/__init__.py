from importlib.metadata import version
from paleoglac import examples
from paleoglac import paleoglac
from paleoglac import ela
from paleoglac import surface

from paleoglac.paleoglac import PaleoGlac


__version__ = version(__name__)

__all__ = [
    'examples',
    'paleoglac',
    'ela',
    'surface',
    'PaleoGlac',
    '__version__'
]

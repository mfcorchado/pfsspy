import distutils.version

import astropy

# Import pfsspy sub-modules to have them available through pfsspy.{name}
import pfsspy.coords
import pfsspy.fieldline
# Import this to register map sources
import pfsspy.map
import pfsspy.sample_data
import pfsspy.tracing
import pfsspy.utils
from .input import Input
from .output import Output
from .pfss import pfss

__all__ = ['Input', 'Output', 'pfss']


from ._version import get_versions

__version__ = get_versions()['version']
del get_versions

# Do a version check for astropy
if (distutils.version.LooseVersion(astropy.__version__) <
        distutils.version.LooseVersion("3")):
    raise RuntimeError('pfsspy requires astropy v3 to run ' +
                       f'(found version {astropy.__version__} installed)')


__citation__ = __bibtex__ = """
@article{Stansby2020,
doi = {10.21105/joss.02732},
url = {https://doi.org/10.21105/joss.02732},
year = {2020},
publisher = {The Open Journal},
volume = {5},
number = {54},
pages = {2732},
author = {David Stansby and Anthony Yeates and Samuel T. Badman},
title = {pfsspy: A Python package for potential field source surface modelling},
journal = {Journal of Open Source Software}
}
"""

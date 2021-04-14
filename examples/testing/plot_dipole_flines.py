"""
Analytic dipole field lines
===========================
"""

###############################################################################
# First, import required modules
import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import numpy as np
import sunpy.map
import pfsspy
import pfsspy.coords as coords
from pfsspy import analytic, tracing

from helpers import pffspy_output

###############################################################################
# Compare the the pfsspy solution to the analytic solutions. Cuts are taken
# on the source surface at a constant phi value to do a 1D comparison.
l = 1
m = 0
nphi = 360
ns = 180
nr = 50
rss = 2


def theta_fline(r, theta0, r0):
    """
    Analytic solution for a field line in a dipole PFSS solution.
    """
    z = r / const.R_sun / rss
    z0 = r0 / const.R_sun / rss
    return np.arccos(np.cos(theta0) *
                     np.sqrt((z0**3 + 2) /
                             (z**3 + 2) *
                             (z / z0))
                     )


###############################################################################
# Calculate PFSS solution
pfsspy_out = pffspy_output(nphi, ns, nr, rss, l, m)


###############################################################################
# Trace some field lines
phi = 180 * u.deg
theta0 = np.linspace(0, 90, 90) * u.deg
r0 = rss * const.R_sun
seeds = SkyCoord(radius=r0, lat=theta0, lon=phi,
                 frame=pfsspy_out.coordinate_frame)

tracer = tracing.FortranTracer(step_size=0.1)
flines = tracer.trace(seeds, pfsspy_out)

###
# Calculate analytical solution
thetas = [theta_fline(fline.coords.radius, theta, r0) for
          fline, theta in zip(flines, theta0)]


fig, ax = plt.subplots()
for fline, theta in zip(flines, thetas):
    ax.scatter(fline.coords.lat[-1],
               (theta[0] - fline.coords.lat[0]).to_value(u.deg))
ax.set_xlabel('Source surface seed latitude')
ax.set_ylabel('Difference in solar suface latitude (deg)')
plt.show()

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
m = 1
nphi = 360
ns = 180
nr = 50
rss = 2


def fr(r):
    rho = r / rss
    return ((rho**l * (2*l + 1)) / ((l * rho**(2*l + 1)) + l + 1))**(1 / (l + 1))


###############################################################################
# Calculate PFSS solution
pfsspy_out = pffspy_output(nphi, ns, nr, rss, l, m)


###############################################################################
# Trace some field lines
n = 10
phi = 75 * u.deg
theta = np.linspace(-90, 90, 180) * u.deg
r0 = rss * const.R_sun
seeds = SkyCoord(radius=r0, lat=theta, lon=phi,
                 frame=pfsspy_out.coordinate_frame)

tracer = tracing.FortranTracer(step_size=0.01)
flines = tracer.trace(seeds, pfsspy_out)
mask = flines.connectivities.astype(bool)
print(mask)
theta_out = flines.open_field_lines.solar_feet.lat
r_out = flines.open_field_lines.solar_feet.radius

###############################################################################
# Calculate analytical solution
theta_analytic = np.arcsin(np.sin(theta[mask]) * fr(r_out.to_value(const.R_sun)))

fig = plt.figure()
ax = fig.add_subplot()
ax.scatter(theta[mask].to_value(u.deg), np.abs(theta_out.to_value(u.deg)) -
                                               np.abs(theta_analytic.to_value(u.deg)))
plt.show()

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
from astropy.visualization import quantity_support

quantity_support()

import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
from matplotlib.gridspec import GridSpec

import pfsspy
import pfsspy.coords as coords
from helpers import fr, pffspy_output, phi_fline_coords, theta_fline_coords
from pfsspy import analytic, tracing

###############################################################################
# Compare the the pfsspy solution to the analytic solutions. Cuts are taken
# on the source surface at a constant phi value to do a 1D comparison.
l = 1
m = 1
nphi = 360
ns = 180
nr = 50
rss = 2


###############################################################################
# Calculate PFSS solution
pfsspy_out = pffspy_output(nphi, ns, nr, rss, l, m)

rss = rss * const.R_sun
###############################################################################
# Trace some field lines
n = 90
# Create 1D theta, phi arrays
phi = np.linspace(0, 360, n * 2)
phi = phi[:-1] + np.diff(phi) / 2
theta = np.arcsin(np.linspace(-1, 1, n, endpoint=False) + 1/n)
# Mesh into 2D arrays
theta, phi = np.meshgrid(theta, phi, indexing='ij')
theta, phi = theta * u.rad, phi * u.deg
seeds = SkyCoord(radius=rss, lat=theta.ravel(), lon=phi.ravel(),
                 frame=pfsspy_out.coordinate_frame)

step_size = 1
dthetas = []
print(f'Tracing {step_size}...')
# Trace
tracer = tracing.FortranTracer(step_size=step_size)
flines = tracer.trace(seeds, pfsspy_out)
# Set a mask of open field lines
mask = flines.connectivities.astype(bool).reshape(theta.shape)

# Get solar surface latitude
phi_solar = np.ones_like(phi) * np.nan
phi_solar[mask] = flines.open_field_lines.solar_feet.lon
theta_solar = np.ones_like(theta) * np.nan
theta_solar[mask] = flines.open_field_lines.solar_feet.lat
r_out = np.ones_like(theta.value) * const.R_sun * np.nan
r_out[mask] = flines.open_field_lines.solar_feet.radius

###########################################################################
# Calculate analytical solution
theta_analytic = theta_fline_coords(r_out, rss, l, m, theta)
dtheta = theta_solar - theta_analytic
phi_analytic = phi_fline_coords(r_out, rss, l, m, theta, phi)
dphi = phi_solar - phi_analytic

fig, axs = plt.subplots(nrows=2)

kwargs = dict(cmap='RdBu')
# Latitude
ax = axs[0]
im = ax.pcolormesh(phi.to_value(u.deg), np.sin(theta).value,
                   dtheta.to_value(u.deg), **kwargs)
ax.set_aspect(360 / 4)
fig.colorbar(im, aspect=10, ax=ax,
             label=r'$\theta_{pfsspy} - \theta_{analytic}$ (deg)')
ax.set_title(f'Error in solar footpoint latitude')

# Longitude
ax = axs[1]
im = ax.pcolormesh(phi.to_value(u.deg), np.sin(theta).value,
                   dphi.to_value(u.deg), **kwargs)
ax.set_aspect(360 / 4)
fig.colorbar(im, aspect=10, ax=ax,
             label=r'$\phi_{pfsspy} - \phi_{analytic}$ (deg)')
ax.set_title(f'Error in solar footpoint longitude')

plt.show()

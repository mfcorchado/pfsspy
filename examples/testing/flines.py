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

from helpers import pffspy_output, fr, theta_fline_coords

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
r0 = rss * const.R_sun
seeds = SkyCoord(radius=r0, lat=theta.ravel(), lon=phi.ravel(),
                 frame=pfsspy_out.coordinate_frame)

fig, ax = plt.subplots()
step_sizes = [10, 5, 2, 1, 0.5, 0.2]
dthetas = []
for step_size in step_sizes:
    print(f'Tracing {step_size}...')
    # Trace
    tracer = tracing.FortranTracer(step_size=step_size)
    flines = tracer.trace(seeds, pfsspy_out)
    # Set a mask of open field lines
    mask = flines.connectivities.astype(bool).reshape(theta.shape)

    # Get source surface latitude
    theta_solar = np.ones_like(theta) * np.nan
    theta_solar[mask] = flines.open_field_lines.solar_feet.lat
    r_out = np.ones_like(theta.value) * const.R_sun * np.nan
    r_out[mask] = flines.open_field_lines.solar_feet.radius

    ###########################################################################
    # Calculate analytical solution
    theta_analytic = theta_fline_coords(r_out.to_value(const.R_sun), rss, l, m, theta)
    dtheta = (theta_solar - theta_analytic).to_value(u.deg)

    fig, ax = plt.subplots()
    im = ax.pcolormesh(phi.to_value(u.deg), np.sin(theta).value, dtheta,
                       cmap='RdBu')
    ax.set_aspect(360 / 4)
    fig.colorbar(im)
    ax.set_title(f'Step size = {step_size}')

    dtheta = dtheta[np.isfinite(dtheta)]
    dthetas.append(dtheta)

fig, ax = plt.subplots()
ax.plot(step_sizes, [np.median(np.abs(dt)) for dt in dthetas],
        marker='o',
        label='Median')
ax.plot(step_sizes, [np.max(np.abs(dt)) for dt in dthetas],
        marker=11, lw=0.5, color='tab:blue',
        label='Max')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Step size')
ax.set_ylabel('$|d\theta|$ (deg)')
ax.legend()

plt.show()

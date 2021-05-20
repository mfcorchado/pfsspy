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
n = 90
# Create 1D theta, phi arrays
phi = np.linspace(0, 360, n * 2)
theta = np.arcsin(np.linspace(-1, 1, n))
# Mesh into 2D arrays
theta, phi = np.meshgrid(theta, phi, indexing='ij')
theta, phi = theta * u.rad, phi * u.deg
r0 = rss * const.R_sun
seeds = SkyCoord(radius=r0, lat=theta.ravel(), lon=phi.ravel(),
                 frame=pfsspy_out.coordinate_frame)

fig, ax = plt.subplots()
step_sizes = [1, 0.5, 0.2, 0.1, 0.05]
# step_sizes = [1, 0.5]
dthetas = []
for step_size in step_sizes:
    print(f'Tracing {step_size}...')
    # Trace
    tracer = tracing.FortranTracer(step_size=step_size)
    flines = tracer.trace(seeds, pfsspy_out)
    # Set a mask of open field lines
    mask = flines.connectivities.astype(bool).reshape(theta.shape)

    # Get source surface latitude
    theta_ss = np.ones_like(theta) * np.nan
    theta_ss[mask] = flines.open_field_lines.solar_feet.lat
    r_out = np.ones_like(theta.value) * const.R_sun * np.nan
    r_out[mask] = flines.open_field_lines.solar_feet.radius

    ###########################################################################
    # Calculate analytical solution
    theta_analytic = np.arcsin(np.sin(theta) * fr(r_out.to_value(const.R_sun)))
    dtheta = (theta_ss - theta_analytic).to_value(u.deg).ravel()
    dtheta = dtheta[np.isfinite(dtheta)]
    dthetas.append(dtheta)

for step_size, dtheta in zip(step_sizes, dthetas):
    pctiles = np.percentile(np.abs(dtheta), [1, 50, 99])
    ax.errorbar(step_size, pctiles[1], yerr=[[pctiles[1] - pctiles[0]],
                                             [pctiles[2] - pctiles[1]]],
                color='tab:blue', marker='o')

ax.set_yscale('log')
'''
im = ax.imshow((theta_ss - theta_analytic).to_value(u.deg), cmap='RdBu', extent=[0, 360, 0, 180])
ax.set_title(r'Solar surface $\theta_{traced} - \theta_{analytic}$')
fig.colorbar(im, label='deg')

print(np.nanmean(np.abs(theta_ss - theta_analytic)).to(u.deg))
'''
plt.show()

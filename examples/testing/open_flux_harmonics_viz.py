"""
Open flux
=========
Comparing total unsigned flux to analytic solutions. This is done on a fixed
number of radial grid points, and plotted as a function of spherical harmonic.
"""

###############################################################################
# First, import required modules
import json

import matplotlib.colors as mcolor
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np


with open("open_flux_harmonics.json", "r") as f:
    results = json.load(f, parse_int=int)
print(results)
###############################################################################
# Plot results
fig, ax = plt.subplots()
norm = mcolor.Normalize(vmin=1, vmax=1.06)
for lstr in results['analytic']:
    l = int(lstr)
    analytic = np.atleast_2d(list(results['analytic'][lstr].values())).T
    numeric = np.atleast_2d(list(results['numeric'][lstr].values())).T
    im = ax.pcolor([l-0.5, l+0.5], np.arange(-0.5-l, l+0.6, 1),
                   numeric / analytic,
                   norm=norm, cmap=plt.get_cmap(name='inferno', lut=12),
                   edgecolors='white', linewidths=2)
ax.set_aspect(1)
fig.colorbar(im, label=r'$\Phi_{pfsspy} / \Phi_{analytic}$')
ax.set_xlim(0.5, l + 0.5)
ax.set_ylim(-0.5 - l, l + 0.5)
ax.xaxis.set_major_locator(mticker.MultipleLocator(1))
ax.yaxis.set_major_locator(mticker.MultipleLocator(1))


def fmt(x, pos):
    return str(int(x))


ax.xaxis.set_major_formatter(fmt)
ax.yaxis.set_major_formatter(fmt)
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
ax.spines.left.set_visible(False)
ax.spines.bottom.set_visible(False)
ax.tick_params(length=0)
ax.yaxis.tick_right()
ax.set_xlabel('l')
ax.set_ylabel('m')
fig.savefig('figs/flux_harmonics.pdf', bbox_inches='tight')
plt.show()

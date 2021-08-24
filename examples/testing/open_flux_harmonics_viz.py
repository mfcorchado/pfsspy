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
for lstr in results:
    l = int(lstr)
    data = np.atleast_2d(list(results[lstr].values())).T
    im = ax.imshow(data, extent=[l-0.5, l+0.5, -0.5, l+0.5],
                   norm=norm)

fig.colorbar(im, label=r'$\Phi_{pfsspy} / \Phi_{analytic}$')
ax.set_xlim(0.5, l + 0.5)
ax.set_ylim(-0.5, l + 0.5)
ax.xaxis.set_major_locator(mticker.MultipleLocator(1))
ax.yaxis.set_major_locator(mticker.MultipleLocator(1))


def fmt(x, pos):
    return str(int(x))


ax.xaxis.set_major_formatter(fmt)
ax.yaxis.set_major_formatter(fmt)
ax.set_xlabel('l')
ax.set_ylabel('m')
fig.savefig('flux_harmonics.pdf', bbox_inches='tight')
plt.show()

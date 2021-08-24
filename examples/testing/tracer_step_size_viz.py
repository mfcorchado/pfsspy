"""
Visualising tracer step size
============================
"""

###############################################################################
# First, import required modules
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import pandas as pd
import numpy as np

dthetas = pd.read_csv('results/dthetas_1-1.csv', index_col=0)
dphis = pd.read_csv('results/dphis_1-1.csv', index_col=0)


def plot_distributions(df, ax, color):
    for pc in ([0, 100], [1, 99], [10, 90]):
        pctiles = np.nanpercentile(np.abs(df), pc, axis=1)
        ax.fill_between(df.index, pctiles[0], pctiles[1],
                        alpha=0.4, color=color, edgecolor='face')
        ax.scatter(df.index, pctiles[0], s=1, marker='.', color='black', alpha=0.5)
        ax.scatter(df.index, pctiles[1], s=1, marker='.', color='black', alpha=0.5)


fig, axs = plt.subplots(nrows=2, sharex=True, sharey=True,
                        gridspec_kw={'hspace': 0, 'wspace': 0},
                        figsize=(4, 3))
plot_distributions(dthetas, axs[0], 'tab:blue')
plot_distributions(dphis, axs[1], 'tab:orange')

ax = axs[1]
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Step size')
ax.set_ylabel('deg')
ax.set_xlim(np.min(dphis.index), np.max(dphis.index))
ax.set_ylim(1e-3, 5e1)


def formatter(x, pos):
    return str(int(x))


ax.xaxis.set_minor_formatter(mticker.FuncFormatter(formatter))
ax.xaxis.set_major_formatter(mticker.FuncFormatter(formatter))
ax.xaxis.set_major_locator(mticker.FixedLocator([2, 4, 6, 8, 10, 20, 40]))
ax.xaxis.set_minor_locator(mticker.FixedLocator([]))
fig.tight_layout()

plt.show()

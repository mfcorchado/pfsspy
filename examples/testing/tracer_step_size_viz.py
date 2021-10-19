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


def plot_distributions(axs, l, m):
    try:
        dphis = pd.read_hdf(f'results/dphis_{l}{m}.hdf', 'table')
        dthetas = pd.read_hdf(f'results/dthetas_{l}{m}.hdf', 'table')
    except FileNotFoundError:
        print(f'❌ Files not found for l={l}, m={m}')
        return
    plot_single_distribution(dphis, axs[1], 'tab:orange')
    plot_single_distribution(dthetas, axs[0], 'tab:blue')

    for ax in axs:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(np.min(dphis.index), np.max(dphis.index))
        ax.set_ylim(1e-2, 5e1)

        def formatter(x, pos):
            return str(int(x))

        ax.xaxis.set_minor_formatter(mticker.FuncFormatter(formatter))
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(formatter))
        ax.xaxis.set_major_locator(mticker.FixedLocator(
            [1, 2, 4, 8, 20]))
        ax.xaxis.set_minor_locator(mticker.FixedLocator([]))
    axs[0].set_title(f'l={l}, m={m}')
    axs[1].set_xlabel('Step size')
    axs[1].set_ylabel('deg')


def plot_single_distribution(df, ax, color):
    for pc in ([0, 100], [1, 99], [10, 90]):
        pctiles = np.nanpercentile(np.abs(df), pc, axis=1)
        ax.fill_between(df.index, 1e-3, pctiles[1],
                        alpha=0.4, color=color, edgecolor='face')
        # ax.scatter(df.index, pctiles[0], s=1, marker='.', color='black', alpha=0.5)
        ax.scatter(df.index, pctiles[1], s=1, marker='.', color='black', alpha=0.5)


ls = [1, 2, 3]
max_m = max(ls)
fig = plt.figure()
outer_grid = fig.add_gridspec(len(ls), 2 * max(ls) + 1, wspace=0, hspace=0)
for li, l in enumerate(ls):
    for mi, m in enumerate(range(-l, l+1)):
        print(l, m)
        inner_grid = outer_grid[li, m+max_m].subgridspec(2, 1, wspace=0, hspace=0)
        axs = inner_grid.subplots()

        plot_distributions(axs, l, m)

        # Axis formatting
        if l != max(ls) or m != -max(ls):
            [ax.xaxis.set_major_formatter(mticker.NullFormatter()) for ax in axs]
            [ax.yaxis.set_major_formatter(mticker.NullFormatter()) for ax in axs]
            [ax.yaxis.set_label_text('') for ax in axs]

# fig.tight_layout()
# fig.savefig(f'figs/step_size_{l}{m}.pdf', bbox_inches='tight')

plt.show()

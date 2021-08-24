"""
Open flux
=========
Comparing total unsigned flux to analytic solutions. This is done as a function
of the number of radial grid cells used in the PFSS calculation.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('results/open_flux_results.csv', index_col=0)

for lm in lms:
    l = lm // 10
    m = lm % 10
    color = {1: 'tab:blue', 2: 'tab:orange', 3: 'tab:green'}[l]
    marker = {0: 'o', 1: 10, 2: 11}[m]
    ax.plot(nrhos, df[str(lm)],
            marker=marker, color=color,
            label=f'l={l}, m={m}')

ax.legend()
ax.set_ylim(1)
ax.set_xlim(8)
ax.yaxis.grid(linestyle='--')
ax.set_xlabel('$n_{r}$')
ax.set_ylabel(r'$\Phi_{pfsspy} / \Phi_{analytic}$')
plt.show()

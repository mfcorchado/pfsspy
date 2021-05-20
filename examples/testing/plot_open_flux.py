"""
Open flux
=========

Comparing pfsspy solutions of total unsigned flux their analytic solutions.
"""

###############################################################################
# First, import required modules
import functools

import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import numpy as np
import scipy.integrate
import sunpy.map
import sympy
import pandas as pd
import pfsspy
import pfsspy.coords as coords
from pfsspy import analytic

from helpers import brss_pfsspy


def open_flux_analytic(l, m, zss):
    Br = analytic.Br(l, m, zss)
    Br = functools.partial(Br, zss)
    absBr = lambda theta, phi: np.abs(Br(theta, phi)) * np.sin(theta)
    res = scipy.integrate.nquad(absBr, ranges=[[0, np.pi], [0, 2 * np.pi]])
    return res[0]


def open_flux_numeric(l, m, zss, nrho):
    print(nrho)
    nphi = 360
    ns = 180
    br = brss_pfsspy(nphi, ns, nrho, zss, l, m)
    return np.sum(np.abs(br)) * (4 * np.pi) / nphi / ns


zss = 2

fig, ax = plt.subplots()
lms = [10, 11, 20, 21, 22, 30]

# nrhos = np.array([5, 50])
nrhos = np.arange(10, 51, 2)

df = pd.DataFrame(index=nrhos, columns=lms)

for lm in lms:
    l = lm // 10
    m = lm % 10
    flux_analytic = open_flux_analytic(l, m, zss)
    flux_numeric = []
    flux_numeric = np.array([open_flux_numeric(l, m, zss, nrho) for nrho in nrhos])
    flux_ratio = flux_numeric / flux_analytic
    df[lm] = flux_ratio

df.to_csv('results.csv')

df = pd.read_csv('results.csv', index_col=0)

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
ax.set_ylabel('$\Phi_{pfsspy} / \Phi_{analytic}$')
plt.show()

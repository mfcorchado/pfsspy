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
lm = [[1, 0], [2, 0], [3, 0]]
for l, m in lm:
    print(l, m)
    flux_analytic = open_flux_analytic(l, m, zss)
    flux_numeric = []
    nrhos = np.array([5, 10, 20, 50])
    nrhos = np.arange(10, 51, 2)
    flux_numeric = np.array([open_flux_numeric(l, m, zss, nrho) for nrho in nrhos])
    ax.plot(nrhos, flux_numeric / flux_analytic, marker='o', label=f'l={l}, m={m}')

ax.legend()
ax.set_ylim(1)
ax.set_xlim(0)
ax.set_xlabel('$n_{r}$')
ax.set_ylabel('$\Phi_{pfsspy} / \Phi_{analytic}$')
plt.show()

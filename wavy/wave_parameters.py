#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from copy import deepcopy

def pseudo_wave_age(Hs, U10):
    # L. L. Fu and R. Glazman,
    # “The effect of the degree of wave development
    #  on the sea state bias in radar altimetry measurement”
    # J. Geophys. Res., vol. 96, no. C1, pp. 829–834, Jan. 1991.
    """
    dimensionless
    """
    return 3.25*(Hs*9.81/U10**2)**2

def altimeter_Tz(Hs, U10, GA=2, wa=None):
    if GA == 1:
        Tz = GA1(Hs, U10)
    elif GA == 2:
        Tz = GA2(Hs, U10, wa)
    return Tz

def GA2(Hs, U10, wa):
    # GA-2 algorithm from Remya G., Kumar, R., Basu, S. & Sarkar, A.
    # Altimeter-derived ocean wave period using genetic algorithm.
    # IEEE Geoscience and Remote Sensing Letters, 8(2), 354–358, 2010.
    # Tz = (((wa-5.78) / (wa+(U10/(Hs*((U10/Hs)+Hs)))))+(Hs+5.7))

    # GA-2 algorithm (re-factored)
    Tz = (wa - 5.78) / (wa + U10 / (U10 + Hs**2)) + Hs + 5.70
    return Tz

def GA1(Hs, U10):
    # GA-1 algorithm (re-factored)
    Tz = Hs + (8.01 * Hs / U10) - (U10 / 4.38) + 5.09
    return Tz

def mean_wave_energy_density(Hs):
    rho_water = 1027  # kg/m3
    g = 9.81
    return 1/16*rho_water*g*Hs**2

def group_velocity_deep_water(Tz):
    g = 9.81
    return g * Tz / (4 * np.pi)

def wave_power(E, cg):
    return E*cg/1000  # /1000 for kW

#!/usr/bin/env python

from larf.drift import mobility
from larf import units

import numpy as np
import matplotlib.pyplot as plt

def test_mobility():
    Emag = np.arange(0, 1000, 1)*units.V/units.cm
    mu = mobility(Emag)
    plt.clf()
    e_units = units.V/units.cm
    mu_units = (units.cm**2 / units.second) / units.V
    plt.loglog(Emag/e_units, mu/mu_units)
    plt.xlabel("Efield (V/cm)")
    plt.ylabel("Mobility (cm^2/s/V)")


def test_speed():
    Emag = np.arange(0, 1000, 1)*units.V/units.cm
    mu = mobility(Emag)
    s = mu*Emag

    e_units = units.V/units.cm
    s_units = units.mm/units.us

    plt.clf()
    plt.plot(Emag / e_units, s / s_units)
    
    plt.xlabel("Efield (V/cm)")
    plt.ylabel("Drift speed (mm/us)")

def test_speed_sou():
    Emag = np.arange(0, 1000, 1)*units.V/units.cm
    mu = mobility(Emag)
    s = mu*Emag

    plt.clf()
    plt.plot(Emag, s)
    
    plt.xlabel("Efield ([V/distance])")
    plt.ylabel("Drift speed ([distance/time])")

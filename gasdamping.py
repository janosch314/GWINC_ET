#
import os
import numpy as np
import const

def S_F_weiss(ifo):
    """
    Returns the freq. independent force PSD on a single test mass due to impinging residual gas particles as given by T0900509
    r = mirror radius [m]
    p = pressure of gas [Pa]
    m = molecular mass of particle [kg] (default = 2.99e-26 for water)
    T= temperaure [K] (default = 300K)
    """
    Lambda = ifo.Laser.Wavelength
    L = ifo.Infrastructure.Length
    kT = ifo.Infrastructure.Temp * const.kB
    P = ifo.Infrastructure.ResidualGas.pressure
    m = ifo.Infrastructure.ResidualGas.mass
    alpha = ifo.Infrastructure.ResidualGas.polarizability
    r=ifo.Materials.MassRadius


    S_F = 8 * P * np.sqrt(kT*m) * np.pi * r**2
    return S_F


def S_F_cavalleri(ifo):
    """
    Returns the freq. independent force PSD on a single test mass due to impinging residual gas particles as given by  Cavalleri et al 2009 (https://doi.org/10.1016/j.physleta.2010.06.041)
    Assumes radius to thickness ratio of mirror to be 1:1
    r = mirror radius [m]
    p = pressure of gas [Pa]
    m = molecular mass of particle [kg] (default = 2.99e-26 for water)
    T= temperaure [K] (default = 300K)
    """
    Lambda = ifo.Laser.Wavelength
    L = ifo.Infrastructure.Length
    kT = ifo.Infrastructure.Temp * const.kB
    P = ifo.Infrastructure.ResidualGas.pressure
    m = ifo.Infrastructure.ResidualGas.mass
    alpha = ifo.Infrastructure.ResidualGas.polarizability
    r=ifo.Materials.MassRadius
    
    S_F = P * np.sqrt((128/np.pi)*m*kT) * np.pi * r**2 * (1 + r/(2*r) + np.pi/4)
    return S_F

def calc_x_noise(f, S_F, ifo):
    """
    Returns the PSD for four test masses for a given force PSD.
    S_F = PSD of force noise
    M = mass of test mass [kg]
    f = frequency [Hz]
    """
    M=ifo.Suspension.Stage[0].Mass
    x = 4* S_F / (M * (2*np.pi*f)**2)**2
    return x

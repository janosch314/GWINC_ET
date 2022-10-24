import os
from gwinc import nb, const
from gwinc.ifo.noises import *
from gwinc.ifo import PLOT_STYLE
from susthnew import STNpy
from gasdamping import S_F_cavalleri
from gasdamping import calc_x_noise

from envnoise import (
        atmospheric_noise,
        cavern_noise,
        body_wave,
        rayleigh_wave,
        seismic_noise
        )
#newtonian_mitigation_factor = 3

class QuantumVacuum(nb.Budget):
    """Quantum Vacuum

    """
    style = dict(
        label='Quantum Vacuum',
        color='#ad03de',
    )

    noises = [
        QuantumVacuumAS,
        QuantumVacuumArm,
        QuantumVacuumSEC,
        QuantumVacuumFilterCavity,
        QuantumVacuumInjection,
        QuantumVacuumReadout,
        QuantumVacuumQuadraturePhase,
    ]
class SusThermal(nb.Noise):
    style = dict(
        label = 'Suspension Thermal',
        color='#0d75f8',
        )
    def calc(self):
        #STNpy return PSD
        #_, noise = STNRmodal(self.freq,self.ifo.Suspension,self.ifo)
        #violin = STNViol(self.freq, self.ifo.Suspension, self.ifo)
        [noise,noise_h,noise_v]=STNpy(self.freq, self.ifo)
        #turn into displacement PSD
        return (noise).real
            
class ExcessGasDampingH2(nb.Noise):
    """Residual gas damping for H2

    """
    style = dict(
        label='H$_2$ damping',
        color='xkcd:red orange',
        linestyle='--',
    )
    def calc(self):
        species = self.ifo.Infrastructure.ResidualGas.H2
        dam=calc_x_noise(self.freq,S_F_cavalleri(self.ifo,species),self.ifo)
        return dam

class ExcessGasDampingN2(nb.Noise):
    """Excess gas damping for N2

    """
    style = dict(
        label='N$_2$ damping',
        color='xkcd:emerald',
        linestyle='--',
    )
    def calc(self):
        species = self.ifo.Infrastructure.ResidualGas.N2
        dam=calc_x_noise(self.freq,S_F_cavalleri(self.ifo,species),self.ifo)
        return dam


class ExcessGasDampingH2O(nb.Noise):
    """Excess gas damping for H2O

    """
    style = dict(
        label='H$_2$O damping',
        color='xkcd:water blue',
        linestyle='--',
    )
    def calc(self):
        species = self.ifo.Infrastructure.ResidualGas.H2O
        dam=calc_x_noise(self.freq,S_F_cavalleri(self.ifo,species),self.ifo)
        return dam


class ExcessGasDampingO2(nb.Noise):
    """Excess gas damping for O2

    """
    style = dict(
        label='O$_2$ damping',
        color='xkcd:grey',
        linestyle='--',
    )
    def calc(self):
        species = self.ifo.Infrastructure.ResidualGas.O2
        dam=calc_x_noise(self.freq,S_F_cavalleri(self.ifo,species),self.ifo)
        return dam

            
class ExcessGas(nb.Budget):
    """Excess Gas

    """
    style = dict(
        label='Residual Gas',
        color='#add00d',
        linestyle='-',
    )

    noises = [
        ExcessGasScatteringH2,
        ExcessGasScatteringN2,
        ExcessGasScatteringH2O,
        ExcessGasScatteringO2,
        ExcessGasDampingH2,
        ExcessGasDampingN2,
        ExcessGasDampingH2O,
        ExcessGasDampingO2,
    ]

class Seismic(nb.Noise):
    style = dict(
        label = 'Seismic',
        color='#855700'
        )
    def calc(self):
        noise = seismic_noise(self.freq,self.ifo.Seismic)**2
        return noise
        
class NewtonianBodyWave(nb.Noise):
    style = dict(
        label = 'Body Wave',
        )
    def calc(self):
        noise = body_wave(self.freq,self.ifo.Seismic)**2
        return noise /  self.ifo.Seismic.Omicron**2

class NewtonianRayleighWave(nb.Noise):
    style = dict(
        label = 'Rayleigh Wave',
        )
    def calc(self):
        noise = rayleigh_wave(self.freq,self.ifo.Seismic)**2
        return noise/ self.ifo.Seismic.Omicron**2

class NewtonianCavern(nb.Noise):
    style = dict(
        label = 'Cavern',
        )
    def calc(self):
        noise = cavern_noise(self.freq,self.ifo.Seismic)**2
        return noise / self.ifo.Seismic.Omicron**2

class NewtonianAtmospheric(nb.Noise):
    style = dict(
        label = 'Atmospheric',
        )
    def calc(self):
        noise = atmospheric_noise(self.freq,self.ifo.Seismic)**2
        return noise/ self.ifo.Seismic.Omicron**2

class Newtonian(nb.Budget):
    """Newtonian"""
    style = dict(
        label = 'Newtonian Gravity',
        color='#15b01a'
        )
    noises = [
            NewtonianBodyWave,
            NewtonianRayleighWave,
            NewtonianCavern,
            NewtonianAtmospheric
            ]


class ETHF(nb.Budget):

    name = 'ETHF'

    noises = [
        QuantumVacuum,
        Seismic,
        Newtonian,
        SusThermal,
        CoatingBrownian,
        CoatingThermoOptic,
        SubstrateBrownian,
        SubstrateThermoElastic,
        ExcessGas,
        ITMThermoRefractive
    ]

    calibrations = [
        Strain,
    ]

    plot_style = PLOT_STYLE

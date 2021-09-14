from gwinc.ifo.noises import *
from gwinc.ifo import PLOT_STYLE
from gwinc.ifo.noises import arm_cavity
from susth import STNpy
from thermoelastic import substratethermoelastic
from envnoise import (
        atmospheric_noise,
        cavern_noise,
        body_wave,
        rayleigh_wave,
        seismic_noise
        )

newtonian_mitigation_factor = 3

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
        noise, _, _ = STNpy(self.freq, self.ifo)
        return noise

class SubThermalElastic(nb.Noise):
    style = dict(
        label = 'Substrate Thermal Elastic',
        color='#ffd700',
        linestyle='--',
        )
    def calc(self):
        cavity = arm_cavity(self.ifo)
        nITM = substratethermoelastic(
            self.freq, self.ifo.Materials, cavity.wBeam_ITM)
        nETM = substratethermoelastic(
            self.freq, self.ifo.Materials, cavity.wBeam_ETM)
        return (nITM + nETM) * 2
        
class Seismic(nb.Noise):
    style = dict(
        label = 'Seismic',
        color='#855700'
        )
    def calc(self):
        noise = seismic_noise(self.freq)**2
        return noise

class NewtonianBodyWave(nb.Noise):
    style = dict(
        label = 'Body Wave',
        )
    def calc(self):
        noise = body_wave(self.freq)**2
        return noise / newtonian_mitigation_factor**2

class NewtonianRayleighWave(nb.Noise):
    style = dict(
        label = 'Rayleigh Wave',
        )
    def calc(self):
        noise = rayleigh_wave(self.freq)**2
        return noise / newtonian_mitigation_factor**2

class NewtonianCavern(nb.Noise):
    style = dict(
        label = 'Cavern',
        )
    def calc(self):
        noise = cavern_noise(self.freq)**2
        return noise / newtonian_mitigation_factor**2

class NewtonianAtmospheric(nb.Noise):
    style = dict(
        label = 'Atmospheric',
        )
    def calc(self):
        noise = atmospheric_noise(self.freq)**2
        return noise / newtonian_mitigation_factor**2

class NewtonianNoise(nb.Budget):
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

class ETLF(nb.Budget):

    name = 'ETLF'

    noises = [
        QuantumVacuum,
        Seismic,
        NewtonianNoise,
        SusThermal,
        CoatingBrownian,
        CoatingThermoOptic,
        SubstrateBrownian,
        SubThermalElastic,
        ExcessGas,
    ]

    calibrations = [
        Strain,
    ]

    plot_style = PLOT_STYLE


    freq = '1:3000:4000'


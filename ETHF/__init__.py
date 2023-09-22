import os
from gwinc import nb, const
from gwinc import noise
from gwinc.ifo.noises import *
from gwinc.ifo import PLOT_STYLE
from susthnew import STNpy
from gasdamping import S_F_cavalleri
from gasdamping import calc_x_noise
from gwinc.ifo.noises import arm_cavity
        
from envnoise import (
        atmospheric_noise,
        cavern_noise,
        body_wave,
        rayleigh_wave,
        seismic_noise,
        HFseismic_noise
        )
#newtonian_mitigation_factor = 3

from gwinc.noise.quantum import (
    Quantum,
    QuantumRelShotNoise,
    QuantumRelGamma,
    QuantumXi,
)


    
class Coating(nb.Budget):
    """Coating Thermal

    """

    name = 'Coating'

    style = dict(
        label='Coating Thermal',
        color='#fe0002',
    )

    noises = [
        noise.coatingthermal.CoatingBrownian,
        noise.coatingthermal.CoatingThermoOptic,
    ]


class Substrate(nb.Budget):
    """Substrate Thermal

    """

    name = 'Substrate'

    style = dict(
        label='Substrate Thermal',
        color='#fb7d07',
    )

    noises = [
        noise.substratethermal.SubstrateBrownian,
        noise.substratethermal.SubstrateThermoElastic,
        noise.substratethermal.ITMThermoRefractive
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
            
class Seismic(nb.Noise):
    style = dict(
        label = 'Seismic',
        color='#855700'
        )
    def calc(self):
        noise = HFseismic_noise(self.freq,self.ifo)
        return noise
  
#class Seismic(nb.Noise):
  #  style = dict(
  #      label = 'Seismic',
   #     color='#855700'
   #     )
  #  def calc(self):
  #      noise = seismic_noise(self.freq,self.ifo.Seismic)**2
  #      return noise
        
class NewtonianBodyWave(nb.Noise):
    style = dict(
        label = 'Body Wave',
        color='#85a3b2'
        )
    def calc(self):
        noise = body_wave(self.freq,self.ifo.Seismic)**2
        return noise /  self.ifo.Seismic.Omicron**2

class NewtonianRayleighWave(nb.Noise):
    style = dict(
        label = 'Rayleigh Wave',
        color='#C20078'
        )
    def calc(self):
        noise = rayleigh_wave(self.freq,self.ifo.Seismic)**2
        return noise/ self.ifo.Seismic.Omicron**2

class NewtonianCavern(nb.Noise):
    style = dict(
        label = 'Cavern',
        color='#650021'
        )
    def calc(self):
        noise = cavern_noise(self.freq,self.ifo.Seismic)**2
        return noise / self.ifo.Seismic.Omicron**2

class NewtonianAtmospheric(nb.Noise):
    style = dict(
        label = 'Atmospheric',
        color='#ffa62b'
        )
    def calc(self):
        noise = atmospheric_noise(self.freq,self.ifo.Seismic)**2
        return noise/ self.ifo.Seismic.Omicron**2

class Newtonian(nb.Budget):
    """Newtonian"""
    style = dict(
        label = 'Newtonian Gravity',
        color='#15b01a',
        alpha=1
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
        Quantum,
        Seismic,
        Newtonian,
        SusThermal,
        Coating,
        Substrate,
        noise.residualgas.ResidualGas
    ]

    calibrations = [
        Strain,
    ]

    plot_style = PLOT_STYLE

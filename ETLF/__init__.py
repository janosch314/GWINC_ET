import os
import copy
from gwinc.ifo.noises import *
from gwinc.ifo import PLOT_STYLE
from gwinc import noise
from gwinc.ifo.noises import arm_cavity
from gwinc.ifo.noises import ifo_power
from Scattering import BackResult
from susthnew import STNpy
from gasdamping import S_F_cavalleri
from gasdamping import calc_x_noise
import coatingth
from gwinc.ifo.noises import arm_cavity

from thermoelastic import substratethermoelastic
from envnoise import (
        atmospheric_noise,
        cavern_noise,
        body_wave,
        rayleigh_wave,
        seismic_noise
        )

DEFAULT_FREQ = np.logspace(np.log10(1),np.log10(1000),3000)

def mirror_struct(ifo, tm):
    """Create a "mirror" Struct for a LIGO core optic

    This is a copy of the ifo.Materials Struct, containing Substrate
    and Coating sub-Structs, as well as some basic geometrical
    properties of the optic.

    """
    # NOTE: we deepcopy this struct since we'll be modifying it (and
    # it would otherwise be stored by reference)
    mirror = copy.deepcopy(ifo.Materials)
    optic = ifo.Optics.get(tm)
    coatingth.build_stacks(mirror, optic)
    mirror.update(optic)
    mirror.MassVolume = pi * mirror.MassRadius**2 * mirror.MassThickness
    mirror.MirrorMass = mirror.MassVolume * mirror.Substrate.MassDensity
    return mirror
    

from gwinc.noise.quantum import (
    Quantum,
    QuantumRelShotNoise,
    QuantumRelGamma,
    QuantumXi,
)


class CoatingBrownian(nb.Noise):
    """Coating Brownian

    """
    style = dict(
        label='Coating Brownian',
        color='#fe0002',
    )

    def calc(self):
        ITM = mirror_struct(self.ifo, 'ITM')
        ETM = mirror_struct(self.ifo, 'ETM')
        cavity = arm_cavity(self.ifo)
        wavelength = self.ifo.Laser.Wavelength
        nITM = coatingth.coating_brownian(
            self.freq, ITM, wavelength, cavity.wBeam_ITM
        )
        nETM = coatingth.coating_brownian(
            self.freq, ETM, wavelength, cavity.wBeam_ETM
        )
        return (nITM + nETM) * 2


class CoatingThermoOptic(nb.Noise):
    """Coating Thermo-Optic

    """
    style = dict(
        label='Coating Thermo-Optic',
        color='#02ccfe',
        linestyle='--',
    )

    def calc(self):
        wavelength = self.ifo.Laser.Wavelength
        materials = self.ifo.Materials
        ITM = mirror_struct(self.ifo, 'ITM')
        ETM = mirror_struct(self.ifo, 'ETM')
        cavity = arm_cavity(self.ifo)
        nITM, junk1, junk2, junk3 = coatingth.coating_thermooptic(
            self.freq, ITM, wavelength, cavity.wBeam_ITM,
        )
        nETM, junk1, junk2, junk3 = coatingth.coating_thermooptic(
            self.freq, ETM, wavelength, cavity.wBeam_ETM,
        )
        return (nITM + nETM) * 2
        
class Coating(nb.Budget):
    """Coating Thermal

    """

    name = 'Coating'

    style = dict(
        label='Coating Thermal',
        color='#fe0002',
    )

    noises = [
        CoatingBrownian,
        CoatingThermoOptic,
    ]

class SubstrateThermoElastic(nb.Noise):
    style = dict(
        label = 'Substrate Thermal Elastic',
        color='#f5bf03',
        linestyle='--',
        )
    def calc(self):
        cavity = arm_cavity(self.ifo)
        nITM = substratethermoelastic(
            self.freq, self.ifo.Materials, cavity.wBeam_ITM)
        nETM = substratethermoelastic(
            self.freq, self.ifo.Materials, cavity.wBeam_ETM)
        return (nITM + nETM) * 2
        
class ITMThermoRefractive(nb.Noise):

    style = dict(
        label='ITM Thermo-Refractive',
        color='#448ee4',
        linestyle='--',
    )

    def calc(self):
        power = ifo_power(self.ifo)
        gPhase = power.finesse * 2/np.pi
        cavity = arm_cavity(self.ifo)
        n = noise.substratethermal.substrate_thermorefractive(
            self.freq, self.ifo.Materials, cavity.wBeam_ITM, exact=True)
        return n * 2 / gPhase**2

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
        SubstrateThermoElastic,
        ITMThermoRefractive
    ]


class SuspensionThermal(nb.Noise):
    style = dict(
        label = 'Suspension Thermal',
        color='#0d75f8',
        )
    def calc(self):
        #_,noise = STNRmodal(self.freq, self.ifo.Suspension, self.ifo)
        #violin = STNViol(self.freq, self.ifo.Suspension, self.ifo)
        [noise,noise_h,noise_v]=STNpy(self.freq, self.ifo)
        return (noise).real
        
class SeismicHR(nb.Noise):
    style = dict(
        label = 'Rayleigh wave Horizental',
        color='#800000',
        linestyle='--'
        )
    def calc(self):
        noiseHR,_,_,_,_ = seismic_noise(self.freq,self.ifo.Seismic)
        return noiseHR
        
class SeismicHB(nb.Noise):
    style = dict(
        label = 'Body wave Horizental',
        color='#000080',
        linestyle='--'
        )
    def calc(self):
        _,noiseHB,_,_,_ = seismic_noise(self.freq,self.ifo.Seismic)
        return noiseHB
        
class SeismicVR(nb.Noise):
    style = dict(
        label = 'Rayleigh wave Vertical',
        color='#808000',
        linestyle='--'
        )
    def calc(self):
        _,_,noiseVR,_,_ = seismic_noise(self.freq,self.ifo.Seismic)
        return noiseVR
        
class SeismicVB(nb.Noise):
    style = dict(
        label = 'Body wave Vertical',
        color='#FFA500',
        linestyle='--'
        )
    def calc(self):
        _,_,_,noiseVB,_ = seismic_noise(self.freq,self.ifo.Seismic)
        return noiseVB

class SeismicTR(nb.Noise):
    style = dict(
        label = 'Rayleigh wave Tilt',
        color='#FF4500',
        linestyle='--'
        )
    def calc(self):
        _,_,_,_,noiseTR = seismic_noise(self.freq,self.ifo.Seismic)
        return noiseTR
        
    


class Seismic(nb.Budget):
    """Seismic"""
    style = dict(
        label = 'Seismic',
        color='#855700'
        )
    noises = [
            SeismicHR,
            SeismicHB,
            SeismicVR,
            SeismicVB,
            SeismicTR
            ]

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
        return noise/  self.ifo.Seismic.Omicron**2

class NewtonianCavern(nb.Noise):
    style = dict(
        label = 'Cavern Acoustic',
        color='#650021'
        )
    def calc(self):
        noise = cavern_noise(self.freq,self.ifo.Seismic)**2
        return noise / (self.ifo.Seismic.Omicron)**2

class NewtonianAtmospheric(nb.Noise):
    style = dict(
        label = 'Atmospheric',
        color='#ffa62b'
        )
    def calc(self):
        noise = atmospheric_noise(self.freq,self.ifo.Seismic)**2
        return noise/  (self.ifo.Seismic.Omicron)**2

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


class StrayLight(nb.Noise):
    """StrayLight"""
    name = 'StrayLight'
    style = dict(
        label='Back Scattering',
        color='#00EEEE',
    )
    def calc(self):
        noise = BackResult(self.freq,self.ifo)**2
        return noise
        
class ETLF(nb.Budget):

    name = 'ETLF'
    noises = [
        Quantum,
        Seismic,
        Newtonian,
        StrayLight,
        SuspensionThermal,
        Coating,
        Substrate,
        noise.residualgas.ResidualGas
    ]

    calibrations = [
        Strain,
    ]

    plot_style = PLOT_STYLE


    #freq = '1:3000:4000'


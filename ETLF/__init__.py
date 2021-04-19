from gwinc.ifo.noises import *
from gwinc.ifo import PLOT_STYLE
from susth import STNpy

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
        label = 'suspension thermal',
        color = 'orange'
        )
    def calc(self):
        noise, noise_h, noise_cv = STNpy(self.freq, self.ifo)
        return noise*(self.ifo.Infrastructure.Length)**2

class ETLF(nb.Budget):

    name = 'ETLF'

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
    ]

    calibrations = [
        Strain,
    ]

    plot_style = PLOT_STYLE
<<<<<<< HEAD
=======
    freq = '1:3000:4000'
>>>>>>> master

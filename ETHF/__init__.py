import os
from gwinc import nb, const
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
        label = 'Suspension Thermal',
        color = 'orange'
        )
    def calc(self):
        #STNpy return PSD
        noise, _, _ = STNpy(self.freq, self.ifo)
        #turn into displacement PSD
        return noise

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
    ]

    calibrations = [
        Strain,
    ]

    plot_style = PLOT_STYLE

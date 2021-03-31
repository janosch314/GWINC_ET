from gwinc.ifo.noises import *
from gwinc import load_budget
import numpy as np

DEFAULT_FREQ = '1:3000:6000'

class BudgetWrapper(nb.Noise):
    ifo_name = ''

    def load(self):
        self.budget = load_budget(self.ifo_name, freq=self.freq).run()
    
    def calc(self):
        return self.budget.psd

class ETLF(BudgetWrapper):
    """ET Low-Frequency"""
    ifo_name = 'ETLF' 

class ETHF(BudgetWrapper):
    """ET High-Frequency"""
    ifo_name = 'ETHF' 


def invsum(traces):
    return 1.0/np.sum(traces, axis=0)

class BudgetEnvelope(nb.Budget):
    def calc_trace(self, calibration=1, calc=True, _precomp=None):
        """Calculate all budget noises and return BudgetTrace object

        The envelope of all budgets will be given as the total noise,
        instead of the sum of all noises.

        `calibration` is always set to 1 for this calculation.
        TODO: is that ok?

        If `calc` is False, the noise will not be calculated and the
        trace PSD will be None.  This is useful for just getting the
        trace style info.

        """
        if _precomp is None:
            _precomp = dict()

        _cals = {}
        if calc:
            for name, cal in self._cal_objs.items():
                _cals[name] = cal._calc(_precomp)
        budget = []
        for name in self._noise_objs:
            trace = self.calc_noise(
                name,
                calibration=1.0, #TODO???
                calc=calc,
                _cals=_cals,
                _precomp=_precomp,
            )
            budget.append(trace)
        
        total = invsum([1.0/trace.psd for trace in budget if trace.name in self._budget_noises])
        return self._make_trace(
            psd=total, budget=budget
        )

class ET(BudgetEnvelope):
    name = 'ET-D'
    freq = DEFAULT_FREQ

    noises = [
        ETLF,
        ETHF
    ]
    
    calibrations = [] # do not do the to-strain conversion
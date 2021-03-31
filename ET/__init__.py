from gwinc.ifo.noises import *
from gwinc import load_budget
import numpy as np
from numpy.core.shape_base import _accumulate

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


def invsum(data):
    return 1.0/np.nansum([1.0/x for x in data], axis=0)

class ET(nb.Budget):
    name = 'Einstein Telescope'
    freq = DEFAULT_FREQ

    noises = [
        ETLF,
        ETHF
    ]
    
    calibrations = [] # do not do the to-strain conversion

    accumulate = invsum # calculate envelope of noise curves instead of sum

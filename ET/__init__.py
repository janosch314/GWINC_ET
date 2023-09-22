import os
import numpy as np
from gwinc.ifo.noises import *
from gwinc import load_budget
import scipy.interpolate

DEFAULT_FREQ = np.logspace(np.log10(1),np.log10(10000),3000)

def get_local_filename(*parts):
    return os.path.join(os.path.dirname(__file__), *parts)

class BudgetWrapper(nb.Noise):
    ifo_name = ''

    def load(self):
        self.budget = load_budget(self.ifo_name, freq=DEFAULT_FREQ).run()
    
    def calc(self):
        return self.budget.psd

class ETLF(BudgetWrapper):
    """ET Low-Frequency"""
    ifo_name = 'ETLF'
    style = dict(
        color='blue',
        linewidth=3
    )

class ETHF(BudgetWrapper):
    """ET High-Frequency"""
    ifo_name = 'ETHF'
    style = dict(
        color='red',
        linewidth=3
    )

class ETDesignReport(nb.Noise):
    style = dict(
        label='ET-D',
        color='black',
        linestyle='--',
        linewidth=2
    )

    def load(self):
        data = get_local_filename('et_d.txt')
        freq, asd = np.loadtxt(data).T
        func = scipy.interpolate.interp1d(freq, asd**2)
        self.psd = func(DEFAULT_FREQ)
        
    def calc(self):
        return self.psd

def invsum(data):
    return 1.0/np.nansum([1.0/x for x in data], axis=0)

class ET(nb.Budget):
    name = 'Einstein Telescope'
    freq = DEFAULT_FREQ
    style = dict(
        color='black',
        alpha=1,
        linewidth=4
    )

    noises = [
        ETLF,
        ETHF
    ]

    # no to-strain conversion, as budgets are already in strain
    calibrations = [] 

    references = [
        ETDesignReport
    ]
    
    accumulate = invsum # calculate envelope of noise curves instead of sum

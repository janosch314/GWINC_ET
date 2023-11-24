To use this package for plotting noise budgets for Einstein Telescope one needs a working installation of Python, a PyGWINC and an inspiral-range packages installed.

More information and a tutorial can be found on the [WP wiki](https://wiki.et-gw.eu/ISB/Interferometer/ObservatoryDesignAndNoiseBudget/WebHome)

**To install PyGWINC**:

_[Option 1:]_ (Recomended) Install with git from the PyGwinc home page: https://git.ligo.org/gwinc/pygwinc, Run "git check out superQK" 

_[Option 2:]_ Run the folloing command:
`pip install gwinc`  (Not adapted to current ET files)

**To install inspiral-range**:

_[Option 1:]_ Run the folloing command: 
`pip install inspiral-range`

_[Option 2:]_ Download from the inspiral-range home page: https://git.ligo.org/gwinc/inspiral-range

**HOW TO USE**

- Download this repository to any folder on your machine
- Install PyGWINC and inspiral_range packages
- Run one of the Jupyter notebook files as described below OR create your own from templates


**To plot the noise budgets** for ET-HF and ET-LF, run the [Run.ipynb](https://gitlab.et-gw.eu/et/isb/interferometer/wpiii.1-observatory-design-and-noise-budget./-/blob/master/Run.ipynb) Jupyter notebook file

**To plot the Redshift vs. Source Mass plot**, run [Run_ horizon.ipynb](https://gitlab.et-gw.eu/et/isb/interferometer/wpiii.1-observatory-design-and-noise-budget./-/blob/master/Run_%20horizon.ipynb) Jupyter notebook file

For an **example of how to tweak individual parameters** and compare the new and the old noise budgets, run [Param_change_example.ipynb](https://gitlab.et-gw.eu/et/isb/interferometer/wpiii.1-observatory-design-and-noise-budget./-/blob/master/Param_change_example.ipynb)

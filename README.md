To use this package for plotting noise budgets for Einstein Telescope one needs a working installation of Python, a PyGWINC and an inspiral-range packages installed.

**To install PyGWINC**:

1. For the full process, follow instructions on the [wiki page](https://wiki.et-gw.eu/ISB/Interferometer/ObservatoryDesignAndNoiseBudget/BeginnerSGuideToInstallingAndRunningPygwinc).

2. If you know what you are doing, run `pip install git+https://git.ligo.org/gwinc/pygwinc.git@superQK` to install the pygwinc with from the correct branch required to compute ET noise budget

**To pull ET model into main colab folder**:

git init
git remote add origin https://github.com/janosch314/GWINC_ET.git
git pull origin master

**To install inspiral-range**:

_[Option 1:]_ Run the folloing command: 
`pip install inspiral-range`

_[Option 2:]_ Download from the inspiral-range home page: https://git.ligo.org/gwinc/inspiral-range

**HOW TO USE**

- Download this repository to any folder on your machine
- Install PyGWINC and inspiral_range packages
- Run one of the Jupyter notebook files as described below OR create your own from templates

**To plot the noise budgets** for ET-HF and ET-LF, run the [Run_example.ipynb](https://gitlab.et-gw.eu/et/isb/interferometer/wpiii.1-observatory-design-and-noise-budget./-/blob/master/Run_example.ipynb) Jupyter notebook file

**To plot the Redshift vs. Source Mass plot**, run [Run_ horizon.ipynb](https://gitlab.et-gw.eu/et/isb/interferometer/wpiii.1-observatory-design-and-noise-budget./-/blob/master/Run_horizon_demo.ipynb) Jupyter notebook file

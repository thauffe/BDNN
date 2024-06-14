# BDNN simulator

This is a snapshot of the [BDNNsim](https://github.com/thauffe/BDNNsim) library for python v3. It simulates lineage diversification under a birth-death process, the evolution of categorical and continuous traits, and a synthetic fossil record of occurrence data. It requires the python libraries numpy, scipy, pandas, and dendropy.

Basic usage to simulate constant diversification without lineage-specific rates in python is:

```
import sys
import os


# Set working directory to the download location of the GitHub repository
#########################################################################
wd = '/.../BDNN'
name = 'ExampleConstant'

# Import the BDNNsim library
########################################
directory_BDNNsim = os.path.join(wd, 'BDNNsim')
sys.path.insert(0, r'%s' % directory_BDNNsim)
from bdnn_simulator import *

# Specify parameters of birth-death simulation
##############################################
bd_sim = bdnn_simulator(rangeSP = [200, 300],                 # min/max size data set
                        root_r = [35., 35.],                  # range root ages (here fixed)
                        rangeL = [0.2, 0.2],                  # range of randomly drawn birth rate (here fixed)
                        rangeM = [0.1, 0.1],                  # range of randomly drawn death rate (here fixed)
                        n_cont_traits = [1, 1],               # number of continuous traits (here fixed to 1)
                        cont_traits_sigma_clado = [0.2, 0.2], # jump parameter for continuous trait at speciation event (here fixed)
                        cont_traits_sigma = [0.02, 0.02],     # evolutionary rates for continuous traits (here fixed)
                        n_cat_traits = [1, 1],                # range for the number of categorical traits (here fixed to 1)
                        n_cat_traits_states = [2, 2],         # range number of states for categorical trait (here fixed to two states)
                        cat_traits_diag = 0.9)                # Probability of no state change at speciation for categorical traits


# Run birth-death simulation
############################
res_bd = bd_sim.run_simulation(verbose = True)


# Sampling simulation
#####################
fossil_sim = fossil_simulator(range_q = [0.5, 1.5]) # Zero shifts in sampling rate over time
sim_fossil = fossil_sim.run_simulation(res_bd['ts_te'])


# Write input files for PyRate analysis
#######################################
output_wd = os.path.join(wd, 'Simulations') # Writes fossil occurrences in the BDNN/Simulations directory
write_PyRate = write_PyRate_files(output_wd = output_wd, name = name)
name_file = write_PyRate.run_writter(sim_fossil, res_bd)
```

The scripts for the nine diversification scenarios [Simulations](https://github.com/thauffe/BDNN/tree/main/Simulations) directory provide the exact setting used in the manuscript.

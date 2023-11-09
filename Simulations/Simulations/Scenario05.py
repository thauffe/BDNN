import copy
import sys
import os
import subprocess
import numpy as np

user = 'thauffe'

sys.path.insert(0, r'/home/%s/BDNN/BDNNsim' % user)
from bdnn_simulator import *

# Specify parameters of birth-death simulation
##############################################
rnd_seed = int(os.getenv('SLURM_ARRAY_TASK_ID', 0)) + 999
scenario = 'Scenario05'
name = str(rnd_seed)

bd_sim = bdnn_simulator(s_species = 1,  # number of starting species
                        rangeSP = [200, 300],  # min/max size data set
                        minEX_SP = 0,  # minimum number of extinct lineages allowed
                        minExtant_SP = 2, # minimum number of extant lineages
                        root_r = [35., 35.],  # range root ages
                        rangeL = [0.1, 0.1],  # range of birth rates
                        rangeM = [0.05, 0.05],  # range of death rates
                        scale = 100.,
                        p_mass_extinction = 0.0,
                        magnitude_mass_ext = [0.0, 0.0],
                        poiL = 0,  # expected number of birth rate shifts
                        poiM = 0,  # expected number of death rate shift
                        range_linL = [0.0, 0.0],
                        range_linM = [0.0, 0.0],
                        fixed_Mtt = np.array([[35., 0.01], [0.0, 0.1]]),
                        n_cont_traits = [1, 1],  # number of continuous traits
                        cont_traits_sigma_clado = [0.2, 0.2],
                        cont_traits_sigma = [0.02, 0.02], # evolutionary rates for continuous traits
                        cont_traits_cor = [-1, 1],  # evolutionary correlation between continuous traits
                        cont_traits_Theta1 = [0.0, 0.0], # morphological optima; 0 is no directional change from the ancestral values
                        cont_traits_alpha = [0.0, 0.0],
                        n_cat_traits = [1, 1],
                        n_cat_traits_states = [2, 2], # range number of states for categorical trait
                        cat_traits_ordinal = [False, False],
                        cat_traits_dir = 2,
                        cat_traits_diag = 0.9,
                        cat_traits_effect = np.array([[5., 5.],[5., 5.]]),
                        cat_traits_effect_decr_incr = np.array([[False, False],[False, False]]),
                        n_areas = [1, 1],
                        dispersal = [0.005, 0.01],
                        extirpation = [0.05, 0.2],
                        seed = rnd_seed)  # if > 0 fixes the random seed to make simulations reproducible


# Run birth-death simulation
############################
res_bd = bd_sim.run_simulation(verbose = False)


# Sampling simulation
#####################
fossil_sim = fossil_simulator(range_q = [0.5, 1.5],
                              range_alpha = [1000.0, 1000.0],
                              poi_shifts = 0,
                              seed = rnd_seed)
sim_fossil = fossil_sim.run_simulation(res_bd['ts_te'])


# Write input files for PyRate analysis
#######################################
write_PyRate = write_PyRate_files(output_wd = '/home/%s/BDNN/Simulations/%s' % (user, scenario),
                                  delta_time = 1.0,
                                  name = name)
name_file = write_PyRate.run_writter(sim_fossil, res_bd)


# Run PyRate
############
# RJMCMC
RJMCMC_run = subprocess.run(['python3', '/home/%s/BDNN/PyRate/PyRate.py' % user,
                             '/home/%s/BDNN/Simulations/%s/%s/%s.py' % (user, scenario, name_file, name_file),
                             '-A', '4',
                             '-mHPP',
                             '-n', '10000000',
                             '-s', '10000',
                             '-p', '100000'])
# BDNN
BDNN_run = subprocess.run(['python3', '/home/%s/BDNN/PyRate/PyRate.py' % user,
                           '/home/%s/BDNN/Simulations/%s/%s/%s.py' % (user, scenario, name_file, name_file),
                           '-trait_file', '/home/%s/BDNN/Simulations/%s/%s/%s_traits.csv' % (user, scenario, name_file, name_file),
                           '-A', '0',
                           '-mHPP',
                           '-BDNNmodel', '1',
                           '-n', '10000000',
                           '-s', '10000',
                           '-p', '100000'])


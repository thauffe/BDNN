import copy
import sys
import os
import subprocess
import numpy as np
import pandas as pd


# Set seed, simulation scenario, and replicate name
###################################################
rnd_seed = int(os.getenv('SLURM_ARRAY_TASK_ID', 0))
scenario = 'Scenario09'
name = str(rnd_seed)
threads = str(1)


# Setup paths
#############
wd = '/home/rcooper/BDNN/'
simulations_wd = '/data/users/rcooper/BDNN/Simulations'
dir_BDNNsim = os.path.join(wd, 'BDNNsim')
dir_scenario = os.path.join(simulations_wd, scenario)
path_PyRate = os.path.join(wd, 'PyRate', 'PyRate.py')


# Import BDNN simulator
#######################
sys.path.insert(0, r'%s' % dir_BDNNsim)
from bdnn_simulator import *


# Specify parameters of birth-death simulation
##############################################
bd_sim = bdnn_simulator(rangeSP=[200, 300],                 # min/max size data set
                        minEX_SP=2,                         # minimum number of extinct lineages allowed
                        minExtant_SP=2,                     # minimum number of extant lineages
                        root_r=[35., 35.],                  # range root ages
                        rangeL=[0.5, 0.5],                  # range of birth rates
                        rangeM=[0.4, 0.4],                  # range of death rates
                        n_cont_traits=[1, 1],               # number of continuous traits
                        cont_traits_sigma_clado=[0.2, 0.2], #
                        cont_traits_sigma=[0.02, 0.02],     # evolutionary rates for continuous traits
                        cont_traits_effect_sp=np.array([[[ [0.8, 0.8] ]]]),
                        cont_traits_effect_ex=np.array([[[ [0.8, 0.8] ]]]),
                        cont_traits_effect_bellu_sp=np.array([[[ [1, 1] ]]]),
                        cont_traits_effect_bellu_ex=np.array([[[ [1, 1] ]]]),
                        n_cat_traits=[1, 1],                # range number of categorical traits
                        n_cat_traits_states=[2, 2],         # range number of states for categorical trait
                        cat_traits_diag=0.9,                # transition probability between categorical states at speciation event
                        cat_traits_min_freq=[0.3],
                        seed=rnd_seed)                      # fixes the random seed to make simulations reproducible


# Run birth-death simulation
############################
res_bd = bd_sim.run_simulation(verbose=False)


# Sampling simulation
#####################
fossil_sim = fossil_simulator(range_q=[0.5, 5.0],
                              range_alpha=[0.5, 5.0],
                              fixed_shift_times=np.array([28.1, 23.03, 20.44, 15.97, 13.82, 11.63, 7.246, 5.33, 3.6, 2.58, 1.8, 0.781, 0.126]),
                              q_loguniform=True,
                              alpha_loguniform=True,
                              seed=rnd_seed)
sim_fossil = fossil_sim.run_simulation(res_bd['ts_te'])


# Write input files for PyRate analysis
#######################################
write_PyRate = write_PyRate_files(output_wd=dir_scenario, name=name)
_ = write_PyRate.run_writter(sim_fossil, res_bd, num_pvr=2, write_tree=True)


# Write occurrences and lineages-through-time to text file
##########################################################
write_occurrence_table(sim_fossil, dir_scenario, name_file = name)

write_ltt(res_bd, dir_scenario, name)


# Plot tree and traits
######################
subprocess.run(['Rscript',
                os.path.join(simulations_wd, 'Scripts', 'PlotPhylo.R'),
                os.path.join(dir_scenario, name)])


# Replace influential trait
###########################
trait_file_in = os.path.join(dir_scenario, name, name + '_traits_pvr.csv')
traits = pd.read_csv(trait_file_in, sep='\t')
traits['cont_trait_0'] = np.random.normal(size=traits.shape[0])
trait_file_out = os.path.join(dir_scenario, name, name + '_traits_pvr_randcont.csv')
traits.to_csv(trait_file_out, header = True, sep = '\t', index = False, na_rep = 'NA')


# Run PyRateBDNN
################
BDNN_run = subprocess.run(['python3', path_PyRate,
                           os.path.join(dir_scenario, name, '%s.py' % name),
                           '-out', '_pvr',
                           '-seed', name,
                           '-trait_file', os.path.join(dir_scenario, name, '%s_traits_pvr.csv' % name),
                           '-A', '0',
                           '-mG',
                           '-qShift', os.path.join(dir_scenario, name, '%s_q_epochs.txt' % name),
                           '-pP', '1.5', '0.0',
                           '-BDNNmodel', '1',
                           '-n', '10000000',
                           '-s', '10000',
                           '-p', '100000',
                           '-thread', threads, '0'])

# Postprocessing
################
burnin = 0.25

# RTT
subprocess.run(['python3', path_PyRate,
                '-plotBDNN', os.path.join(dir_scenario, name, 'pyrate_mcmc_logs', '%s_pvr_reg1_G_BDS_BDNN_16_8Tc' % name),
                '-b', str(burnin)])


# PDP
subprocess.run(['python3', path_PyRate,
                '-plotBDNN_effects', os.path.join(dir_scenario, name, 'pyrate_mcmc_logs', '%s_pvr_reg1_G_BDS_BDNN_16_8Tc_mcmc.log' % name),
                '-plotBDNN_transf_features', os.path.join(dir_scenario, name, '%s_backscale_cont_traits.csv' % name),
                '-b', str(burnin),
                '-thread', threads, '0'])


# Predictor importance
subprocess.run(['python3', path_PyRate,
                '-BDNN_pred_importance', os.path.join(dir_scenario, name, 'pyrate_mcmc_logs', '%s_pvr_reg1_G_BDS_BDNN_16_8Tc_mcmc.log' % name),
                '-plotBDNN_transf_features', os.path.join(dir_scenario, name, '%s_backscale_cont_traits.csv' % name),
                '-b', str(burnin),
                '-thread', threads, '0'])


# Run PyRateBDS
###############
RJMCMC_run = subprocess.run(['python3', path_PyRate,
                             os.path.join(dir_scenario, name, '%s.py' % name),
                             '-seed', name,
                             '-A', '4',
                             '-mG',
                             '-qShift', os.path.join(dir_scenario, name, '%s_q_epochs.txt' % name),
                             '-pP', '1.5', '0.0',
                             '-n', '10000000',
                             '-s', '10000',
                             '-p', '100000'])


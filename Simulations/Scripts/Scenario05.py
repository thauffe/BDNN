import copy
import sys
import os
import subprocess
import numpy as np


# Set seed, simulation scenario, and replicate name
###################################################
rnd_seed = int(os.getenv('SLURM_ARRAY_TASK_ID', 0))
scenario = 'Scenario05'
name = str(rnd_seed)
threads = str(1)


# Setup paths
#############
wd = '/home/thauffe/BDNN/'
simulations_wd = '/data/users/thauffe/BDNN/Simulations'
dir_BDNNsim = os.path.join(wd, 'BDNNsim')
dir_scenario = os.path.join(simulations_wd, scenario)
path_PyRate = os.path.join(wd, 'PyRate', 'PyRate.py')


# Import BDNN simulator
#######################
sys.path.insert(0, r'%s' % dir_BDNNsim)
from bdnn_simulator import *


# Specify parameters of birth-death simulation
##############################################
bd_sim = bdnn_simulator(rangeSP=[200, 300],                            # min/max size data set
                        minEX_SP=2,                                    # minimum number of extinct lineages allowed
                        minExtant_SP=2,                                # minimum number of extant lineages
                        root_r=[35., 35.],                             # range root ages
                        rangeL=[0.1, 0.1],                             # range of birth rates
                        rangeM=[0.05, 0.05],                           # range of death rates
                        fixed_Mtt=np.array([[35., 0.01], [0.0, 0.1]]), #
                        n_cont_traits=[1, 1],                          # number of continuous traits
                        cont_traits_sigma_clado=[0.2, 0.2],
                        cont_traits_sigma=[0.02, 0.02],                # evolutionary rates for continuous traits
                        n_cat_traits=[1, 1],
                        n_cat_traits_states=[2, 2],                    # range number of states for categorical trait
                        cat_traits_diag=0.9,
                        cat_traits_effect=np.array([[5., 5.],[5., 5.]]),
                        cat_traits_effect_decr_incr=np.array([[False, False], [False, False]]),
                        cat_traits_min_freq=[0.3],
                        seed=rnd_seed)


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




# Run PyRateBDNN
################
BDNN_run = subprocess.run(['python3', path_PyRate,
                           os.path.join(dir_scenario, name, '%s.py' % name),
                           '-out', '_pvr_reg1',
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


# PDP for influential feature
subprocess.run(['python3', path_PyRate,
                '-BDNN_interaction', os.path.join(dir_scenario, name, 'pyrate_mcmc_logs', '%s_pvr_reg1_G_BDS_BDNN_16_8Tc_mcmc.log' % name),
                '-BDNN_groups', '{\"cat_trait_0\": [\"cat_trait_0\"]}',
                '-b', str(burnin),
                '-thread', threads, '0'])
subprocess.run(['python3', path_PyRate,
                '-BDNN_interaction', os.path.join(dir_scenario, name, 'pyrate_mcmc_logs', '%s_pvr_reg1_G_BDS_BDNN_16_8Tc_mcmc.log' % name),
                '-BDNN_groups', '{\"cat_trait_0\": [\"cat_trait_0\"], \"time\": [\"time\"]}',
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


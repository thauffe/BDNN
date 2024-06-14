#!/usr/bin/env bash

#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=10-00:00:00
#SBATCH --job-name=BDNN_S10
#SBATCH --mail-user=torsten.hauffe@unifr.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/home/thauffe/output/output_%j.o
#SBATCH --error=/home/thauffe/error/error_%j.e
#SBATCH --array=1-100%100

user="thauffe"
burnin="0.25"
regprior="1"
scenario="10"


python3 /home/${user}/BDNN/ClusterScripts/Scenario${scenario}.py


python3 /home/${user}/BDNN/PyRate/PyRate.py \
/data/users/${user}/BDNN/Simulations/Scenario${scenario}/${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_ID}.py \
-out _pvr_reg${regprior} \
-seed ${SLURM_ARRAY_TASK_ID} \
-trait_file /data/users/${user}/BDNN/Simulations/Scenario${scenario}/${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_ID}_traits_pvr_randcat.csv \
-A 0 \
-mG \
-qShift /data/users/${user}/BDNN/Simulations/Scenario${scenario}/${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_ID}_q_epochs.txt \
-pP 1.5 0.0 \
-BDNNmodel 1 \
-BDNNreg ${regprior} \
-n 10000000 -s 10000 -p 100000


python3 /home/${user}/BDNN/PyRate/PyRate.py \
-plotBDNN /data/users/${user}/BDNN/Simulations/Scenario${scenario}/${SLURM_ARRAY_TASK_ID}/pyrate_mcmc_logs/${SLURM_ARRAY_TASK_ID}_pvr_reg${regprior}_G_BDS_BDNN_16_8Tc \
-b ${burnin}


python3 /home/${user}/BDNN/PyRate/PyRate.py \
-plotBDNN_effects /data/users/${user}/BDNN/Simulations/Scenario${scenario}/${SLURM_ARRAY_TASK_ID}/pyrate_mcmc_logs/${SLURM_ARRAY_TASK_ID}_pvr_reg${regprior}_G_BDS_BDNN_16_8Tc_mcmc.log \
-plotBDNN_transf_features /data/users/${user}/BDNN/Simulations/Scenario${scenario}/${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_ID}_backscale_cont_traits.csv \
-b ${burnin}


python3 /home/${user}/BDNN/PyRate/PyRate.py \
-BDNN_pred_importance /data/users/${user}/BDNN/Simulations/Scenario${scenario}/${SLURM_ARRAY_TASK_ID}/pyrate_mcmc_logs/${SLURM_ARRAY_TASK_ID}_pvr_reg${regprior}_G_BDS_BDNN_16_8Tc_mcmc.log \
-plotBDNN_transf_features /data/users/${user}/BDNN/Simulations/Scenario${scenario}/${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_ID}_backscale_cont_traits.csv \
-b ${burnin}

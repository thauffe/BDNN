#!/usr/bin/env bash

#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=10-00:00:00
#SBATCH --job-name=BDNN_S11
#SBATCH --mail-user=torsten.hauffe@unifr.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/home/thauffe/output/output_%j.o
#SBATCH --error=/home/thauffe/error/error_%j.e
#SBATCH --array=1-100%100

user="thauffe"
scenario="11"



python3 /home/${user}/BDNN/PyRate/PyRate.py \
/data/users/${user}/BDNN/Simulations/Scenario${scenario}/${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_ID}.py \
-seed ${SLURM_ARRAY_TASK_ID} \
-A 4 \
-mG \
-qShift /data/users/${user}/BDNN/Simulations/Scenario${scenario}/${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_ID}_q_epochs.txt \
-pP 1.5 0.0 \
-n 10000000 -s 10000 -p 100000




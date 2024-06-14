#!/usr/bin/env bash

#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=10-00:00:00
#SBATCH --job-name=BDNN_S01
#SBATCH --mail-user=torsten.hauffe@unifr.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/home/thauffe/output/output_%j.o
#SBATCH --error=/home/thauffe/error/error_%j.e
#SBATCH --array=1-100%100

python3 /home/thauffe/BDNN/ClusterScripts/Scenario01.py

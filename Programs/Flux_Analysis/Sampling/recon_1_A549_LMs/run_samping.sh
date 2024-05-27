#!/bin/bash
#SBATCH --job-name=example_job
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
module load miniconda
conda activate flux_env
python "/home/adg66/Flux_Code/Programs/Flux_Analysis/Sampling/recon_1_A549_LMs/all_sampling_sbatch_McCleary.py"
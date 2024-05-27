#!/bin/sh
#SBATCH --job-name=multithreadSSA		# Job name
#SBATCH --mail-type=NONE            		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=agoetz@ufl.edu 		# Where to send mail	
#SBATCH --nodes=1                     	# Use one node
#SBATCH --ntasks=1                    	# Run a single task
#SBATCH --cpus-per-task=1				# Number of CPU cores to use
#SBATCH --mem=5gb                   	# Memory limit
#SBATCH --time=200:00:00               	# Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.out   	# Standard output and error log
#SBATCH --array=0-319

pwd; hostname; date

python /blue/pdixit/agoetz/Flux_Code_New/Programs/Flux_Analysis/Sampling/recon_1b_t_cells_gene_Day_2_single_size/all_sampling_sbatch.py $SLURM_ARRAY_TASK_ID 


date

Run the files in the following order
1. gene_HR_warmup_all_mice.py
2. all_sampling_sbatch.py
Number 2 will likely take some work to run in an efficient manner. Ideally a job batch script should be run on a cluster.
The job should pass an array of 0-319. An example can be found in the "example.sh" file contained in the same directory as this one.

The gurobi key will likely cause the greatest issue, for the hipergator where this code was run, the line:
recon_flux_model.add_gp_key_env_to_model('grb-ts.ufhpc')
is used
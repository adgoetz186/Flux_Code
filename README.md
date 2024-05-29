# A Bayesian framework for integrative analysis of cellular metabolism
This code contains a flux balance model class which allows for fast and easy handling of flux balance models along with many tools for probing the high dimensional flux polytope. It is used to obtain results for the publication titled, "A global view of T cell metabolism in systemic lupus erythematosus". DOI: 10.3389/fimmu.2024.1371708

# Reproducing workflow of paper
While this code is used in publications it is still being worked on and updated. As a result unique instructions for reproducing published results will be placed in relevant folders. This section details the order to follow to reproduce results for a given paper.

## A global view of T cell metabolism in systemic lupus erythematosus
Data file can be downloaded from https://doi.org/10.7910/DVN/CQBIXG. The "Data.7z" file should be unzipped, the "_SLE" removed from its name, and saved to the "Flux_Code" folder. The processed sample file is located in Data/HR/HR_fin/mat_file. The actual sampling generates far more points, many of which are dropped to insure less correlated samples. These points are also availible being contained in the "B6_mice.7z" and "TC_mice.7z" files of the repository. 

While I try to not break any past versions with updates, the commit "Final version for a global view of T cell metabolism in systemic lupus erythematosus" contains a version which should be used to be safe.

### Loading flux balance model for sampling

Go navagate to Flux_Code/Programs/Flux_Analysis/Model_Loading/Loading_recon_1_t_cells_with_lupus_gene_support_12_11_2023
Run the files in the following order:
1. from_mat.py
2. curate_model.py (This compiles gene data after rxns are added)
3. flux_array_to_experimental_align_file.py
3. experimental_align_t_cell_initial_fit.py
4. experimental_align_t_cell_single_size.py
5. essential_flux_dataframe.py
6. consensus_model_from_dataframes.py
7. post_min_model_alterations_single_size.py
8. generate_internal_S_mat.py
The output of 8 is a matrix "S" which should be processed in matlab with the command "NS = null(S,"rational")"
The output of this should be saved as NS in the same folder.
This can be all done in python but any method which does not generate a sparse nullspace will lead to significant increases in testing for thermo feasibility
(numpy does not natively seem to be able to do this and it is prohibitively costly in sympy)

### Sampling

Go navagate to Flux_Code/Programs/Flux_Analysis/Sampling/recon_1_t_cells_gene_integration_12_11_23
Run the files in the following order
1. gene_HR_warmup_all_mice.py
2. all_sampling_sbatch.py
Number 2 will likely take some work to run in an efficient manner. Ideally a job batch script should be run on a cluster.
The job should pass an array of 0-319. An example can be found in the "example.sh" file contained in the same directory as this one.

The gurobi key will likely cause the greatest issue, for the hipergator where this code was run, the line:
recon_flux_model.add_gp_key_env_to_model('grb-ts.ufhpc')
is used

# Possible issues
Errors with setting working directory
There can be issues in getting the right path, this seems to often fix the issue
1. go to the python "site-packages" folder
2. create a new text file
3. add the absolute path to "Flux_Code"
4. save with the .pth extension
Alternatively there can be issues if your cwd starts at Flux_Code, adding 'if not Path(os.getcwd()).parts[-1] == "Flux_Code":' right above 'path_to_FC = ""' can fix this

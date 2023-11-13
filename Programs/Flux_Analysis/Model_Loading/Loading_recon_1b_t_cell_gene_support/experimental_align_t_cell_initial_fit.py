import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import scipy.io as sio
from pathlib import Path
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model



# _____ Setting the CWD to be Flux_Code BEGIN _____
# Flux_Code path here
path_to_FC = ""
if path_to_FC == "":
	try:
		# Obtains the location of the Cell_signaling_information folder if it is in the cwd parents
		path_to_FC = Path.cwd().parents[[Path.cwd().parents[i].parts[-1] for i in range(len(Path.cwd().parts) - 1)].index(
			"Flux_Code")]
	except ValueError:
		print("Flux_Code not found in cwd parents, trying sys.path")
		try:
			# Obtains the location of the Cell_signaling_information folder if it is in sys.path
			path_to_CSI = Path(sys.path[[Path(i).parts[-1] for i in sys.path].index("Flux_Code")])
		except ValueError:
			print("Flux_Code not found in sys.path "
			      "consult 'Errors with setting working directory' in README")
else:
	path_to_CSI = Path(path_to_FC)
os.chdir(path_to_FC)
# _____ Setting the CWD to be Flux_Code END _____

model_file_location = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/Complete_prior.json")
gene_data_file_location = Path("Data/scRNA/Bulk/")
cfp = Path("Data/experimental_alignment_data/recon_1b_t_cells_gene_support/Individual_Mice")
# Sets the day to read the scRNA file
Day = 2

# Dictionary of mice names
mouse_number_filename_dict = {1:"B6-1",2:"B6-2",3:"B6-3",4:"B6-4",5:"TC-5",6:"TC-6",7:"TC-7",8:"TC-8"}

# specifies paths for:
# allowed exchange reactions, measured rates, biomass datasets, experimental alignment results, experimental
# alignment results, and cell sizes for each biological replicate using the keys given above
expt_allow_ex_files_dict = {i:cfp / mouse_number_filename_dict[i]/"Allowed_Exchange_Reactions" for i in mouse_number_filename_dict.keys()}
expt_measured_rates_dict = {i:cfp / mouse_number_filename_dict[i]/"Measured_Rates.csv" for i in mouse_number_filename_dict.keys()}
biomass_rates_dict = {i:cfp / mouse_number_filename_dict[i]/"Biomass_Data.csv" for i in mouse_number_filename_dict.keys()}
output_opt_data_file_location_a_dict = {i:cfp / mouse_number_filename_dict[i]/f"experimental_alignment_result_data_a_Day_{Day}" for i in mouse_number_filename_dict.keys()}
cell_size_file_dict = {i:cfp / mouse_number_filename_dict[i]/f"cell_size.txt" for i in mouse_number_filename_dict.keys()}


list_of_mrna = []
list_of_median_cell_mrna = []
list_of_cell_count = []
weighted_average_denom = 0
error_array = ""
for mouse_number in mouse_number_filename_dict.keys():
    experimental_Allowed_Exchange_file_location = expt_allow_ex_files_dict[mouse_number]
    experimental_measured_rates = expt_measured_rates_dict[mouse_number]
    biomass_rates = biomass_rates_dict[mouse_number]
    output_opt_data_file_location = output_opt_data_file_location_a_dict[mouse_number]
    cell_size_file = cell_size_file_dict[mouse_number]
    
    cell_size_array = np.loadtxt(cell_size_file,delimiter=",")
    mass_of_cell = 10**np.linspace(np.log10(cell_size_array[0]),np.log10(cell_size_array[1]),int(cell_size_array[2]))
    print(f"mass of cell = {mass_of_cell}")
    mass_of_cell *= 10**(-12)
    alpha_list = list(10**(-12)*(1/mass_of_cell))
    if error_array == "":
        error_array = np.zeros((len(mouse_number_filename_dict),len(alpha_list)))
    recon1_rxn_added = Flux_Balance_Model()
    recon1_rxn_added.load_fast_key_json_model(model_file_location)

    print(recon1_rxn_added.test_feasibility())
    recon1_rxn_added.pinch_restricted_exchange_reactions(experimental_Allowed_Exchange_file_location,
                                                         restore_essential=False, true_exchange_signifier="(e)",
                                                         restore_previously_pinched=True)

    # add biomass
    # recon1_rxn_added.print_nonpinched_exchange()

    tol = 1e-7
    br = pd.read_csv(biomass_rates)["value"].to_numpy()[0]
    recon1_rxn_added.update_reaction_bounds("biomass_reaction", br - tol, br + tol)


    for folder in os.listdir(gene_data_file_location):
       if "D" + str(Day) == folder.split("-")[1] and mouse_number_filename_dict[mouse_number].replace("-","") == folder.split("-")[0]:
           features = np.load(gene_data_file_location / folder / "features.npy")
           bulk_gene_mat = np.load(gene_data_file_location / folder / "matrix.npy")
    print(bulk_gene_mat)
    print(np.average(bulk_gene_mat))
    print(np.shape(bulk_gene_mat))
    print(recon1_rxn_added.model_dict.keys())
    print(recon1_rxn_added.model_dict["gene_dict"])
    input()
    RAS_dict = recon1_rxn_added.create_RAS_dict(bulk_gene_mat,features)
    for rxn_name in RAS_dict.keys():
        if RAS_dict[rxn_name] == 0:
            del recon1_rxn_added.model_dict["rxn_dict"][rxn_name]
    

    print(recon1_rxn_added.test_feasibility())

    
    # There did seem to be some risk of not allowed fluxes being just very close to 0
    # Pinching wont fix this
    # Write option to pinch to zero if zero falls between lb and ub
    #flux_modelb.pinch_reactions(10**(-5))
    #alpha_list = 10**np.linspace(20,22,30)
    
    # A549 has 244 pg of protein
    # the fraction of biomass protein is 0.706
    recon1_rxn_added.fit_to_experimental_data(experimental_measured_rates,alpha_list,search_save_path = output_opt_data_file_location)



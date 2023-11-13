import os
import numpy as np
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model


model_file_location = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Models/pkl_models/recon2_2_A549/rxn_added"

experimental_Allowed_Exchange_file_location = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/experimental_alignment_data/recon2_2_A549/Allowed_Exchange_Reactions"

experimental_measured_rates = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/experimental_alignment_data/recon2_2_A549/Measured_Rates.txt"

output_model_file_location = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Models/pkl_models/recon2_2_A549/exp_aligned_alpha_given"

output_opt_data_file_location = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/experimental_alignment_data/recon2_2_A549/experimental_alignment_result_data_alpha_given"

cell_size_file = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/experimental_alignment_data/recon2_2_A549/cell_size.txt"




cell_size_array = np.loadtxt(cell_size_file,delimiter=",")
mass_of_cell = 10**np.linspace(np.log10(cell_size_array[0]),np.log10(cell_size_array[1]),int(cell_size_array[2]))
print(f"mass of cell = {mass_of_cell}")
mass_of_cell *= 10**(-12)
alpha_list = list(10**(-12)*(1/mass_of_cell))

recon2_2_rxn_added = Flux_Balance_Model()
recon2_2_rxn_added.load_pkl_model(model_file_location)

recon2_2_rxn_added.pinch_restricted_exchange_reactions(experimental_Allowed_Exchange_file_location, restore_essential=False)

tol = 1e-7
# There did seem to be some risk of not allowed fluxes being just very close to 0
# Pinching wont fix this
# Write option to pinch to zero if zero falls between lb and ub
#flux_modelb.pinch_reactions(10**(-5))
#alpha_list = 10**np.linspace(20,22,30)

# A549 has 244 pg of protein
# the fraction of biomass protein is 0.706

recon2_2_rxn_added.fit_to_experimental_data(experimental_measured_rates,alpha_list,search_save_file = output_opt_data_file_location)

if recon2_2_rxn_added.test_feasibility():
    recon2_2_rxn_added.save_model_as_pkl(output_model_file_location)
else:
    print("Model Infeasible")


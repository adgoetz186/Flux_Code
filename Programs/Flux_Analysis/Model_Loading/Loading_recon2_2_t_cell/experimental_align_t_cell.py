import os
import numpy as np
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
from pathlib import Path


model_file_location = Path("/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Models/pkl_models/recon2_2_t_cell/rxn_added.pkl")

mouse_number_filename_dict = {1:"B6-1",2:"B6-2",3:"B6-3",4:"B6-4",5:"TC-5",6:"TC-6",7:"TC-7",8:"TC-8"}
error_array = ""
for mouse_number in mouse_number_filename_dict.keys():

    experimental_Allowed_Exchange_file_location = Path(f"/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/experimental_alignment_data/recon2_2_t_cell/Individual_Mice/{mouse_number_filename_dict[mouse_number]}/Allowed_Exchange_Reactions")
    
    experimental_measured_rates = Path(f"/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/experimental_alignment_data/recon2_2_t_cell/{mouse_number_filename_dict[mouse_number]}/Measured_Rates.csv")
    
    output_model_file_location = Path(f"/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Models/pkl_models/recon2_2_t_cell/exp_aligned/{mouse_number_filename_dict[mouse_number]}_exp_aligned")
    
    output_opt_data_file_location = Path(f"/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/experimental_alignment_data/recon2_2_t_cell/Individual_Mice/{mouse_number_filename_dict[mouse_number]}/experimental_alignment_result_data")
    
    cell_size_file = Path(f"/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/experimental_alignment_data/recon2_2_t_cell/Individual_Mice/{mouse_number_filename_dict[mouse_number]}/cell_size.txt")

    
    cell_size_array = np.loadtxt(cell_size_file,delimiter=",")
    mass_of_cell = 10**np.linspace(np.log10(cell_size_array[0]),np.log10(cell_size_array[1]),int(cell_size_array[2]))
    print(f"mass of cell = {mass_of_cell}")
    mass_of_cell *= 10**(-12)
    alpha_list = list(10**(-12)*(1/mass_of_cell))
    if error_array == "":
        error_array = np.zeros((len(mouse_number_filename_dict.keys()),len(alpha_list)))
    
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
    
    recon2_2_rxn_added.fit_to_experimental_data(experimental_measured_rates,alpha_list,search_save_path = output_opt_data_file_location)
    
    
    
    if recon2_2_rxn_added.test_feasibility():
        recon2_2_rxn_added.save_model_as_pkl(output_model_file_location)
    else:
        print("Model Infeasible")


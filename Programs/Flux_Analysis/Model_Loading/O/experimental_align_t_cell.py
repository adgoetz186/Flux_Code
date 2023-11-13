import os
import numpy as np
import matplotlib.pyplot as plt
import pickle
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

model_file_location = Path("Data/Models/pkl_models/recon1b_t_cell/rxn_added.pkl")

mouse_number_filename_dict = {1:"B6-1",2:"B6-2",3:"B6-3",4:"B6-4",5:"TC-5",6:"TC-6",7:"TC-7",8:"TC-8"}
error_array = ""
for mouse_number in mouse_number_filename_dict.keys():
    experimental_Allowed_Exchange_file_location = Path(f"Data/experimental_alignment_data/recon1b_t_cell/Individual_Mice/{mouse_number_filename_dict[mouse_number]}/Allowed_Exchange_Reactions")
    
    experimental_measured_rates = Path(f"Data/experimental_alignment_data/recon1_b_t_cell/Individual_Mice/{mouse_number_filename_dict[mouse_number]}/Measured_Rates.csv")
    
    output_opt_data_file_location = Path(f"Data/experimental_alignment_data/recon1_b_t_cell/Individual_Mice/{mouse_number_filename_dict[mouse_number]}/experimental_alignment_result_data_a")
    
    cell_size_file = Path(f"Data/experimental_alignment_data/recon1_b_t_cell/Individual_Mice/{mouse_number_filename_dict[mouse_number]}/cell_size.txt")
    
    cell_size_array = np.loadtxt(cell_size_file,delimiter=",")
    mass_of_cell = 10**np.linspace(np.log10(cell_size_array[0]),np.log10(cell_size_array[1]),int(cell_size_array[2]))
    print(f"mass of cell = {mass_of_cell}")
    mass_of_cell *= 10**(-12)
    alpha_list = list(10**(-12)*(1/mass_of_cell))
    if error_array == "":
        error_array = np.zeros((len(mouse_number_filename_dict),len(alpha_list)))
    recon1_rxn_added = Flux_Balance_Model()
    recon1_rxn_added.load_pkl_model(model_file_location)
    print(recon1_rxn_added.test_feasibility())
    input()
    recon1_rxn_added.pinch_restricted_exchange_reactions(experimental_Allowed_Exchange_file_location, restore_essential=False,true_exchange_signifier="(e)",unpinch_allowed = True)
    recon1_rxn_added.dicts_to_mats()
    print(recon1_rxn_added.test_feasibility())
    input()
    tol = 1e-7
    # There did seem to be some risk of not allowed fluxes being just very close to 0
    # Pinching wont fix this
    # Write option to pinch to zero if zero falls between lb and ub
    #flux_modelb.pinch_reactions(10**(-5))
    #alpha_list = 10**np.linspace(20,22,30)
    
    # A549 has 244 pg of protein
    # the fraction of biomass protein is 0.706
    
    recon1_rxn_added.fit_to_experimental_data(experimental_measured_rates,alpha_list,search_save_path = output_opt_data_file_location)

    with open(output_opt_data_file_location.parent/(output_opt_data_file_location.stem+".pkl"), "rb") as readfile:
        model_dict = pickle.load(readfile)
    X_array = model_dict["model_alignment"]
    exp = model_dict["experimental_fluxes"]
    ind = model_dict["experimental_flux_index"]

    labels = [model_dict["saved_rxn_name_order"][i] for i in ind]
    
    # remove valine flux from consideration in error term
    del ind[labels.index('EX_val_L(e)')]

    for alpha_index in range(np.size(alpha_list)):
        alpha = alpha_list[alpha_index]
        expt = np.array([exp[i] for i in ind])
        mdl = np.array([X_array[alpha_index][i] for i in ind])
        error_array[mouse_number - 1, alpha_index] = np.sum((mdl - expt * alpha) ** 2 / (expt * alpha) ** 2)
print(error_array)
error_vec = np.average(error_array,axis=0)
min_ind = np.argmin(error_vec)
plt.scatter(mass_of_cell[min_ind]*10**12,error_vec[min_ind])
plt.plot(mass_of_cell*10**12,error_vec)
plt.show()
for i in error_array:
    plt.plot(mass_of_cell*10**12,i)
plt.scatter(mass_of_cell[min_ind]*10**12,error_vec[min_ind])
plt.show()
print("done")
mass_of_cell = mass_of_cell[min_ind]



for mouse_number in mouse_number_filename_dict.keys():
    
    experimental_Allowed_Exchange_file_location = Path(f"/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/experimental_alignment_data/recon1_b_t_cell/Individual_Mice/{mouse_number_filename_dict[mouse_number]}/Allowed_Exchange_Reactions")
    
    experimental_measured_rates = Path(f"/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/experimental_alignment_data/recon1_b_t_cell/Individual_Mice/{mouse_number_filename_dict[mouse_number]}/Measured_Rates.csv")
    
    output_model_file_location = Path(f"/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Models/pkl_models/recon1_b_t_cell/exp_aligned/{mouse_number_filename_dict[mouse_number]}_exp_aligned")
    
    output_opt_data_file_location = Path(f"/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/experimental_alignment_data/recon1_b_t_cell/Individual_Mice/{mouse_number_filename_dict[mouse_number]}/experimental_alignment_result_data_b")
    
    print(f"mass of cell = {mass_of_cell}")
    alpha_list = [10 ** (-12) * (1 / mass_of_cell)]
    if error_array == "":
        error_array = np.zeros((len(mouse_number_filename_dict), len(alpha_list)))
    recon1_rxn_added = Flux_Balance_Model()
    recon1_rxn_added.load_pkl_model(model_file_location)
    recon1_rxn_added.model_name = mouse_number_filename_dict[mouse_number]

    recon1_rxn_added.pinch_restricted_exchange_reactions(experimental_Allowed_Exchange_file_location,
                                                           restore_essential=False, true_exchange_signifier="(e)")
    
    tol = 1e-7
    # There did seem to be some risk of not allowed fluxes being just very close to 0
    # Pinching wont fix this
    # Write option to pinch to zero if zero falls between lb and ub
    # flux_modelb.pinch_reactions(10**(-5))
    # alpha_list = 10**np.linspace(20,22,30)
    
    # A549 has 244 pg of protein
    # the fraction of biomass protein is 0.706
    
    recon1_rxn_added.fit_to_experimental_data(experimental_measured_rates, alpha_list,
                                                search_save_path=output_opt_data_file_location)
    
    if recon1_rxn_added.test_feasibility():
        recon1_rxn_added.save_model_as_pkl(output_model_file_location)
    else:
        print("Model Infeasible")


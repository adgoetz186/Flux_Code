import os
import pickle
import gurobipy as gp
import copy
import copy as cp
import scipy.optimize as so
from gurobipy import GRB
from pathlib import Path
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import numpy as np
from matplotlib import pyplot as plt

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

mouse_number_filename_dict = {1:"B6-1",2:"B6-2",3:"B6-3",4:"B6-4",5:"TC-5",6:"TC-6",7:"TC-7",8:"TC-8"}
#mouse_number_filename_dict = {1:"B6-1",5:"TC-5"}
rel_error_all_mice = []

for mouse_number in mouse_number_filename_dict.keys():
    file_location = Path(f"Data/experimental_alignment_data/recon_1b_t_cells_gene_support/Individual_Mice/{mouse_number_filename_dict[mouse_number]}/experimental_alignment_result_data_a_Day_2.pkl")
    with open(file_location, "rb") as readfile:
        model_dict = pickle.load(readfile)
    model_file_location = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/exp_aligned_Day_2/")
    model_file_name = f"{mouse_number_filename_dict[mouse_number]}.json"
    
    X_array = model_dict["model_alignment"]
    best_ind = model_dict["best_alpha_index"]
    alpha_list = model_dict["alpha_list"]
    
    exp = model_dict["experimental_fluxes"]
    Obj_array = model_dict["model_alignment_performances"]
    ind = model_dict["experimental_flux_index"]
    recon2_2_rxn_added = Flux_Balance_Model()
    recon2_2_rxn_added.load_fast_key_json_model(model_file_location/model_file_name)
    lb,ub,S,b,ordered_rxn,ordered_met = recon2_2_rxn_added.dicts_to_mats()

    labels = [ordered_rxn[i] for i in ind]


    del ind[labels.index('EX_val_L(e)')]


    if rel_error_all_mice == []:
        rel_error_all_mice = np.zeros((8, np.size(alpha_list)))
    
    for alpha_index in range(np.size(alpha_list)):
        alpha = alpha_list[alpha_index]
        expt = np.array([exp[i] for i in ind])
        mdl = np.array([X_array[alpha_index][i] for i in ind])
        rel_error_all_mice[mouse_number-1,alpha_index] = np.sum((mdl-expt*alpha)**2/(expt*alpha)**2)
    print(rel_error_all_mice[mouse_number-1,np.argmin(rel_error_all_mice[mouse_number-1])])
   
rel_error = np.average(rel_error_all_mice,axis=0)
rel_std = np.var(rel_error_all_mice,axis=0)**.5
best_ind = np.argmin(rel_error)
print("start")
for i in range(8):
    
    print(i, rel_error_all_mice[i,best_ind])

print(1/np.array(alpha_list[best_ind]))
print(rel_error[best_ind])
print(rel_std[best_ind])
plt.scatter(1/np.array(alpha_list[best_ind]),rel_error[best_ind])
plt.plot(1/np.array(alpha_list),rel_error)
#plt.hlines(np.min(Obj_array),0,1000)
#plt.hlines(1.05*np.min(Obj_array), 0, 1000)
#plt.hlines(1.10 * np.min(Obj_array), 0, 1000)
plt.vlines(45,0,100,color = 'black')
plt.vlines(79, 0, 100,color = 'black',linestyles="--")
plt.vlines(125, 0, 100,color = 'black')
plt.fill_between([45,125],[2,2],color = "grey",alpha = 0.3)
plt.xlabel("Cell Size (pg)")
plt.ylabel("Error")
plt.xscale("log")
#plt.ylim([0,2])
plt.xlim([0,1000])
plt.title(f"All mice\nProblematic Valine Terms Removed")
plt.show()

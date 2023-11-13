import os
import pickle
import gurobipy as gp
import copy
import copy as cp
import scipy.optimize as so
from gurobipy import GRB
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import numpy as np
from pathlib import Path
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
    file_location = Path(f"Data/experimental_alignment_data/recon_1b_t_cells/Individual_Mice/{mouse_number_filename_dict[mouse_number]}/experimental_alignment_result_data_b.pkl")
    with open(file_location, "rb") as readfile:
        model_dict = pickle.load(readfile)
    model_file_location = Path("Data/Models/pkl_models/recon_1b_t_cells/exp_aligned/")
    model_file_name = f"{mouse_number_filename_dict[mouse_number]}_exp_aligned.pkl"
    recon2_2_rxn_added = Flux_Balance_Model()
    recon2_2_rxn_added.load_pkl_model(model_file_location/model_file_name)
    X_array = model_dict["model_alignment"]
    best_ind = model_dict["best_alpha_index"]
    alpha_list = model_dict["alpha_list"]
    ordered_rxn = model_dict["saved_rxn_name_order"]
    exp = model_dict["experimental_fluxes"]
    Obj_array = model_dict["model_alignment_performances"]
    ind = model_dict["experimental_flux_index"]
    print(1/np.array(alpha_list))
    labels = [model_dict["saved_rxn_name_order"][i] for i in ind]
    

    X = X_array[best_ind]
    alpha = alpha_list[best_ind]

    plot_experiment = [exp[i] for i in ind]
    plot_model = [X[i] / alpha for i in ind]
    Xlb_mdl = [recon2_2_rxn_added.rxn_dict[ordered_rxn[i]]['lb'] / alpha for i in ind]
    Xub_mdl = [recon2_2_rxn_added.rxn_dict[ordered_rxn[i]]['ub'] / alpha for i in ind]

    
    ind_val_removed = cp.deepcopy(ind)
    del ind_val_removed[labels.index('EX_val_L(e)')]
    rel_error = np.zeros(np.size(alpha_list))
    for alpha_index in range(np.size(alpha_list)):
        
        alpha = alpha_list[alpha_index]
        expt = np.array([exp[i] for i in ind_val_removed])
        mdl = np.array([X_array[alpha_index][i] for i in ind_val_removed])
        rel_error[alpha_index] = np.sum((mdl-expt*alpha)**2/(expt*alpha)**2)
        print(expt*alpha,mdl,rel_error[alpha_index])
    print(rel_error)
    best_ind = np.argmin(rel_error)


    
    plt.scatter(1/np.array(alpha_list[best_ind]),rel_error[best_ind])
    plt.plot(1/np.array(alpha_list),rel_error,label = "V penalty removed")
    print(1/np.array(alpha_list[best_ind]))
    plt.scatter(1 / np.array(alpha_list[best_ind]), Obj_array[best_ind])
    plt.plot(1/np.array(alpha_list),Obj_array,label = "Optimized Objective")
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
    plt.title(f"{mouse_number_filename_dict[mouse_number]}")
    plt.legend()
    plt.show()
    recon2_2_rxn_added = Flux_Balance_Model()
    recon2_2_rxn_added.load_pkl_model(model_file_location/model_file_name)
    
    recon2_2_rxn_added.reaction_info("EX_val_L(e)")
    #input()
    exchange_reaction_dict = recon2_2_rxn_added.get_exchanged_reaction_info()
    
    #exchange_reaction_index_list = [recon2_2_rxn_added.reaction_names.index(exchange_reaction_dict[i]) for i in list(exchange_reaction_dict.keys())]
    
    #select_reaction_names = ["ARGNm","AGMTm","UREAtm","ENO","AKGDm","CSm","FUM","FUMm","ICDHxm"]
    #select_reaction_index = [recon2_2_rxn_added.reaction_names.index(i) for i in select_reaction_names]
    
    #select_reactions = np.zeros((len(select_reaction_index),np.size(alpha_list)))
    #exchange_reactions = np.zeros((len(exchange_reaction_index_list),np.size(alpha_list)))

    X = X_array[best_ind]
    alpha = alpha_list[best_ind]
    #labels = [recon2_2_rxn_added.reaction_names[i] for i in ind]
    plot_experiment = [exp[i] for i in ind]
    plot_model = [X[i]/alpha for i in ind]
    Xlb_mdl = [recon2_2_rxn_added.rxn_dict[ordered_rxn[i]]['lb'] / alpha for i in ind]
    Xub_mdl = [recon2_2_rxn_added.rxn_dict[ordered_rxn[i]]['ub'] / alpha for i in ind]
    
    etplot_experiment = [alpha*exp[i] for i in ind]
    etplot_model = [X[i] for i in ind]
    
    error_term = 0
    for i in range(len(etplot_model)):
        error_term+= (etplot_experiment[i]-etplot_model[i])**2/etplot_experiment[i]**2
    print(error_term)
    x = np.arange(len(labels))  # the label locations
    
    width = 0.25
    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width / 2, plot_experiment, width/2, label='Expt')
    rects1 = ax.bar(x - width / 4, plot_model, width/2, label='Mdl')
    rects1 = ax.bar(x + width / 4, Xub_mdl, width/2, label='ub')
    rects2 = ax.bar(x + width / 2, Xlb_mdl, width/2, label='lb')
    ax.set_ylabel('Excretion Flux (fmol/(cell * hr))')
    ax.set_title(f'{mouse_number_filename_dict[mouse_number]}\nModel and Experimental Agreement (ornithine intake)')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90)
    ax.legend()
    # ax.bar_label(rects1, padding=3)
    # ax.bar_label(rects2, padding=3)
    
    fig.tight_layout()
    plt.hlines(0, -1, 23, color="black")
    plt.show()
    
    
    
    #counter = 0
    #for i in select_reaction_index:
    #    plt.plot(np.ones_like(alpha_list)/alpha_list,select_reactions[counter])
    #    plt.title(f"{recon2_2_rxn_added.reaction_names[i]}")
    #    plt.xlabel("Cell Size (pg)")
    #    plt.xscale("log")
    #    plt.ylabel("Model Flux")
    #    plt.show()
    #    counter += 1
    
    #counter = 0
    #for i in exchange_reaction_index_list:
    #    if np.sum(exchange_reactions[counter]) != 0:
    #        plt.plot(np.ones_like(alpha_list)/alpha_list,exchange_reactions[counter])
    #        plt.title(f"{recon2_2_rxn_added.reaction_names[i]}")
    #        plt.xlabel("Cell Size (pg)")
    #        plt.xscale("log")
    #        plt.ylabel("Model Flux")
    #        plt.show()
    #    counter += 1
    #print("done")
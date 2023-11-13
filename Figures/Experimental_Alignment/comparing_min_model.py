import os
import pickle
import gurobipy as gp
import copy
import copy as cp
import scipy.optimize as so
from gurobipy import GRB
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import numpy as np
from matplotlib import pyplot as plt

os.chdir("../../../")
mouse_number_filename_dict = {1: "B6-1", 2: "B6-2", 3: "B6-3", 4: "B6-4", 5: "TC-5", 6: "TC-6", 7: "TC-7", 8: "TC-8"}
# mouse_number_filename_dict = {1:"B6-1",5:"TC-5"}
rel_error_all_mice = []
for mouse_number in mouse_number_filename_dict.keys():
	file_location = f"/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/experimental_alignment_data/recon1_t_cell/Individual_mice/{mouse_number_filename_dict[mouse_number]}/experimental_alignment_result_data_b.pkl"
	with open(file_location, "rb") as readfile:
		model_dict = pickle.load(readfile)
	model_file_location = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Models/pkl_models/recon1_t_cell/exp_aligned/"
	model_file_name = f"{mouse_number_filename_dict[mouse_number]}_exp_aligned.pkl"
	recon2_2_rxn_added = Flux_Balance_Model()
	recon2_2_rxn_added.load_pkl_model(f"{model_file_location}{model_file_name}")
	
	ind = model_dict["experimental_flux_index"]
	size_big = len(recon2_2_rxn_added.lb)
	Xlb_mdl_big = [recon2_2_rxn_added.lb[i] for i in ind]
	Xub_mdl_big = [recon2_2_rxn_added.ub[i] for i in ind]

	labels = [recon2_2_rxn_added.reaction_names[i] for i in ind]
	
	model_file_location = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Models/pkl_models/recon1_t_cell/min_models/"
	model_file_name = f"{mouse_number_filename_dict[mouse_number]}_min_model.pkl"
	recon2_2_rxn_added = Flux_Balance_Model()
	recon2_2_rxn_added.load_pkl_model(f"{model_file_location}{model_file_name}")
	
	print(recon2_2_rxn_added.reaction_names)
	small_ind = [recon2_2_rxn_added.reaction_names.index(i) for i in labels]
	
	size_small = len(recon2_2_rxn_added.lb)
	Xlb_mdl = [recon2_2_rxn_added.lb[i] for i in small_ind]
	Xub_mdl = [recon2_2_rxn_added.ub[i] for i in small_ind]
	
	print(np.max(np.abs(np.array(Xlb_mdl_big) - np.array(Xlb_mdl))))
	print(np.max(np.abs(np.array(Xub_mdl_big) - np.array(Xub_mdl))))
	print(size_big,size_small)
	
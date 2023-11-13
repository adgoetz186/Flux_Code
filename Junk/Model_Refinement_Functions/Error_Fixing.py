import os
import pickle
import numpy as np
import copy as cp

os.chdir("../..")
all_essential = True
nonessential_value = 0
nonessential_list = []
essential_flux_matrix = np.loadtxt("Data/Intermediate/Essential_Flux_Matrix/EFM_A549_recon2_2d.txt",delimiter=",")
for row in essential_flux_matrix:
	if nonessential_value == 9:
		flux_model = pickle.load(open("Data/Intermediate/Experimentally_Aligned_Model/A549_recon2_2d.mdl", "rb"))
		for essential_matrix_index in range(len(row)):
			#if the reaction is essential then lb and ub will be the same, if it is not then lb and ub will be 0
			flux_model.update_reaction_bounds(essential_matrix_index,row[essential_matrix_index]*flux_model.lb[essential_matrix_index],row[essential_matrix_index]*flux_model.ub[essential_matrix_index])
		if not flux_model.test_feasibility():
			all_essential = False
			if nonessential_value not in nonessential_list:
				nonessential_list.append(nonessential_value)
			print("not feasible enough")
		for essential_matrix_index in range(len(row)):
			if essential_matrix_index%1000 ==0:
				print(essential_matrix_index/len(row))
			if row[essential_matrix_index] == 1:
				if not (flux_model.lb[essential_matrix_index]>0 or flux_model.ub[essential_matrix_index] < 0):
					save_ub = cp.deepcopy(flux_model.ub[essential_matrix_index])
					save_lb = cp.deepcopy(flux_model.lb[essential_matrix_index])
					flux_model.update_reaction_bounds(essential_matrix_index,0,0)
					if flux_model.test_feasibility():
						all_essential = False
						if nonessential_value not in nonessential_list:
							nonessential_list.append(nonessential_value)
						print(save_lb, save_ub)
						print(flux_model.reaction_info(essential_matrix_index))
						print("too feasible")
					else:
						flux_model.update_reaction_bounds(essential_matrix_index, save_lb, save_ub)
		print(nonessential_list)
	nonessential_value+=1
print(all_essential)
print(nonessential_list)

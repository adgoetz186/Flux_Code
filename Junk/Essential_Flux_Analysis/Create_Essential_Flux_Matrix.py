import os
import pickle
import numpy as np
import copy as cp

os.chdir("../..")
with open("Data/Intermediate/Essential_Flux_Lists/exhaust_EF_A549_recon2_2_GLUN_and_SERtm_added.txt","r") as eff:
	flux_lists_string = eff.readlines()
flux_lists = [eval(i) for i in flux_lists_string]
#print(flux_lists)

flux_model = pickle.load(open("Data/Intermediate/Experimentally_Aligned_Model/A549_recon2_2_GLUN_and_SERtm_added.mdl", "rb"))

essential_flux_matrix = np.zeros((len(flux_lists),len(flux_model.reaction_names)))
for i in range(len(flux_lists)):
	for j in range(len(flux_lists[i])):
		essential_flux_matrix[i,flux_model.reaction_names.index(flux_lists[i][j])] = 1
print(essential_flux_matrix)

np.savetxt("Data/Intermediate/Essential_Flux_Matrix/EFM_A549_recon2_2_GLUN_and_SERtm_added.csv",essential_flux_matrix,delimiter=",")


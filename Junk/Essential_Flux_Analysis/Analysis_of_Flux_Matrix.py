import os
import pickle
import numpy as np
import copy as cp

os.chdir("../..")
flux_model = pickle.load(open("Data/Intermediate/Experimentally_Aligned_Model/A549_recon2_2_GLUN_and_SERtm_added.mdl", "rb"))
#essential_flux_matrix = np.loadtxt("Data/Intermediate/Essential_Flux_Matrix/EFM_A549_recon2_2_GLUN_and_SERtm_added.csv",delimiter=",")
essential_flux_model = pickle.load(open("Data/Intermediate/Essential_Flux_Model/A549_recon2_2_GLUN_and_SERtm_added.mdl", "rb"))
#essential_flux_count = np.average(essential_flux_matrix,axis=1)*len(flux_model.reaction_names)

#print(np.shape(essential_flux_count))
#print(np.average(essential_flux_count))
#print(np.median(essential_flux_count))
#print(np.min(essential_flux_count))
#print(np.max(essential_flux_count))
print(np.shape(essential_flux_model.S))
essential_matrix = essential_flux_model.S.todense()
essential_flux_model.purge_metabolites()
print(np.shape(essential_flux_model.S))
input()

full_matrix = flux_model.S.todense()
print(flux_model.pinched_reactions)
print(np.shape(essential_matrix))
print(np.shape(essential_matrix),np.shape(full_matrix))
cnt = 0
for i in np.sum(np.abs(essential_matrix),axis = 1):
	if i == 0:
		cnt +=1
		print(i)
print(cnt)
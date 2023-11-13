import os
import pickle
import numpy as np
import copy as cp
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Save_And_Load import save_object

os.chdir("../..")
flux_model = pickle.load(open("Data/Intermediate/Experimentally_Aligned_Model/A549_recon2_2_GLUN_and_SERtm_added.mdl", "rb"))
essential_flux_matrix = np.loadtxt("Data/Intermediate/Essential_Flux_Matrix/EFM_A549_recon2_2_GLUN_and_SERtm_added.csv",delimiter=",")
essential_flux_count = np.average(essential_flux_matrix,axis=1)*len(flux_model.reaction_names)


essential_flux_scores = np.average(essential_flux_matrix,axis=0)
least_essential_flux_indecies = np.argsort(essential_flux_scores)


most_essential_flux_indecies = np.flip(least_essential_flux_indecies)

saved_ub = cp.deepcopy(flux_model.ub)
saved_lb = cp.deepcopy(flux_model.lb)

for i in most_essential_flux_indecies:
	if (flux_model.ub[i] < 0 or flux_model.lb[i] > 0):
		continue
	else:
		flux_model.update_reaction_bounds(i,0,0)
essential_flux_main = []
for i in most_essential_flux_indecies:
	if (flux_model.ub[i] < 0 or flux_model.lb[i] > 0):
		essential_flux_main.append(i)
		continue
	else:
		flux_model.update_reaction_bounds(i,saved_lb[i],saved_ub[i])
		if flux_model.test_feasibility():
			essential_flux_main.append(i)
			break
		else:
			essential_flux_main.append(i)
count = 0
most_essential_flux_names = [flux_model.reaction_names[i] for i in most_essential_flux_indecies]
essential_flux_main_names = [flux_model.reaction_names[i] for i in essential_flux_main]
not_essential_fluxes = list(set(most_essential_flux_names) - set(essential_flux_main_names))
print(flux_model.test_feasibility())
flux_model.del_reaction(not_essential_fluxes)
flux_model.purge_metabolites()
flux_limits = flux_model.fva()
for i in range(len(flux_limits)):
	flux_model.lb[i] = flux_limits[i][0]
	flux_model.ub[i] = flux_limits[i][1]
save_object(flux_model, "Data/Intermediate/Essential_Flux_Model/A549_recon2_2_GLUN_and_SERtm_added.mdl")


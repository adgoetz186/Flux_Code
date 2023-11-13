import os
import pickle
import numpy as np
import copy as cp
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Save_And_Load import save_object

os.chdir("../../..")
flux_model = pickle.load(open("Data/Intermediate/Essential_Flux_Model/A549_recon2_2_GLUN_and_SERtm_added.mdl", "rb"))

count = 0
print(np.shape(flux_model.S))
print(len(flux_model.lb))
print(len(flux_model.b))
reactions_to_kill = []
for i in range(len(flux_model.lb)):
	if abs(flux_model.ub[i]- flux_model.lb[i]) < 1e-9:
		print(flux_model.reaction_names[i], flux_model.lb[i], flux_model.ub[i])
		count+=1
		reactions_to_kill.append(flux_model.reaction_names[i])
flux_model.del_reaction(reactions_to_kill)
flux_model.purge_metabolites()
print(flux_model.test_feasibility())
print(np.shape(flux_model.S))
print(count)
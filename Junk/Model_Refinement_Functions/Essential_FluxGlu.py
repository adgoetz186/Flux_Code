import numpy as np
import os
import pickle
import copy
import time

from Flux_Code.Programs.python_functions import File_handling

#Get the single cache run to return an array with size equal to the number of max reactions
#Should have a 1 if reaction at the index is essential and a 0 if it isnt
#It might be useful to have each worker do more than 1 run in which case the allowed values of the return array will be intagers and not binary





#Sets cwd as the main file directory


File_handling.up_n_files(3)
flux_model = pickle.load(open("Data/Intermediate/Raw_MDL_Models/High_Oxygen_Ecoli_Positive_Flux_10_percent.mdl", "rb"))
initial_flux_counts = dict(zip(flux_model.reaction_names,[0 for i in range(len(flux_model.reaction_names))]))
initial_flux_counts_zeros = copy.deepcopy(initial_flux_counts)
biomass_obj = flux_model.objective_function_generator(["Biomass_Ecoli_core_w_GAM"],[1])
max_val = flux_model.find_objective_value(biomass_obj)

print(max_val)
input()

flux_obj = flux_model.objective_function_generator(np.arange(np.size(flux_model.ub)),-1*np.ones_like(flux_model.ub))
min_flux = flux_model.find_objective_value(flux_obj)

#3 values are important
#Flush count, number of runs before the list is cleared and tallied
#Single monte carlo approximation, the first level approximation
#Multi Monte carlo approximation



saved_model = copy.deepcopy(pickle.load(open("Data/Intermediate/Raw_MDL_Models/High_Oxygen_Ecoli_Positive_Flux_10_percent.mdl", "rb")))
Start = time.time()
flush_count = 500
single_monte_carlo = 100
multi_monte_carlo = 10
error_count = 0
meta_mc_storage = np.zeros((multi_monte_carlo,len(initial_flux_counts.keys())))
for mmc in range(multi_monte_carlo):
    initial_flux_counts = copy.deepcopy(initial_flux_counts_zeros)
    np.array(list(initial_flux_counts.values()))
    this_feels_wrong = []
    for smc in range(single_monte_carlo):
        flux_model = copy.deepcopy(saved_model)
        efl = flux_model.find_essential_fluxes(extra_constraints=[f"np.array({biomass_obj.tolist()}) @ react_flux >= {max_val*.9}",f"np.array({flux_obj.tolist()}) @ react_flux >= {min_flux*1.1}"])
        for j in efl:
            this_feels_wrong.append(j)
    for key in this_feels_wrong:
        initial_flux_counts[key]+=1
    for key in initial_flux_counts.keys():
        initial_flux_counts[key]/=single_monte_carlo
    meta_mc_storage[mmc] = np.array(list(initial_flux_counts.values()))
    print(mmc,1)
print(error_count/1000)
print(np.mean(meta_mc_storage,axis=0))
print(np.sum(np.mean(meta_mc_storage,axis=0))/np.size(flux_model.lb))
sv = np.mean(meta_mc_storage,axis=0)






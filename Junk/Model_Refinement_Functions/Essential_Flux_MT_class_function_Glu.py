import json

import numpy as np
import os
import pickle
import copy
from functools import partial
import multiprocessing as mp
import time

from Flux_Code.Programs.python_functions import File_handling

#Get the single cache run to return an array with size equal to the number of max reactions
#Should have a 1 if reaction at the index is essential and a 0 if it isnt
#It might be useful to have each worker do more than 1 run in which case the allowed values of the return array will be intagers and not binary

def single_cache_run(saved_model,repeat,no,**kwargs):
    initial_flux_counts = dict(zip(saved_model.reaction_names, [0 for i in range(len(saved_model.reaction_names))]))
    zero_list = np.zeros_like(saved_model.lb)
    for i in range(repeat):
        flux_model = copy.deepcopy(saved_model)
        zero_list+=np.array(list(flux_model.find_essential_fluxes(dict_count=initial_flux_counts,seed=no,**kwargs).values()))
    return zero_list
#Sets cwd as the main file directory

if __name__ == '__main__':
    File_handling.up_n_files(3)
    flux_model = pickle.load(open("Data/Intermediate/Raw_MDL_Models/High_Oxygen_Ecoli_Positive_Flux_10_percent.mdl", "rb"))
    biomass_obj = flux_model.objective_function_generator(["Biomass_Ecoli_core_w_GAM"],[1])
    max_val = flux_model.find_objective_value(biomass_obj)
    flux_obj = flux_model.objective_function_generator(np.arange(np.size(flux_model.ub)),-1 * np.ones_like(flux_model.ub))
    min_flux = flux_model.find_objective_value(flux_obj)


    #3 values are important
    #Flush count, number of runs before the list is cleared and tallied
    #Single monte carlo approximation, the first level approximation
    #Multi Monte carlo approximation

    # Important note: Dont update reactions on the fly, only update after a full run
    # You dont want to pinch a reaction and change the index values and make the biomass objective refer to the wrong
    # flux

    Start = time.time()
    flush_count = 16
    single_monte_carlo = 2
    multi_monte_carlo = 10
    meta_mc_storage = np.zeros((multi_monte_carlo,len(flux_model.lb)))
    for mmc in range(multi_monte_carlo):
        pool = mp.Pool(processes=mp.cpu_count())
        repeat = 1
        func = partial(single_cache_run,flux_model,repeat,extra_constraints=[f"np.array({biomass_obj.tolist()}) @ react_flux >= {max_val*.9}",f"np.array({flux_obj.tolist()}) @ react_flux >= {min_flux*1.1}"])
        ResB_P = pool.map(func, [i for i in range(flush_count)])
        pool.close()
        sum_of_output = np.zeros_like(ResB_P[0])
        for i in ResB_P:
            sum_of_output+= i
        print(sum_of_output/(flush_count*repeat))
        meta_mc_storage[mmc] = sum_of_output/(repeat*flush_count)
        print((time.time()-Start)/(mmc+1))
    essential_flux_dictionary = dict(zip(flux_model.reaction_names,np.mean(meta_mc_storage, axis=0)))
    
    with open("Data/Intermediate/Essential_Flux_Values/EF_High_Oxygen_Ecoli_Positive_Flux_10_percent.json","w") as EFO:
        json.dump(essential_flux_dictionary,EFO)



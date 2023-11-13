import numpy as np
import os
import json
import pickle
import gurobipy as gp
import copy
import docplex.mp.model as cpx
from gurobipy import GRB
from functools import partial
import multiprocessing as mp
import time

import Flux_Code.Programs.gene_rules_regex_testing


#Get the single cache run to return an array with size equal to the number of max reactions
#Should have a 1 if reaction at the index is essential and a 0 if it isnt
#It might be useful to have each worker do more than 1 run in which case the allowed values of the return array will be intagers and not binary

def single_cache_run(saved_model,repeat,no,**kwargs):
    initial_flux_counts = dict(zip(saved_model.reaction_names, [0 for i in range(len(saved_model.reaction_names))]))
    zero_list = np.zeros_like(saved_model.lb)
    for i in range(repeat):
        flux_model = copy.deepcopy(saved_model)
        zero_list+=np.array(list(flux_model.find_essential_fluxes_fast(dict_count=initial_flux_counts,seed=no,**kwargs).values()))
    return zero_list
#Sets cwd as the main file directory

if __name__ == '__main__':
    #File_handling.up_n_files(3)
    os.chdir("..")
    flux_model = pickle.load(open("Data/Intermediate/Experimentally_Aligned_Model/A549_recon2_2d.mdl", "rb"))
    initial_flux_counts = dict(zip(flux_model.reaction_names,[0 for i in range(len(flux_model.reaction_names))]))
    initial_flux_counts_zeros = copy.deepcopy(initial_flux_counts)


    print(flux_model.reaction_names)
    saved_model = pickle.load(open("Data/Intermediate/Experimentally_Aligned_Model/A549_recon2_2d.mdl", "rb"))
    Start = time.time()
    flush_count = 8
    single_monte_carlo = 2
    multi_monte_carlo = 10
    meta_mc_storage = np.zeros((multi_monte_carlo,len(initial_flux_counts.keys())))
    #You are using seed wrong I think, add time values and uniquely give all values seeds
    for mmc in range(multi_monte_carlo):
        print(mmc)
        print([mmc * flush_count + i for i in range(flush_count)])
        initial_flux_counts = copy.deepcopy(initial_flux_counts_zeros)
        pool = mp.Pool(processes=mp.cpu_count())
        repeat = 1
        func = partial(single_cache_run,saved_model,repeat)
        ResB_P = pool.map(func, [(mmc*flush_count + i)*time.time() for i in range(flush_count)])
        pool.close()
        sum_of_output = np.zeros_like(ResB_P[0])
        for i in ResB_P:
            sum_of_output+= i
        print(sum_of_output/(flush_count*repeat))
        meta_mc_storage[mmc] = sum_of_output/(repeat*flush_count)
        print((time.time()-Start)/(mmc+1))
    essential_flux_dictionary = dict(zip(flux_model.reaction_names, np.mean(meta_mc_storage, axis=0)))

    with open("Data/Intermediate/Essential_Flux_Values/EF_A549_recon2_2d.json", "w") as EFO:
        json.dump(essential_flux_dictionary, EFO)
    print(np.mean(meta_mc_storage, axis=0))
    print((time.time() - Start) / 10)



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

os.chdir("../..")
t_start_first =time.time()
for i in range(0,300):
    flux_model = pickle.load(open("Data/Intermediate/Experimentally_Aligned_Model/A549_recon2_2_GLUN_and_SERtm_added.mdl", "rb"))
    t_start =time.time()
    efl = flux_model.find_essential_fluxes_fast_nps(seed=i)
    print(time.time() - t_start)
    print((time.time() - t_start_first)/(i+1)*(250-i-1)/60)
    works,err_type = flux_model.test_essential_list_nps(efl)
    if not works:
        print("DOES NOT WORK")
        print(i, err_type)
        with open("Data/Intermediate/Essential_Flux_Lists/exhauset_EF_A549_recon2_2_GLUN_and_SERtm_added_Error_List_test.txt", "a+") as EFO:
            EFO.write(str(efl)+f"\n{i, err_type}, seed = {i}, nps\n")
    if works:
        with open("Data/Intermediate/Essential_Flux_Lists/exhaust_EF_A549_recon2_2_GLUN_and_SERtm_added.txt", "a+") as EFO:
            EFO.write(str(efl)+"\n")




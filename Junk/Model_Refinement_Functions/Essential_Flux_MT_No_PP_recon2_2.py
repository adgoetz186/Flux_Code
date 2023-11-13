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
for i in range(0,250):
    flux_model = pickle.load(open("Data/Intermediate/Experimentally_Aligned_Model/A549_recon2_2_GLUN_and_SERtm_added.mdl", "rb"))
    t_start =time.time()
    seed_val = time.time()*i
    efl = flux_model.find_essential_fluxes_fast(seed=seed_val)
    print(time.time() - t_start)
    works,err_type = flux_model.test_essential_list(efl)
    if not works:
        print("DOES NOT WORK")
        print(i, err_type)
        with open("Data/Intermediate/Essential_Flux_Lists/exhauset_EF_A549_recon2_2_GLUN_and_SERtm_added_Error_List.txt", "a+") as EFO:
            EFO.write(str(efl)+f"\n{i, err_type}, seed = {seed_val}, ps\n")
    if works:
        with open("Data/Intermediate/Essential_Flux_Lists/exhaust_EF_A549_recon2_2_GLUN_and_SERtm_added.txt", "a+") as EFO:
            EFO.write(str(efl)+"\n")




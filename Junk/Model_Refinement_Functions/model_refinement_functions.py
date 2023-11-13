import numpy as np
import os
import pickle
import copy
from functools import partial
import multiprocessing as mp
import time

def single_cache_run(saved_model,repeat,no,**kwargs):
    initial_flux_counts = dict(zip(saved_model.reaction_names, [0 for i in range(len(saved_model.reaction_names))]))
    zero_list = np.zeros_like(saved_model.lb)
    for i in range(repeat):
        flux_model = copy.deepcopy(saved_model)
        zero_list+=np.array(list(flux_model.find_essential_fluxes(dict_count=initial_flux_counts,seed=no,**kwargs).values()))
    return zero_list
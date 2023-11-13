import os
import time
import numpy as np
import pandas as pd
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model


model_file_location = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Models/pkl_models/recon2_2_t_cell/exp_aligned/B6-1_exp_aligned"

output_model_file_location = "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Minimal_Model_Data/recon2_2_t_cell/"
filename_base = "B6-1"
recon2_2_exp_align = Flux_Balance_Model()
recon2_2_exp_align.load_pkl_model(model_file_location)
print(recon2_2_exp_align.reaction_info(""))
input()
print("start")
recon2_2_exp_align.generate_essential_flux_dataframe(2,output_model_file_location,filename_base,seed_list=[0,1],print_progress = True)
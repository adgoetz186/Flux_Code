import numpy as np
import os
import pickle
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Save_And_Load import save_object
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Data_Analysis
#Sets cwd as the main file directory
os.chdir("../../..")
flux_model = pickle.load(open("Data/Intermediate/Warmup_Points/High_Oxygen_Ecoli_Positive_Flux_10_percent_wup.mdl", "rb"))
print(flux_model.warmup_points)
biomass_obj = flux_model.objective_function_generator(["Biomass_Ecoli_core_w_GAM"],[1])
max_val = flux_model.find_objective_value(biomass_obj)

flux_obj = flux_model.objective_function_generator(np.arange(np.size(flux_model.ub)),-1*np.ones_like(flux_model.ub))
min_flux = flux_model.find_objective_value(flux_obj)

flux_model.ACHRSampler(biomass_obj,flux_model.warmup_points,"test",100,100,True,objective_constraint=[">", max_val*.90,f"np.array({flux_obj.tolist()}) @ react_flux >= {min_flux*1.1}"])
save_object(flux_model,"Data/Output/Sampled_Points/High_Oxygen_Ecoli_Positive_Flux_10_percent_sfs.mdl")
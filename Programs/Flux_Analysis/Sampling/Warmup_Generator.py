import numpy as np
import os
import pickle
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Save_And_Load import save_object
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Data_Analysis

#Sets cwd as the main file directory
os.chdir("../../..")

#Takes in a file containing a instance of the flux balance model and generates warmuppoints
#Resaves the file

flux_model = pickle.load(open("Data/Intermediate/Raw_MDL_Models/High_Oxygen_Ecoli_Positive_Flux_10_percent.mdl", "rb"))
biomass_obj = flux_model.objective_function_generator(["Biomass_Ecoli_core_w_GAM"],[1])
max_val = flux_model.find_objective_value(biomass_obj)

flux_obj = flux_model.objective_function_generator(np.arange(np.size(flux_model.ub)),-1*np.ones_like(flux_model.ub))
min_flux = flux_model.find_objective_value(flux_obj)

flux_model.generate_warmup_spanning(1e-9, 1e4,True,extra_constraints=[f"np.array({biomass_obj.tolist()}) @ react_flux >= {max_val*.9}",f"np.array({flux_obj.tolist()}) @ react_flux >= {min_flux*1.1}"],save=True)
print(flux_model.warmup_points)
save_object(flux_model,"Data/Intermediate/Warmup_Points/High_Oxygen_Ecoli_Positive_Flux_10_percent_wup.mdl")
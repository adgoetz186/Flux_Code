import pickle
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Save_And_Load import save_object
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import os
import numpy as np
os.chdir("../..")
flux_model = pickle.load(open("Data/Intermediate/Raw_MDL_Models/ecoli_core_model.mdl", "rb"))


exc = flux_model.get_exchanged_reaction_info()
print(exc)
print(flux_model.species_comp[flux_model.species_names.index("ac[e]")])
for i in exc.keys():
    exchange_composition = flux_model.species_comp[flux_model.species_names.index(i)]
    if "C" in exchange_composition and i != "glc-D[e]":
        flux_model.update_reaction_bounds(exc[i],0,"keep")
    elif i == "glc-D[e]":
        flux_model.update_reaction_bounds(exc[i], -10, -9)
flux_model.update_reaction_bounds("EX_o2(e)",-100,"keep")
biomass_obj = flux_model.objective_function_generator(["Biomass_Ecoli_core_w_GAM"],[1])
max_val = flux_model.find_objective_value(biomass_obj)
flux_model.reaction_info("Biomass_Ecoli_core_w_GAM")
flux_model.initalize_model(shrink_S=0.001,extra_constraints=[f"np.array({biomass_obj.tolist()}) @ react_flux >= {max_val*.9}"])
flux_model.reaction_info("Biomass_Ecoli_core_w_GAM")
flux_model.convert_model_to_positive_flux()
biomass_obj = flux_model.objective_function_generator(["Biomass_Ecoli_core_w_GAM"],[1])
flux_obj = flux_model.objective_function_generator(np.arange(np.size(flux_model.ub)),-1*np.ones_like(flux_model.ub))
min_flux = flux_model.find_objective_value(flux_obj)
print(max_val)
print(12)
print(flux_model.fva()[12])
input()
flux_model.initalize_model(shrink_S=0.001,extra_constraints=[f"np.array({biomass_obj.tolist()}) @ react_flux >= {max_val*.9}",f"np.array({flux_obj.tolist()}) @ react_flux >= {min_flux*1.1}"])
biomass_obj2 = flux_model.objective_function_generator(["Biomass_Ecoli_core_w_GAM"],[1])
#flux_obj = flux_model.objective_function_generator(np.arange(np.size(flux_model.ub)),-1*np.ones_like(flux_model.ub))
#min_flux2 = flux_model.find_objective_value(flux_obj)
#max_val2 = flux_model.find_objective_value(biomass_obj)
print(flux_model.fva()[12])
input()
warmup = flux_model.generate_warmup_spanning(0,1000,False,extra_constraints=[f"np.array({biomass_obj.tolist()}) @ react_flux >= {max_val*.9}",f"np.array({flux_obj.tolist()}) @ react_flux >= {min_flux*1.1}"])
print(warmup)
sample_points = flux_model.ACHRSampler([0],warmup,"test",100,1000,False)
print(sample_points)
np.savetxt("Programs/Flux_Analysis/Model_Loading/Glucose_Sampler/high_oxy_test",sample_points,delimiter=",",header=str(flux_model.reaction_names))

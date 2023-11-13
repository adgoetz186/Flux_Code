import pickle
import gurobipy as gp
import copy
import scipy.optimize as so
from gurobipy import GRB
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Save_And_Load import save_object
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
from Flux_Code.Programs.python_functions import File_handling
import numpy as np
File_handling.up_n_files(4)


flux_modelb = pickle.load(open("Data/Intermediate/Raw_MDL_Models/RECON1b.mdl", "rb"))
flux_modelb.reaction_info("RDH3a")

flux_modelb.pinch_restricted_exchange_reactions("Data/Input/Experimental_Data/RECON1/Allowed_Exchange_Reactions")

for i in range(len(flux_modelb.reaction_names)):
    if "maintenance_ATP" in flux_modelb.reaction_names[i]:
        flux_modelb.update_reaction_bounds(flux_modelb.reaction_names[i],2,2)
    if "biomass" in flux_modelb.reaction_names[i]:
        flux_modelb.update_reaction_bounds(flux_modelb.reaction_names[i],0.03025,0.03025)

#In the original setup:
#CITtbm is the 692nd reaction
#CK is the 693rd reaction
#CLFORtex is the 695th reaction
flux_modelb.update_reaction_bounds("CITtbm","keep",0)
flux_modelb.update_reaction_bounds("CK","keep",0)
flux_modelb.update_reaction_bounds("CLFORtex",0,0)
#flux_modelb.pinch_reactions(10**-100)
#flux_modelb.pinch_reactions(10**-100)
#flux_modelb.update_reaction_bounds("EX_gly(e)",-100,100)

"Data/Input/Experimental_Data/Measured_Rates.txt"
flux_modelb.fit_to_experimental_data("Data/Input/Experimental_Data/RECON1/Measured_Rates.txt","Data/Intermediate/Optimal_Size_Factor/A549.txt",10.0**np.linspace(-2.75,-1,500),1)
#flux_modelb.reaction_info('EX_asp_L(e)')
#flux_modelb.reaction_info('EX_met_L(e)')
#flux_modelb.reaction_info('EX_thr_L(e)')
#flux_modelb.pinch_reactions(10**-100)
save_object(flux_modelb, "Data/Intermediate/Experimentally_Aligned_Model/A549_RECON1.mdl")

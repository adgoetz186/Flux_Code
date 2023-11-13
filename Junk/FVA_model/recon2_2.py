import os
import pickle
import gurobipy as gp
import copy
import scipy.optimize as so
from gurobipy import GRB
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import numpy as np
import cplex
os.chdir("../..")

flux_modelb = pickle.load(open("Data/Intermediate/Raw_MDL_Models/recon2_2b.mdl", "rb"))

for i in range(len(flux_modelb.reaction_names)):
    if "DM_atp_c_" in flux_modelb.reaction_names[i]:
        flux_modelb.update_reaction_bounds(flux_modelb.reaction_names[i],2,2)
        #flux_modelb.update_reaction_bounds(flux_modelb.reaction_names[i], 2, 2)
    if "biomass_reaction" in flux_modelb.reaction_names[i]:
        flux_modelb.update_reaction_bounds(flux_modelb.reaction_names[i],0.03025, 0.03025)
        #flux_modelb.update_reaction_bounds(flux_modelb.reaction_names[i], 0.03025, 0.03025 )
flux_modelb.update_reaction_bounds("CITtbm","keep",0)
flux_modelb.update_reaction_bounds("CK","keep",0)
flux_modelb.update_reaction_bounds("CLFORtex",0,0)
for i in range(len(flux_modelb.lb)):
    if flux_modelb.ub[i] == np.inf:
        flux_modelb.update_reaction_bounds(i,"keep",1000)
    if flux_modelb.lb[i] == -np.inf:
        flux_modelb.update_reaction_bounds(i,-1000,"keep")
flux_modelb.pinch_reactions(10**-100)
print(flux_modelb.test_feasibility())
flux_modelb.initalize_model(shrink_S = 10**-100)

save_object(flux_modelb, "Data/Intermediate/Raw_MDL_Models/recon2_2b_FVA.mdl")
import os
import pickle
import gurobipy as gp
import copy
import copy as cp
import scipy.optimize as so
from gurobipy import GRB
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Save_And_Load import save_object
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import numpy as np
from matplotlib import pyplot as plt
os.chdir("../../..")


flux_modelb = pickle.load(open("Data/Intermediate/Raw_MDL_Models/recon2_2_GLUN_and_SERtm_added.mdl", "rb"))
print(flux_modelb.test_feasibility())
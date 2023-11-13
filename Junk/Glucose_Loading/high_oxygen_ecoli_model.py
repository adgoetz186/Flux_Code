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
os.chdir("../..")


flux_model = pickle.load(open("Data/Output/Sampled_Points/High_Oxygen_Ecoli_Positive_Flux_10_percent_sfs.mdl", "rb"))
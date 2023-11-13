import os
import pickle
import gurobipy as gp
import copy
import scipy.sparse as sp
import copy as cp
import scipy.optimize as so
from gurobipy import GRB
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Save_And_Load import save_object
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions import Flux_Model_Class

import numpy as np
from matplotlib import pyplot as plt

with open("high_oxy_test", "r") as file:
	names = file.readline()
names = names.replace(" ","").replace("#","").replace("'","")[1:-2].split(",")
print(names)
points = np.loadtxt("high_oxy_test",delimiter=",")
print(points[:,12])
plt.hist(points[:,12])
plt.show()
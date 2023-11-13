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

os.chdir("../../../..")
S_val = np.zeros((3,3))
S_val[0,0] = 1
S_val[0,1] = 1
S_val[0,2] = -2
S_val = sp.csr_matrix(S_val)
flux_modelb = Flux_Model_Class.Flux_Balance_Model("Cube",S_val, np.zeros(3),np.ones(3),np.zeros(3))
print(flux_modelb.S)
size = 0
rank = 0
for i in range(100):
	warmup = flux_modelb.generate_warmup_spanning(0,1000,False)
	size+=np.shape(warmup)[0]
	rank+=np.linalg.matrix_rank(warmup)
	center = np.average(warmup, axis=0)
print(size/100)
print(warmup)
print(rank/100)
input()
print(os.getcwd())
sample_points = flux_modelb.HRSampler(np.zeros(3),warmup,"test",100000,2,False, display = ["save","Programs/Flux_Analysis/Model_Loading/Simple_Geometric_Model/cube","cube"])

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')

sequence_containing_x_vals = list(sample_points[:,0])
sequence_containing_y_vals = list(sample_points[:,1])
sequence_containing_z_vals = list(sample_points[:,2])
ax.scatter(sequence_containing_x_vals, sequence_containing_y_vals, sequence_containing_z_vals)
ax.scatter(warmup[:,0], warmup[:,1], warmup[:,2],s = 50)
ax.set_xlabel('$X$', fontsize=20)
ax.set_ylabel('$Y$', fontsize=20)
ax.set_zlabel('$Z$', fontsize=20)
plt.show()

fig, (ax1,ax2,ax3) = plt.subplots(nrows=3)
ax1.hist(list(sample_points[:,0]),density = "pdf",bins=50)
ax1.set_xlabel("X value")
ax2.hist(list(sample_points[:,1]),density = "pdf",bins=50)
ax2.set_xlabel("Y value")
ax2.set_ylabel("Density of Sampling")
ax3.hist(list(sample_points[:,2]),density = "pdf",bins=50)
ax3.set_xlabel("Z value")

plt.show()
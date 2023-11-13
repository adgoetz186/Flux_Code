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
size = 3
S_val = np.zeros((size,size))
S_val = np.round(np.random.random((size,size)),0)
S_val[0] = np.zeros((size))
print(S_val)
S_val = sp.csr_matrix(S_val)
flux_modelb = Flux_Model_Class.Flux_Balance_Model("Cube",S_val, -1*np.arange(size)-1.0,np.arange(size)+1.0,np.zeros(size))
print(flux_modelb.S)
print(1)
flux_modelb.save_model("test2",dir = "Programs/Flux_Analysis/Model_Loading/Simple_Geometric_Model/")

flux_model_c = Flux_Model_Class.Flux_Balance_Model()
flux_model_c.load_model("test2",dir = "Programs/Flux_Analysis/Model_Loading/Simple_Geometric_Model/")

warmup = flux_modelb.generate_warmup_cp(1000)


sample_points = flux_modelb.HRSampler(warmup,10000,200)
plt.hist(sample_points[:,0])
plt.show()
plt.hist(sample_points[:,1])
plt.show()
plt.hist(sample_points[:,2])
plt.show()
distance = 0
for i in range(len(sample_points)-1):
	distance += np.sqrt(np.sum((sample_points[i]-sample_points[i+1])**2))
print(distance/(len(sample_points)-1))
input()

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
import importlib.util
import sys
import os
import pickle
import gurobipy as gp
import copy
import copy as cp
import scipy.optimize as so
from gurobipy import GRB

print(importlib.util.find_spec("Save_And_Load"))
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import numpy as np
import time
from matplotlib import pyplot as plt



flux_model = pickle.load(open("Data/Intermediate/Essential_Flux_Model/A549_recon2_2_GLUN_and_SERtm_added.mdl", "rb"))

warmup = flux_model.generate_warmup_gb(1000)
#print(warmup)
#start = time.time()
point_skip = [1,2,4,5,8,16]
for i in point_skip:
	sample_points = flux_model.HRSampler([0],warmup,"test",1000,i,False)
	np.savetxt("Programs/Flux_Analysis/Model_Loading/Recon_Sampler/"+str(i)+"_sec_2",sample_points,delimiter=",",header=str(flux_model.reaction_names))
#print(time.time()-start,1)
#sample_points = flux_model.ACHRSampler([0],warmup,"test",10000,2,False)
#np.savetxt("Programs/Flux_Analysis/Model_Loading/Recon_Sampler/2_sec",sample_points,delimiter=",",header=str(flux_model.reaction_names))
#print(time.time()-start,2)
#sample_points = flux_model.ACHRSampler([0],warmup,"test",10000,4,False)
#np.savetxt("Programs/Flux_Analysis/Model_Loading/Recon_Sampler/4_sec",sample_points,delimiter=",",header=str(flux_model.reaction_names))
#print(time.time()-start,4)
#sample_points = flux_model.ACHRSampler([0],warmup,"test",10000,8,False)
#np.savetxt("Programs/Flux_Analysis/Model_Loading/Recon_Sampler/8_sec",sample_points,delimiter=",",header=str(flux_model.reaction_names))
#print(time.time()-start,8)
#sample_points = flux_model.ACHRSampler([0],warmup,"test",10000,16,False)
#np.savetxt("Programs/Flux_Analysis/Model_Loading/Recon_Sampler/16_sec",sample_points,delimiter=",",header=str(flux_model.reaction_names))
#print(time.time()-start,16)
#sample_points = flux_model.ACHRSampler([0],warmup,"test",10000,32,False)
#np.savetxt("Programs/Flux_Analysis/Model_Loading/Recon_Sampler/32_sec",sample_points,delimiter=",",header=str(flux_model.reaction_names))
#print(time.time()-start,32)
#sample_points = flux_model.ACHRSampler([0],warmup,"test",10000,64,False)
#np.savetxt("Programs/Flux_Analysis/Model_Loading/Recon_Sampler/64_sec",sample_points,delimiter=",",header=str(flux_model.reaction_names))
#print(time.time()-start,64)
#sample_points = flux_model.ACHRSampler([0],warmup,"test",10000,128,False)
#np.savetxt("Programs/Flux_Analysis/Model_Loading/Recon_Sampler/128_sec",sample_points,delimiter=",",header=str(flux_model.reaction_names))
#print(time.time()-start,128)
#sample_points = flux_model.ACHRSampler([0],warmup,"test",10000,256,False)
#np.savetxt("Programs/Flux_Analysis/Model_Loading/Recon_Sampler/256_sec",sample_points,delimiter=",",header=str(flux_model.reaction_names))
#print(time.time()-start,256)
#sample_points = flux_model.ACHRSampler([0],warmup,"test",10000,512,False)
#np.savetxt("Programs/Flux_Analysis/Model_Loading/Recon_Sampler/512_sec",sample_points,delimiter=",",header=str(flux_model.reaction_names))
#print(time.time()-start,512)
#sample_points = flux_model.ACHRSampler([0],warmup,"test",10000,1042,False)
#np.savetxt("Programs/Flux_Analysis/Model_Loading/Recon_Sampler/1024_sec",sample_points,delimiter=",",header=str(flux_model.reaction_names))
#print(time.time()-start,1042)
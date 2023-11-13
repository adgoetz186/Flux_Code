import numpy as np
import copy as cp
import os
import pickle

os.chdir("../..")
with open("Data/Intermediate/Essential_Flux_Lists/EF_A549_recon2_2d.txt","r") as eff:
	flux_lists_string = eff.readlines()
flux_lists = [eval(i) for i in flux_lists_string]
flux_model = pickle.load(open("Data/Intermediate/Experimentally_Aligned_Model/A549_recon2_2d.mdl", "rb"))
count = 0
for i in flux_lists:
	print((count)/len(flux_lists))
	if not flux_model.test_essential_list(i):
		print("PAIN")
		input()
	count+=1
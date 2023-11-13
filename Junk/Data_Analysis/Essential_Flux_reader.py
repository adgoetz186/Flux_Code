import json
import os
import numpy as np

os.chdir("../..")
with open("Data/Intermediate/Essential_Flux_Values/EF_A549_RECON1_pinched.json") as EF:
	flux_dict = json.load(EF)
print(flux_dict)
efv = np.array([flux_dict[i] for i in flux_dict.keys()])

print(np.mean(efv))
nonzero_count = 0
for i in efv:
	if i>0:
		nonzero_count+=1
print(nonzero_count/len(efv))
print(len(efv))

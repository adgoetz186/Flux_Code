import os
import sys
import pandas as pd
import numpy as np
import scipy.io as sio
from pathlib import Path
import matplotlib.pyplot as plt
import time
import scipy.stats as st
import gurobipy as gp

# _____ Setting the CWD to be Flux_Code BEGIN _____
# Flux_Code path here
if not Path(os.getcwd()).parts[-1] == "Flux_Code":
	path_to_FC = ""
	if path_to_FC == "":
		try:
			# Obtains the location of the Cell_signaling_information folder if it is in the cwd parents
			path_to_FC = Path.cwd().parents[
				[Path.cwd().parents[i].parts[-1] for i in range(len(Path.cwd().parts) - 1)].index(
					"Flux_Code")]
		except ValueError:
			print("Flux_Code not found in cwd parents, trying sys.path")
			try:
				# Obtains the location of the Cell_signaling_information folder if it is in sys.path
				path_to_CSI = Path(sys.path[[Path(i).parts[-1] for i in sys.path].index("Flux_Code")])
			except ValueError:
				print("Flux_Code not found in sys.path "
					  "consult 'Errors with setting working directory' in README")
	else:
		path_to_CSI = Path(path_to_FC)
	os.chdir(path_to_FC)

from Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
# _____ Setting the CWD to be Flux_Code END _____


def consecutive_points_to_distances(point_array,max_range,points_to_average,toss = 0,toss_dist = 0):
	dist_values = []
	for dist in range(1, max_range+1):
		average_dist = 0
		for point_in in range(points_to_average):
			# avoid the total distance point
			average_dist += np.linalg.norm(point_array[dist * point_in+toss*toss_dist] - point_array[dist * (point_in + 1)+toss*toss_dist])
			print(dist * (point_in+toss),dist * (point_in + 1+toss),toss)
		average_dist /= points_to_average
		dist_values.append(average_dist)
	return dist_values

def folder_to_dist_plot(path_to_folder,rxns,max_range,points_to_average,drop_points=0):
	rxns.append("Total_Flux")
	all_points = np.empty(0)
	all_distances_zero_removed = np.empty((len(os.listdir(path_to_folder)),max_range))
	all_distances_w_zero = np.empty((len(os.listdir(path_to_folder)),max_range ))


	file_ind = 0
	for file in os.listdir(path_to_folder):
		if np.size(all_points) == 0:
			points = np.load(path_to_folder/Path(file))
			all_points = points
			print(np.shape(points))
			input()
			all_distances_zero_removed[file_ind] = consecutive_points_to_distances(points[:,:-1],max_range,points_to_average,toss=drop_points,toss_dist = 100)
		else:
			points = np.load(path_to_folder/Path(file))
			all_points = np.vstack((all_points,points))
			all_distances_zero_removed[file_ind] = consecutive_points_to_distances(points[:, :-1], max_range, points_to_average, toss=drop_points,toss_dist = 100)
		file_ind +=1
	print(all_distances_zero_removed)
	plt.plot(np.transpose(all_distances_zero_removed))
	plt.title(f"{path_to_folder.split('/')[-1]}\nTosses first {drop_points} points\n{int(len(os.listdir(path_to_folder))*points_to_average)} points to use")
	plt.xlabel("Tossed points between sampled points (x5000)")
	plt.ylabel("Distance between concurrent sample points")
	plt.tight_layout()
	plt.show()
	#plt.plot(np.transpose(all_distances_w_zero))
	#plt.show()




def folder_to_dtf(path_to_folder,rxns):
	all_points = np.empty(0)
	for file in os.listdir(path_to_folder):
		if np.size(all_points) == 0:
			points = np.load(path_to_folder/Path(file))
			all_points = points[1600::100,:-1]
		else:
			points = np.load(path_to_folder/Path(file))[1600::100,:-1]
			all_points = np.vstack((all_points,points))
	return pd.DataFrame(data=all_points, columns=rxns)
recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_fast_key_json_model("Data/Models/json_models/fast_key_format/recon_1_t_cells_12_11_23/HR_Ready_Day_2_single_size/B6-1.json")
lb, ub, S, b, rxn_list, met_list = recon_flux_model.dicts_to_mats()
for rxn in rxn_list:
	if "got" in rxn.lower():
		print(rxn)
print("done")

print(np.arange(2601)[1100::100])

#folder_to_dist_plot("Data/HR/HR_fin/B6_4_gene_therm",rxn_list)
dataframe_list = ["B6_1","B6_2","B6_3","B6_4","TC_5","TC_6","TC_7","TC_8","B6_1_gene_therm","B6_2_gene_therm","B6_3_gene_therm","B6_4_gene_therm","TC_5_gene_therm","TC_6_gene_therm","TC_7_gene_therm","TC_8_gene_therm"]
B6_dtf_list = []
TC_dtf_list = []
for condition_name in dataframe_list:
	print(condition_name)
	folder_to_dist_plot(f"Data/HR/HR_fin/{condition_name}",rxn_list,10,points_to_average=1000)
	#if "B6" in condition_name:
	#	points_dtf = folder_to_dtf(f"Data/HR/HR_fin/{condition_name}", rxn_list)
	#	B6_dtf_list.append(points_dtf.mean())
	#if "TC" in condition_name:
	#	points_dtf = folder_to_dtf(f"Data/HR/HR_fin/{condition_name}",rxn_list)
	#	TC_dtf_list.append(points_dtf.mean())
B6_array = pd.concat(B6_dtf_list,axis=1).transpose()
TC_array = pd.concat(TC_dtf_list,axis=1).transpose()
for rxn in B6_array.columns:
	if "_inverted" in rxn:
		if rxn.replace("_inverted","") in B6_array.columns:
			B6_array[rxn.replace("_inverted","")] = B6_array[rxn.replace("_inverted","")]-B6_array[rxn]
			B6_array.drop(columns=[rxn],inplace=True)
		else:
			B6_array[rxn] = -1*B6_array[rxn]
			B6_array.rename(columns = {rxn:rxn.replace("_inverted","")},inplace=True)

for rxn in TC_array.columns:
	if "_inverted" in rxn:
		if rxn.replace("_inverted","") in TC_array.columns:
			TC_array[rxn.replace("_inverted","")] = TC_array[rxn.replace("_inverted","")]-TC_array[rxn]
			TC_array.drop(columns=[rxn],inplace=True)
		else:
			TC_array[rxn] = -1*TC_array[rxn]
			TC_array.rename(columns = {rxn:rxn.replace("_inverted","")},inplace=True)
print(B6_array.columns)
print(TC_array)

rxns_of_interest = ["biomass_reaction","EX_glc(e)","EX_lac_L(e)","PYRt2m","CITtam","G6PDH2r","CITL","ACITL","ACONTm","PDHm","CSm","EX_gln_L(e)","EX_glu_L(e)","ASPTA"]
for rxn in rxns_of_interest:
	print(rxn)
	print("B6",B6_array[rxn].mean(),B6_array[rxn].std())
	print("TC",TC_array[rxn].mean(),TC_array[rxn].std())
	#print("B6",B6_array[rxn])
	#print("TC", TC_array[rxn])
	print(st.ttest_ind(B6_array[rxn].to_numpy(),TC_array[rxn].to_numpy(),equal_var = False))
	#print(B6_array[rxn])
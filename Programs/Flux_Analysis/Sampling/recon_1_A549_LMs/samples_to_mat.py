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

def folder_to_dist_plot(path_to_folder,rxns):
	rxns.append("Total_Flux")
	all_points = np.empty(0)
	max_range = 100
	all_distances_zero_removed = np.empty((len(os.listdir(path_to_folder)),max_range))
	all_distances_w_zero = np.empty((len(os.listdir(path_to_folder)),max_range ))
	points_to_average = 10
	points_to_drop = 16
	file_ind = 0
	for file in os.listdir(path_to_folder):
		if np.size(all_points) == 0:
			points = np.load(path_to_folder/Path(file))
			all_points = points
			print(np.shape(points))
			all_distances_zero_removed[file_ind] = consecutive_points_to_distances(points[:,:-1],100,points_to_average,toss=points_to_drop,toss_dist = 100)
		else:
			points = np.load(path_to_folder/Path(file))
			all_points = np.vstack((all_points,points))
			all_distances_zero_removed[file_ind] = consecutive_points_to_distances(points[:, :-1], 100, points_to_average, toss=points_to_drop,toss_dist = 100)
		file_ind +=1
	print(all_distances_zero_removed)
	plt.plot(np.transpose(all_distances_zero_removed))
	plt.title(f"Tosses first {points_to_drop} points\n{int(len(os.listdir(path_to_folder))*points_to_average)} points to use")
	plt.xlabel("Tossed points between sampled points (x500)")
	plt.ylabel("Distance between concurrent sample points")
	plt.show()
	#plt.plot(np.transpose(all_distances_w_zero))
	#plt.show()




def folder_to_dtf(path_to_folder,rxns,start,skip,stop=-1):
	all_points = np.empty(0)
	for file in os.listdir(path_to_folder):
		if np.size(all_points) == 0:
			points = np.load(path_to_folder/Path(file))
			all_points = points[start:stop:skip,:-1]
		else:
			points = np.load(path_to_folder/Path(file))[start:stop:skip,:-1]
			all_points = np.vstack((all_points,points))
	return pd.DataFrame(data=all_points, columns=rxns)
recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_fast_key_json_model("Data/Models/json_models/fast_key_format/recon_1_t_cells_12_11_23/HR_Ready_Day_2_single_size/B6-1.json")
lb, ub, S, b, rxn_list, met_list = recon_flux_model.dicts_to_mats()
for rxn in rxn_list:
	if "got" in rxn.lower():
		print(rxn)
print(rxn_list)
print(len(rxn_list))
print("done")
input()
print(np.arange(2601)[1100::100])

#folder_to_dist_plot("Data/HR/HR_fin/B6_4_gene_therm",rxn_list)
dataframe_list = ["B6_1_gene_therm","B6_2_gene_therm","B6_3_gene_therm","B6_4_gene_therm","TC_5_gene_therm","TC_6_gene_therm","TC_7_gene_therm","TC_8_gene_therm","B6_1","B6_2","B6_3","B6_4","TC_5","TC_6","TC_7","TC_8"]
dtf_list = []
TC_dtf_list = []
for condition_name in dataframe_list:
	if "therm" in condition_name:
		points_dtf = folder_to_dtf(f"Data/HR/HR_fin/{condition_name}", rxn_list,1700,100,stop = 2601)
	else:
		points_dtf = folder_to_dtf(f"Data/HR/HR_fin/{condition_name}", rxn_list, 0, 10)
	for rxn in points_dtf.columns:
		if "_inverted" in rxn:
			if rxn.replace("_inverted", "") in points_dtf.columns:
				points_dtf[rxn.replace("_inverted", "")] = points_dtf[rxn.replace("_inverted", "")] - points_dtf[rxn]
				points_dtf.drop(columns=[rxn], inplace=True)
			else:
				points_dtf[rxn] = -1 * points_dtf[rxn]
				points_dtf.rename(columns={rxn: rxn.replace("_inverted", "")}, inplace=True)
	rxn_names = points_dtf.columns
	dtf_list.append(points_dtf.to_numpy())
	print(np.shape(points_dtf.to_numpy()))


point_dict = dict(zip(dataframe_list,dtf_list))
print(point_dict)
point_dict["rxn_names"] = rxn_names
print(len(point_dict["rxn_names"]))
print(point_dict)
for key in point_dict.keys():
	if key != "rxn_names":
		print(key,np.shape(point_dict[key]))
sio.savemat("Data/HR/HR_fin/mat_file/sample_points.mat",point_dict)

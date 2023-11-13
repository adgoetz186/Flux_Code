import os
import pandas as pd
import numpy as np
from pathlib import Path
import time
import scipy.stats as st
import matplotlib.pyplot as plt
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model

# _____ Setting the CWD to be Flux_Code BEGIN _____
# Flux_Code path here
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
# _____ Setting the CWD to be Flux_Code END _____

HR_model_location = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/HR_Ready_Day_2_single_size_B6_1_A")

HR_Points_location_B6_1_path = Path("Data/HR/recon_1b_t_cells_gene_Day_2_single_B6_1_gene_therm")
HR_Points_location_B6_1_all = np.array([])
distaces = np.array([])
start = 5
max_distance = 100
values = 5

point_sample_seperation = np.arange(1, max_distance + 1)
for file in os.listdir(HR_Points_location_B6_1_path):
	
	new_points = np.load(HR_Points_location_B6_1_path / file)
	points_for_distances = new_points[start:start + max_distance * values + 1:max_distance]
	if np.size(HR_Points_location_B6_1_all) == 0:
		HR_Points_location_B6_1_all = points_for_distances
	else:
		HR_Points_location_B6_1_all = np.vstack((HR_Points_location_B6_1_all, points_for_distances))
len_total = np.shape(HR_Points_location_B6_1_all)
print(len_total[0] // 2)
HR_Points_location_B6_1_A = HR_Points_location_B6_1_all[len_total[0] // 2:]
HR_Points_location_B6_1_B = HR_Points_location_B6_1_all[:len_total[0] // 2]
print(np.shape(HR_Points_location_B6_1_A))
print(np.shape(HR_Points_location_B6_1_B))
input()

flux_model_dict = {}
for filename in os.listdir(HR_model_location):
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_fast_key_json_model(HR_model_location / filename)
	flux_model_dict[filename.split(".")[0]] = recon_flux_model

lb, ub, S, b, rxn_list, met_list = flux_model_dict['B6-1'].dicts_to_mats()

cd_list = []
for i in range(np.shape(HR_Points_location_B6_1_B)[1]):
	variance_B6_1_A = np.var(HR_Points_location_B6_1_A[:, i])
	mean_B6_1_A = np.mean(HR_Points_location_B6_1_A[:, i])
	
	variance_B6_1_B = np.var(HR_Points_location_B6_1_B[:, i])
	mean_B6_1_B = np.mean(HR_Points_location_B6_1_B[:, i])
	# print(st.ttest_ind(points_A[:, i], points_B[:, i]))
	cd_list.append((mean_B6_1_A - mean_B6_1_B) / np.sqrt((variance_B6_1_A + variance_B6_1_B)))
	print(mean_B6_1_A, mean_B6_1_B, np.sqrt((variance_B6_1_A + variance_B6_1_B)))
	print(i)
	print(variance_B6_1_A / mean_B6_1_A)
# plt.title(rxn_list[i])
# plt.hist(points_B[:,i],bins = np.linspace(np.min(np.vstack((points_B[:,i],points_A[:,i]))),np.max(np.vstack((points_B[:,i],points_A[:,i]))),50),alpha = 0.5)
# plt.hist(points_A[:, i], bins = np.linspace(np.min(np.vstack((points_B[:,i],points_A[:,i]))),np.max(np.vstack((points_B[:,i],points_A[:,i]))),50), alpha=0.5)
# plt.show()
print(cd_list)
cd_list = np.array(cd_list)
plt.scatter(np.log10(np.abs(cd_list)),
            np.log10(np.mean(HR_Points_location_B6_1_A, axis=0) + np.mean(HR_Points_location_B6_1_B, axis=0)))
plt.show()
cd_list_ind = np.flip(np.argsort(np.abs(cd_list)))
print(cd_list)
for i in cd_list_ind:
	print(rxn_list[i])
	print(cd_list[i], i)
	if not np.isnan(cd_list[i]):
		plt.title(rxn_list[i])
		print()
		print(HR_Points_location_B6_1_B[:, i])
		print(HR_Points_location_B6_1_A[:, i])
		plt.hist(HR_Points_location_B6_1_B[:, i],
		         bins=np.linspace(np.min(np.hstack((HR_Points_location_B6_1_B[:, i], HR_Points_location_B6_1_A[:, i]))),
		                          np.max(np.hstack((HR_Points_location_B6_1_B[:, i], HR_Points_location_B6_1_A[:, i]))),
		                          50), alpha=0.5,
		         label="B6_1 first half")
		plt.hist(HR_Points_location_B6_1_A[:, i],
		         bins=np.linspace(np.min(np.hstack((HR_Points_location_B6_1_B[:, i], HR_Points_location_B6_1_A[:, i]))),
		                          np.max(np.hstack((HR_Points_location_B6_1_B[:, i], HR_Points_location_B6_1_A[:, i]))),
		                          50), alpha=0.5,
		         label='B6_1 second half')
		plt.legend()
		plt.show()
input()



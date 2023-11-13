import os
import pandas as pd
import scipy as sc
import numpy as np
import scipy.io as sio
from pathlib import Path
import matplotlib.pyplot as plt
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model

# _____ Setting the CWD to be Flux_Code BEGIN _____
# Flux_Code path here
path_to_FC = ""
if path_to_FC == "":
	try:
		# Obtains the location of the Cell_signaling_information folder if it is in the cwd parents
		path_to_FC = Path.cwd().parents[[Path.cwd().parents[i].parts[-1] for i in range(len(Path.cwd().parts) - 1)].index(
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

HR_model_location = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/HR_Ready_Day_2_single_size_B6_1_A")

HR_Points_location_high_A = Path(
	"Data/HR/HR_Points_Distances/recon_1b_t_cells_gene_Day_2_single_B6_1_gene_med_lin_comb_A")

#HR_Points_location_high_A = Path("Data/HR/HR_Points_Distances/recon_1b_t_cells_gene_Day_2_single_B6_1_gene_high_lin_comb_A")

#HR_Points_location_high_A = Path("Data/HR/HR_Warmup/recon_1b_t_cells_gene_Day_2_single_size_B6_1_A")

HR_Points_location_high_B = Path(
	"Data/HR/HR_Points_Distances/recon_1b_t_cells_gene_Day_2_single_B6_1_gene_med_lin_comb_A")

for file in os.listdir(HR_Points_location_high_A):
	points_A = np.load(HR_Points_location_high_A / file)

for file in os.listdir(HR_Points_location_high_B):
	points_B = np.load(HR_Points_location_high_B / file)
points_A[np.abs(points_A)<1e-3] = 0

print(np.min(np.abs(points_A)[np.nonzero(np.abs(points_A))]))

flux_model_dict = {}
for filename in os.listdir(HR_model_location):
	print(filename)
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_fast_key_json_model(HR_model_location / filename)
	flux_model_dict[filename.split(".")[0]] = recon_flux_model
name_gene_dict = {"B6-1":"B61-D2-GEX-B5"}
lb, ub, S, b, rxn_list, met_list = flux_model_dict['B6-1'].dicts_to_mats()
print(rxn_list)
print(rxn_list.index("TKT2"))

rxn_list.append("Total_Flux")
print(len(rxn_list))

point_sample_seperation = np.arange(1, 501)
point_distances_orig = []
point_distances_A = []
point_distances_B = []
# distance_matrix = np.zeros((100-1,np.shape(points_A)[1]))
for flux_model_name in flux_model_dict.keys():
	gene_array = np.load(f"Data/scRNA/Bulk/{name_gene_dict[flux_model_name]}/matrix.npy")
	feature_array = np.load(f"Data/scRNA/Bulk/{name_gene_dict[flux_model_name]}/features.npy")
	gene_array = gene_array/np.average(gene_array)
	print(np.average(gene_array))

	flux_model = flux_model_dict[flux_model_name]
	RAS = np.array(flux_model.create_RAS_values(gene_array,feature_array,False))
	RPS = flux_model.RAS_to_RPS_mat(RAS)
	simp_neg_ind, comp_neg_ind, comp_perm, cutter, exchange_cutter = flux_model.precomp_positive_flux_point_to_bi_dir(
		rxn_list[:-1], prune_specific=["biomass_reaction"])
	test_bi_S = flux_model.positive_S_to_bi_dir(S, simp_neg_ind, cutter, exchange_cutter=exchange_cutter)
	
	NS = sc.linalg.null_space(test_bi_S)
	NS = sio.loadmat("Data/null_space_S/recon_1b_t_cells_single_size/Null_S_internal.mat")['Null_S_internal']
	print(NS)
	trans_NS = np.transpose(NS)

	count = 0
	feas_count = 0
	infeas_count = 0
	avg_val = 0
	feas_list = []
	for i in points_A:
		ent_flux_point = flux_model.positive_flux_point_to_bi_dir(i[:-1], simp_neg_ind, comp_neg_ind, comp_perm, cutter,
		                                                          exchange_cutter=exchange_cutter)
		val = flux_model.generate_therm_model_new(trans_NS, ent_flux_point)
		count+=1
		avg_val += val
		if count % 100 == 0:
			feas_list.append(avg_val/100)
			print(count,avg_val/100)
			avg_val = 0
		feas_list.append(val)
	plt.plot(np.arange(feas_list),feas_list)
	plt.show()
	error_list = []
	for point in points_A:
		val = flux_model.gene_point_penalty(point[:-1],RPS,1)
		error_list.append(val)
		int_point = self.positive_flux_point_to_bi_dir(curPoint[:-1], simp_neg_ind, comp_neg_ind, comp_perm, cutter, exchange_cutter=exchange_cutter)

	plt.plot(np.arange(len(error_list)),error_list)
	plt.show()
	print(RPS)
	input()
print(time_list)
plt.ylabel("Average Euclidean Distance")
plt.xlabel("Iteration difference")
plt.plot(point_sample_seperation * 500, point_distances_A)
plt.plot(point_sample_seperation * 500, point_distances_B)
plt.show()

# plt.hist(points_orig[:,-1],alpha = 0.5,bins = np.linspace(149,149.5,20))
# plt.hist(points_lc[:,-1],alpha = 0.5,bins = np.linspace(149,149.5,20))
# plt.show()


# plt.plot(np.arange(np.size(points_orig[:,-1])),points_orig[:,-1])
# plt.plot(np.arange(np.size(points_lc[:,-1])),points_lc[:,-1])
plt.plot((np.arange(np.size(points_A[:, -1]))), points_A[:, -1])
plt.plot((np.arange(np.size(points_B[:, -1]))), points_B[:, -1])
plt.show()

print(np.shape(points_A))
points_A = points_A[1000::1000]
points_B = points_B[1000::1000]
print(np.shape(points_A))
print(points_A[:, -1])

plt.hist(points_A[:, -1], bins=np.linspace(min(np.min(points_A[:, -1]), np.min(points_A[:, -1])),
                                           max(np.max(points_A[:, -1]), np.max(points_A[:, -1])), 50), alpha=0.5)
plt.hist(points_B[:, -1], bins=np.linspace(min(np.min(points_A[:, -1]), np.min(points_A[:, -1])),
                                           max(np.max(points_B[:, -1]), np.max(points_B[:, -1])), 50), alpha=0.5)
# plt.hist(points_ss[:,-1],bins = np.linspace(149,max(np.max(points_orig[:,-1]), np.max(points_lc[:,-1]),np.max(points_ss[:,-1])),20),alpha = 0.5)
plt.show()
cd_list = []
for i in range(np.shape(points_B)[1]):
	variance_A = np.var(points_A[:, i])
	mean_A = np.mean(points_A[:, i])
	
	variance_B = np.var(points_B[:, i])
	mean_B = np.mean(points_B[:, i])
	print(st.ttest_ind(points_A[:, i], points_B[:, i]))
	cd_list.append((mean_A - mean_B) / np.sqrt((variance_A + variance_B)))
	print(mean_A, mean_B, np.sqrt((variance_A + variance_B)))
	print(i)
	print(variance_A / mean_A)
# plt.title(rxn_list[i])
# plt.hist(points_B[:,i],bins = np.linspace(np.min(np.vstack((points_B[:,i],points_A[:,i]))),np.max(np.vstack((points_B[:,i],points_A[:,i]))),50),alpha = 0.5)
# plt.hist(points_A[:, i], bins = np.linspace(np.min(np.vstack((points_B[:,i],points_A[:,i]))),np.max(np.vstack((points_B[:,i],points_A[:,i]))),50), alpha=0.5)
# plt.show()
print(cd_list)
cd_list = np.array(cd_list)
plt.scatter(np.log10(np.abs(cd_list)), np.log10(np.mean(points_A, axis=0) + np.mean(points_B, axis=0)))
plt.show()
cd_list_ind = np.flip(np.argsort(np.abs(cd_list)))
print(cd_list)
for i in cd_list_ind:
	print(cd_list[i])
	if not np.isnan(cd_list[i]):
		plt.title(rxn_list[i])
		plt.hist(points_B[:, i], bins=np.linspace(np.min(np.vstack((points_B[:, i], points_A[:, i]))),
		                                          np.max(np.vstack((points_B[:, i], points_A[:, i]))), 50), alpha=0.5)
		plt.hist(points_A[:, i], bins=np.linspace(np.min(np.vstack((points_B[:, i], points_A[:, i]))),
		                                          np.max(np.vstack((points_B[:, i], points_A[:, i]))), 50), alpha=0.5)
		plt.show()
input()
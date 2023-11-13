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

HR_Points_location_high_A = Path("Data/HR/recon_1b_t_cells_gene_Day_2_single_B6_1_gene_therm")
HR_Points_location_high_B = Path("Data/HR/recon_1b_t_cells_gene_Day_2_single_TC_5")
points_A = np.array([])
distaces = np.array([])
start = 0
max_distance = 100
values = 10

point_sample_seperation = np.arange(1,max_distance+1)
for file in os.listdir(HR_Points_location_high_A):
	print(file)
	new_points = np.load(HR_Points_location_high_A / file)
	print(np.arange(100)[:11:2])
	print(np.shape(new_points))
	#input()
	#print(test[-1])
	
	print(np.shape(new_points))
	point_distances_A = []
	for i in point_sample_seperation:
		test_vals = np.arange(np.shape(new_points)[0])[:i * values+1:i]
		#print(test_vals)
		points_for_distances = new_points[start:start+i * values+1:i]
		#print(np.shape(points_for_distances))
		#input()
		dist_B = 0
		dist_A = 0
		dists_A = np.zeros(values)
		for j in range(np.size(dists_A)):
			# dist_orig += np.linalg.norm(n_point_orig[j]-n_point_orig[j+1])
			# dist_lc += np.linalg.norm(n_point_lc[j] - n_point_lc[j + 1])
			dists_A[j] += np.linalg.norm(points_for_distances[j] - points_for_distances[j + 1])
		print(dists_A,np.sum(dists_A),point_sample_seperation[i-1])
		point_distances_A.append(np.average(dists_A[:]))
	print(point_distances_A)
	print(np.shape(points_for_distances))
	plt.plot(np.arange(len(point_distances_A)),point_distances_A)
plt.show()
print("A")
input()
for file in os.listdir(HR_Points_location_high_B):
	points_B = np.load(HR_Points_location_high_B / file)
	
flux_model_dict = {}
for filename in os.listdir(HR_model_location):
	print(filename)
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_fast_key_json_model(HR_model_location / filename)
	flux_model_dict[filename.split(".")[0]] = recon_flux_model

lb, ub, S, b, rxn_list, met_list = flux_model_dict['B6-1'].dicts_to_mats()
print(rxn_list)
print(rxn_list.index("TKT2"))

rxn_list.append("Total_Flux")
print(len(rxn_list))

point_sample_seperation = np.arange(1,2000)
point_distances_orig = []
point_distances_A = []
point_distances_B = []
#distance_matrix = np.zeros((100-1,np.shape(points_A)[1]))
for i in point_sample_seperation:
	#n_point_orig = points_orig[0:i * 100:i]
	n_points_A = points_A[:i*5:i]
	n_points_B = points_B[:i * 5:i]
	#n_point_lc = points_lc[0:i*100:i]
	dist_orig = 0
	dist_B = 0
	dist_A = 0
	for j in range(5-1):
		#dist_orig += np.linalg.norm(n_point_orig[j]-n_point_orig[j+1])
		#dist_lc += np.linalg.norm(n_point_lc[j] - n_point_lc[j + 1])
		dist_A += np.linalg.norm(n_points_A[j]-n_points_A[j+1])
		dist_B += np.linalg.norm(n_points_B[j] - n_points_B[j + 1])
		#distance_matrix[j] = np.abs(n_point_orig[j] - n_point_lc[j + 1])
	#avg_distances = np.average(distance_matrix,axis=0)
	#print(np.shape(avg_distances))
	#print(np.max(avg_distances))
	#input()
	#point_distances_orig.append(dist_orig/(100-1))
	#point_distances_lc.append(dist_lc / (100 - 1))
	point_distances_A.append(dist_A/(10-1))
	point_distances_B.append(dist_B / (10 - 1))
#plt.plot(point_sample_seperation,point_distances_orig)
#plt.plot(point_sample_seperation,point_distances_lc)
plt.ylabel("Average Euclidean Distance")
plt.xlabel("Iteration difference")
plt.plot(point_sample_seperation*500,point_distances_A)
plt.plot(point_sample_seperation*500,point_distances_B)
plt.show()

#plt.hist(points_orig[:,-1],alpha = 0.5,bins = np.linspace(149,149.5,20))
#plt.hist(points_lc[:,-1],alpha = 0.5,bins = np.linspace(149,149.5,20))
#plt.show()


#plt.plot(np.arange(np.size(points_orig[:,-1])),points_orig[:,-1])
#plt.plot(np.arange(np.size(points_lc[:,-1])),points_lc[:,-1])
plt.plot((np.arange(np.size(points_A[:,-1]))),points_A[:,-1])
plt.plot((np.arange(np.size(points_B[:,-1]))),points_B[:,-1])
plt.show()

print(np.shape(points_A))
points_A = points_A[::500]
points_B = points_B[::500]
print(np.shape(points_A))
print(points_A[:,-1])


plt.hist(points_A[:,-1],bins = np.linspace(min(np.min(points_A[:,-1]), np.min(points_A[:,-1])),max(np.max(points_A[:,-1]), np.max(points_A[:,-1])),50),alpha = 0.5)
plt.hist(points_B[:,-1],bins = np.linspace(min(np.min(points_A[:,-1]), np.min(points_A[:,-1])),max(np.max(points_B[:,-1]), np.max(points_B[:,-1])),50),alpha = 0.5)
#plt.hist(points_ss[:,-1],bins = np.linspace(149,max(np.max(points_orig[:,-1]), np.max(points_lc[:,-1]),np.max(points_ss[:,-1])),20),alpha = 0.5)
plt.show()
cd_list = []
for i in range(np.shape(points_B)[1]):
	variance_A = np.var(points_A[:, i])
	mean_A = np.mean(points_A[:,i])
	
	variance_B = np.var(points_B[:,i])
	mean_B = np.mean(points_B[:,i])
	print(st.ttest_ind(points_A[:,i],points_B[:,i]))
	cd_list.append((mean_A - mean_B) / np.sqrt((variance_A+variance_B)))
	print(mean_A,mean_B,np.sqrt((variance_A+variance_B)))
	print(i)
	print(variance_A / mean_A)
	#plt.title(rxn_list[i])
	#plt.hist(points_B[:,i],bins = np.linspace(np.min(np.vstack((points_B[:,i],points_A[:,i]))),np.max(np.vstack((points_B[:,i],points_A[:,i]))),50),alpha = 0.5)
	#plt.hist(points_A[:, i], bins = np.linspace(np.min(np.vstack((points_B[:,i],points_A[:,i]))),np.max(np.vstack((points_B[:,i],points_A[:,i]))),50), alpha=0.5)
	#plt.show()
print(cd_list)
cd_list = np.array(cd_list)
plt.scatter(np.log10(np.abs(cd_list)),np.log10(np.mean(points_A,axis=0)+np.mean(points_B,axis=0)))
plt.show()
cd_list_ind = np.flip(np.argsort(np.abs(cd_list)))
print(cd_list)
for i in cd_list_ind:
	print(cd_list[i],i)
	if not np.isnan(cd_list[i]):
		plt.title(rxn_list[i])
		plt.hist(points_B[:,i],bins = np.linspace(np.min(np.hstack((points_B[:,i],points_A[:,i]))),np.max(np.hstack((points_B[:,i],points_A[:,i]))),50),alpha = 0.5,label="therm A")
		plt.hist(points_A[:, i], bins = np.linspace(np.min(np.hstack((points_B[:,i],points_A[:,i]))),np.max(np.hstack((points_B[:,i],points_A[:,i]))),50), alpha=0.5,label = 'gene therm')
		plt.legend()
		plt.show()
input()
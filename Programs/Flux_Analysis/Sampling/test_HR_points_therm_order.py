import os
import pandas as pd
import numpy as np
from pathlib import Path
import scipy.io as sio
import time
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

HR_model_location = Path("Data/Models/pkl_models/recon_1b_t_cells/HR_ready_models")
HR_Points_location = Path("Data/HR/HR_Points_Distances/recon_1b_t_cells")
Internal_null_S_path = Path("Data/null_space_S/recon_1b_t_cells/Null_S_int.mat")
Therm_feas_mat_location = Path("Data/HR/Therm_test_mat")

flux_model_dict = {}
for filename in os.listdir(HR_model_location):
	print(filename)
	recon_flux_model = Flux_Balance_Model()
	
	recon_flux_model.load_pkl_model(HR_model_location / filename)
	flux_model_dict[filename.split("_")[0]] = recon_flux_model
	break

first_key = list(flux_model_dict.keys())[0]
lb, ub, pos_S, b, pos_rxn_list, met_list = flux_model_dict[first_key].dicts_to_mats()
simp_neg_ind, comp_neg_ind, comp_perm, cutter, exchange_cutter = flux_model_dict[first_key].precomp_positive_flux_point_to_bi_dir(
		pos_rxn_list, prune_specific=["biomass_reaction"])

point_dict = {}
for filename in os.listdir(HR_Points_location):
	if ".npy" in filename:
		point_dict[filename.split("_")[0].split(".")[0]] = np.load(HR_Points_location / filename,allow_pickle=True)[:,:-1]
		print(point_dict[filename.split("_")[0].split(".")[0]].dtype)
		print(point_dict.keys())
		break


NS = np.transpose(sio.loadmat(Internal_null_S_path)["NS"])
time_list = []
for flux_model_name in flux_model_dict.keys():
	feas_list = []
	start = time.time()
	flux_model = flux_model_dict[flux_model_name]
	flux_points = point_dict[flux_model_name]
	count = 0
	for i in flux_points:
		count+=1
		ent_flux_point = flux_model.positive_flux_point_to_bi_dir(i, simp_neg_ind, comp_neg_ind, comp_perm, cutter,
		                                                          exchange_cutter=exchange_cutter)
		feas_list.append(3-flux_model.generate_therm_model(NS,ent_flux_point))
		if count % 100 == 0:
			if count == 50000:
				print(time.time()-start)
				input()
			print(count/np.shape(flux_points)[0],count)
	feas_array = np.array(feas_list)
	print(np.average(feas_array))
	if np.max(feas_array) > 1:
		print("ERROR")
	elif np.min(feas_array) < 0:
		print("ERROR")
	time_list.append(time.time() - start)
	np.save(Therm_feas_mat_location/flux_model_name,feas_array)
print(time_list)
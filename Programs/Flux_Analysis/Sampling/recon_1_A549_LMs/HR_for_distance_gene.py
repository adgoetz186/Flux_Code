import os
import sys
import pandas as pd
import scipy.io as sio
import numpy as np
from pathlib import Path
import time


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
from Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
# _____ Setting the CWD to be Flux_Code END _____

HR_model_location = Path("Data/Models/json_models/fast_key_format/recon_1_t_cells_12_11_23/HR_Ready_Day_2_single_size")
warmup_sample_location = Path("Data/HR/HR_Warmup/recon_1_t_cells_gene_Day_2_single_size_12_11_23")
HR_Points_location = Path("Data/HR/HR_Points_Distances/recon_1_t_cells_gene_Day_2_single_size_12_11_23")

Internal_nullspace_S = sio.loadmat(Path("Data/null_space_S/recon_1b_t_cells_single_size/Null_S_internal.mat"))['Null_S_internal']

name_gene_dict = {"B6-1":"B61-D2-GEX-B5","TC-5":"TC5-D2-GEX-B1"}

flux_model_dict = {}
for filename in os.listdir(HR_model_location):
	print(filename)
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_fast_key_json_model(HR_model_location / filename)
	flux_model_dict[filename.split(".")[0]] = recon_flux_model

warmup_dict = {}
for filename in os.listdir(warmup_sample_location):
	warmup_dict[filename.split("_")[-1].split(".")[0]] = np.load(warmup_sample_location / filename)
	print(np.shape(np.load(warmup_sample_location / filename)))
	print(np.load(warmup_sample_location / filename))


print(flux_model_dict)
time_list = []
for flux_model_name in name_gene_dict.keys():
	gene_array = np.load(f"Data/scRNA/Bulk/{name_gene_dict[flux_model_name]}/matrix.npy")
	feature_array = np.load(f"Data/scRNA/Bulk/{name_gene_dict[flux_model_name]}/features.npy")
	gene_array = gene_array/np.average(gene_array)
	print(np.average(gene_array))
	start = time.time()
	flux_model = flux_model_dict[flux_model_name]
	RAS = np.array(flux_model.create_RAS_values(gene_array,feature_array,False))
	RPS = flux_model.RAS_to_RPS_mat(RAS)
	thermo_const = {"prune_specific": ["biomass_reaction"], "NS_internal": Internal_nullspace_S}
	#point_load = np.load(HR_Points_location_load / (flux_model_name + "_HR_points.npy"))
	# might wanna add a test here to make sure names line up, right now just test for complete
	# 100000, 500
	#HR_samples = flux_model.HRSampler_gene_bias_lincomb_pinch(warmup_dict[flux_model_name], 100000, 500,RPS,1.1,gene_penalty_mod=1e-10)
	HR_samples = flux_model.HRSampler_gene_bias_lincomb_pinch(warmup_dict[flux_model_name], 100000, 500,RPS,1.1,gene_penalty_mod=1e-10)

	#HR_samples = flux_model.HRSampler_lincom(warmup_dict[flux_model_name], 10000, 500)
	time_list.append(time.time() - start)
	print(time_list)
	np.save(HR_Points_location / (flux_model_name + "_HR_points"), HR_samples)
print(time_list)
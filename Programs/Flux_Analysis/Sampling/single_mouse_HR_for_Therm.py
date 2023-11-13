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
# _____ Setting the CWD to be Flux_Code END _____

HR_model_location = Path("Data/Models/pkl_models/recon_1b_t_cells/HR_ready_models")
warmup_sample_location = Path("Data/HR/HR_Warmup/recon_1b_t_cells")
HR_Points_location = Path("Data/HR/HR_Therm_sample/recon_1b_t_cells")
Internal_null_S_path = Path("Data/null_space_S/recon_1b_t_cells/Null_S_int.mat")

flux_model_dict = {}
for filename in os.listdir(HR_model_location):
	print(filename)
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_pkl_model(HR_model_location / filename)
	flux_model_dict[filename.split("_")[0]] = recon_flux_model
	
warmup_dict = {}
for filename in os.listdir(warmup_sample_location):
	warmup_dict[filename.split("_")[-1].split(".")[0]] = np.load(warmup_sample_location/filename)
first_key = list(flux_model_dict.keys())[0]
lb, ub, pos_S, b, pos_rxn_list, met_list = flux_model_dict[first_key].dicts_to_mats()
simp_neg_ind, comp_neg_ind, comp_perm, cutter, exchange_cutter = flux_model_dict[first_key].precomp_positive_flux_point_to_bi_dir(
		pos_rxn_list, prune_specific=["biomass_reaction"])
NS = np.transpose(sio.loadmat(Internal_null_S_path)["NS"])

print(warmup_dict)


print(flux_model_dict)
time_list = []
for flux_model_name in flux_model_dict.keys():
	start = time.time()
	flux_model = flux_model_dict[flux_model_name]
	# might wanna add a test here to make sure names line up, right now just test for complete
	# 100000, 500
	HR_samples = flux_model.HRSampler_therm(warmup_dict[flux_model_name],1000, 15000,10,NS,simp_neg_ind,comp_neg_ind,comp_perm,cutter,exchange_cutter=exchange_cutter)
	time_list.append(time.time()-start)

	np.save(HR_Points_location/(flux_model_name+"_HR_points"),HR_samples)
print(time_list)
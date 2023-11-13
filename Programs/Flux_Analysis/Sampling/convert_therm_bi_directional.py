import copy

import numpy as np
import matplotlib.pyplot as plt
import time
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import os
import pickle
import pandas as pd
import scipy.io as sio
from pathlib import Path
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

input_path = Path("Data/HR/HR_Therm_sample/recon_1b_t_cells")
HR_model_location = Path("Data/Models/pkl_models/recon_1b_t_cells/HR_ready_models")
output_npy_path = Path("Data/HR/HR_Therm_bi_directional/recon_1b_t_cells")
output_mat_path = Path("Data/HR/HR_mat/HR_Therm_bi_directional/recon_1b_t_cells/therm_model.mat")





# loads flux models

recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_pkl_model(HR_model_location / os.listdir(HR_model_location)[0])

p_lb, p_ub, p_S, p_b, p_rxn_list, p_met_list = recon_flux_model.dicts_to_mats()
simp_neg_ind, comp_neg_ind, comp_perm, cutter, exchange_cutter = recon_flux_model.precomp_positive_flux_point_to_bi_dir(p_rxn_list,prune_exchange = False)

recon_flux_model.convert_model_to_bidirectional_flux()

# loads reaction names
lb, ub, S, b, rxn_list, met_list = recon_flux_model.dicts_to_mats()
rxn_list.append("total_flux")
names = rxn_list
print(names)






flux_dict = {}
flux_dict["rxn_names"] = names[:-1]

B6_list = []
TC_list = []
# generates dtf
for file_name in os.listdir(input_path):
	flux_samples = np.load(input_path / file_name)
	flux_sample_list = []
	for i in flux_samples:
		flux_sample_list.append(recon_flux_model.positive_flux_point_to_bi_dir(i[:-1],simp_neg_ind, comp_neg_ind, comp_perm, cutter))
	flux_samples = np.vstack(flux_sample_list)
	#flux_samples = np.average(flux_samples,axis = 0)
	np.save(output_npy_path/file_name,flux_samples)
	flux_dict[file_name.split(".")[0].split("_")[0].replace("-","_")] = flux_samples
sio.savemat(output_mat_path,flux_dict)
print(flux_dict)





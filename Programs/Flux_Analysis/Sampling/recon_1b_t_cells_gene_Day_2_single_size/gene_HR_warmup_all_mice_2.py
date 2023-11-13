import os
import pandas as pd
import numpy as np
from pathlib import Path
import scipy.io as sio
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

# Remove
Internal_nullspace_S = sio.loadmat(Path("Data/null_space_S/recon_1b_t_cells_single_size/Null_S_internal.mat"))['Null_S_internal']
trans_NS = np.transpose(Internal_nullspace_S)


positive_min_model_location = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/edited_pos_min_models_Day_2_single_size")
# it should be possible to eventually just have HR ready models after exp_aligned
# This way breaks things up to save intermediates
HR_ready_model_location = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/HR_Ready_Day_2_single_size_all")
warmup_sample_location = Path("Data/HR/HR_Warmup/recon_1b_t_cells_gene_Day_2_single_size_all")

flux_model_dict = {}
for filename in os.listdir(positive_min_model_location):
	print(filename)
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_fast_key_json_model(positive_min_model_location / filename)
	flux_model_dict[recon_flux_model.model_dict["model_name"]] = recon_flux_model

print(flux_model_dict)

for flux_model_name in flux_model_dict.keys():
	
	print(flux_model_name.split("-"))
	if (int(flux_model_name.split("-")[-1]) > 2) and (int(flux_model_name.split("-")[-1]) <= 4):
		flux_model = flux_model_dict[flux_model_name]
		# 10000
		# 5
		warmup = flux_model.generate_warmup_gb_pinch(10000,[1,1.005,1.01], repeat_basis= 5,points_save_path = warmup_sample_location,model_save_path = HR_ready_model_location,trans_NS = trans_NS)

	
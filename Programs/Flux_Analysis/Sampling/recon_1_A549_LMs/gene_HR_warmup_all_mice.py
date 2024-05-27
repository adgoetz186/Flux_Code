import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
import scipy.io as sio


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
from Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
# _____ Setting the CWD to be Flux_Code END _____

# Remove
Internal_nullspace_S = sio.loadmat(Path("Data/null_space_S/recon_1_A549_4_30_2024/NS.mat"))['NS']
trans_NS = np.transpose(Internal_nullspace_S)


positive_min_model_location = Path("Data/Models/json_models/fast_key_format/recon_1_A549_4_30_2024/CurVar_pos_min_model.json")
# it should be possible to eventually just have HR ready models after exp_aligned
# This way breaks things up to save intermediates
HR_ready_model_location = Path("Data/Models/json_models/fast_key_format/recon_1_A549_4_30_2024")
warmup_sample_location = Path("Data/HR/HR_Warmup/Recon_1_A549_4_30_2024")


recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_fast_key_json_model(positive_min_model_location)
recon_flux_model.model_dict["model_name"] = "HR_ready_model"
lb, ub, S, b, rxn_list, met_list = recon_flux_model.dicts_to_mats()
recon_flux_model.reaction_info('biomass_reaction')
recon_flux_model.reaction_info('EX_val_L(e)_inverted')

warmup = recon_flux_model.generate_warmup_gb_pinch(10000,[1,2,4], repeat_basis= 2,points_save_path = warmup_sample_location,model_save_path = HR_ready_model_location,trans_NS = trans_NS)
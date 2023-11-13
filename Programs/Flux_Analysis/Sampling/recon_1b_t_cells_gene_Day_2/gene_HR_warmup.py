import os
import pandas as pd
import numpy as np
from pathlib import Path
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



positive_min_model_location = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/pos_min_model_Day_2")
# it should be possible to eventually just have HR ready models after exp_aligned
# This way breaks things up to save intermediates
HR_ready_model_location = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/HR_Ready_Day_2")
warmup_sample_location = Path("Data/HR/HR_Warmup/recon_1b_t_cells_gene_Day_2")
penalty_location = Path("Data/scRNA/Processed/RPS_npy")

flux_model_dict = {}
for filename in os.listdir(positive_min_model_location):
	print(filename)
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_fast_key_json_model(positive_min_model_location / filename)
	flux_model_dict[recon_flux_model.model_dict["model_name"]] = recon_flux_model
	
print(flux_model_dict)

for flux_model_name in flux_model_dict.keys():
	print(flux_model_name)
	penalty_array = np.zeros(1)
	for i in os.listdir(penalty_location):
		if flux_model_name.replace("-","") in i and "D2" in i:
			penalty_array = np.load(penalty_location/i)
			print(penalty_array)
	flux_model = flux_model_dict[flux_model_name]
	warmup = flux_model.generate_warmup_gb(10000,[1,1.1,1.2],warmup_sample_location,HR_ready_model_location,rxn_penalty_vec = penalty_array)


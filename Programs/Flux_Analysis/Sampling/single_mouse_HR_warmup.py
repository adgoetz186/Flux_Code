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



positive_min_model_location = Path("Data/Models/pkl_models/recon_1b_t_cells/edited_pos_min_models")
# it should be possible to eventually just have HR ready models after exp_aligned
# This way breaks things up to save intermediates
HR_ready_model_location = Path("Data/Models/pkl_models/recon_1b_t_cells/HR_ready_models")
warmup_sample_location = Path("Data/HR/HR_Warmup/recon_1b_t_cells")

flux_model_list = []
for filename in os.listdir(positive_min_model_location):
	print(filename)
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_pkl_model(positive_min_model_location / filename)
	flux_model_list.append(recon_flux_model)

for flux_model_ind in range(len(flux_model_list)):
	flux_model = flux_model_list[flux_model_ind]
	warmup = flux_model.generate_warmup_gb(10000,[1,1.005,1.01],warmup_sample_location,HR_ready_model_location)


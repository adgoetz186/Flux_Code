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

min_model_location = Path("Data/Models/pkl_models/recon1_t_cell/min_models/")
positive_min_model_location = Path("Data/Models/pkl_models/recon1_t_cell/p_min_models")

# Potentially make sure all reactions the same for all cells
# check to make sure updates dont impact any reaction or metabolite infos


flux_model_list = []
for filename in os.listdir(min_model_location):
	print(filename)
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_pkl_model(min_model_location / filename)
	recon_flux_model.convert_model_to_positive_flux()
	recon_flux_model.save_model_as_pkl(positive_min_model_location / (filename.split("_")[0]+"_positive_min_model"))

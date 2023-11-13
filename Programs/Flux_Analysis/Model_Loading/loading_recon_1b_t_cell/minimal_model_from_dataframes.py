import os
import pandas as pd
import numpy as np
from pathlib import Path
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.FMC_Object_Functions as fof

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

exp_align_folder = Path("Data/Models/pkl_models/recon_1b_t_cells/exp_aligned/")
min_model_matrix_location = Path("Data/Minimal_Model_Data/recon_1b_t_cells/")
min_model_location = Path("Data/Models/pkl_models/recon_1b_t_cells/pos_min_models/")




flux_model_list = []
for filename in os.listdir(exp_align_folder):
	print(filename)
	
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_pkl_model(exp_align_folder / filename)
	print(recon_flux_model.model_name)
	flux_model_list.append(recon_flux_model)

min_list = fof.minimal_flux_list_multimodel(flux_model_list,min_model_matrix_location)
pos_min_list = fof.make_uniform_pos(min_list)
for model_ind in range(len(pos_min_list)):
	pos_min_list[model_ind].save_model_as_pkl(min_model_location/(pos_min_list[model_ind].model_name+"_pos_min_model"))

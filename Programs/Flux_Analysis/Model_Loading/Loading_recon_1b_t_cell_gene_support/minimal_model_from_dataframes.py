import os
import random

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

exp_align_folder = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/exp_aligned_Day_2_single_size/")
min_model_matrix_location = Path("Data/Minimal_Model_Data/recon_1b_t_cells_gene_support_single_size/")
min_model_location = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/pos_min_model_Day_2_single_size/")




flux_model_list = []
for filename in os.listdir(exp_align_folder):
	print(filename)
	
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_fast_key_json_model(exp_align_folder / filename)
	flux_model_list.append(recon_flux_model)
print([i.model_dict["model_name"] for i in flux_model_list])
random.shuffle(flux_model_list)
print([i.model_dict["model_name"] for i in flux_model_list])
print(flux_model_list)

min_list = fof.minimal_flux_list_multimodel_model_differences(flux_model_list,min_model_matrix_location)
pos_min_list = fof.make_uniform_pos(min_list)
for model_ind in range(len(pos_min_list)):
	pos_min_list[model_ind].save_model_as_fast_key_json(min_model_location)

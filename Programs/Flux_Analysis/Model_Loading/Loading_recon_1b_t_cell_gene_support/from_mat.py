from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import os
from pathlib import Path
import numpy as np
import sys
import gurobipy as gp

# File for loading the model from an external datasource and saving it as a pkl file
# The new saved model can be loaded with load_pkl_model()
# Since save formats differ this will not be universal

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

mat_file_location = Path("Data/Models/mat_models/zeromodel.mat")

save_file_location = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/")

# You might should add subsystem converter to model

model_comp={"S":"S", "lb" : "lb", "ub" : "ub", "b" : "b","pinched_reactions" : None, "metabolite_names" : "mets", "reaction_names" : "rxns", "grRules" : "grRules", "genes" : "genes", "metabolite_comp" : "metFormulas"}
read_model_name = "model"
model_name_to_save = "recon_1b_raw"

recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_recon1_mat_model(mat_file_location,new_model_name = "Raw")

recon_flux_model.save_model_as_fast_key_json(save_file_location)

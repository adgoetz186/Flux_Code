from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import os
from pathlib import Path
import numpy as np
import sys
import json

print(os.listdir())
# Loads the recon model from the specified mat file and saves it as a dictionary of pickled objects.
# The new saved model can be loaded with load_pkl_model()

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

# File locations should end with "/"
# File names should not end with file specifier (.txt, .py, etc)
mat_file_location = Path("Data/Models/mat_models/recon2_2.mat")


output_file_location = Path("Data/Models/json_models/fast_key_format/recon2_2/raw.json")

# This program assumes certain standard names for various parts of the model:
#   standard names: Stochiometric matrix = "S", lower bounds = "lb", upper bounds = "ub",
#   rhs of Sv = b constraint = "b", dictionary of pinched reactions = "pinched_reactions",
#   list of reaction names = "reaction_names", gene reaction rule list = "grRules",
#   list of gene names = "genes", list of metabolite names = "metabolite_names"
#   list of metabolite compositions = "metabolite_comp"
# If a loaded model uses different names, the loader wont know how to handle the alternative names
# This dictionary serves as a translator, the keys should be the standard names given above and the values should be
# the corresponding keys in the mat file. Running the program with no model_header_dict will provide a list of
# mat file keys.
model_name_to_save = "Recon2"

recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_recon2_2_mat_model(mat_file_location)

for key in recon_flux_model.rxn_dict.keys():
	if "bio" in key:
		print(key)
recon_flux_model.reaction_info("biomass_reaction")
biomass_names = ["biomass_protein","biomass_DNA","biomass_RNA","biomass_carbohydrate","biomass_lipid","biomass_other"]
for bmn in biomass_names:
	recon_flux_model.reaction_info(bmn)

recon_flux_model.save_model_as_fast_key_json(output_file_location)
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import os
from pathlib import Path
import numpy as np
import sys

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
json_file_location = Path("Data/Models/json_models/Bigg_format/Recon3D.json")


output_file_location = Path("Data/Models/json_models/fast_key_format/Recon3D/")

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
model_comp={"S":"S", "lb" : "lower_bound", "ub" : "upper_bound", "b" : "b","pinched_reactions" : None,"subsystem":"subsystem", "metabolite_names" : "mets", "rxn_names" : "id", "grRules" : "gene_reaction_rule", "metabolite_comp" : "metFormulas"}
model_name = "Raw"

recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_bigg_json_model(json_file_location,new_model_name="Raw")
recon_flux_model.save_model_as_fast_key_json(output_file_location)




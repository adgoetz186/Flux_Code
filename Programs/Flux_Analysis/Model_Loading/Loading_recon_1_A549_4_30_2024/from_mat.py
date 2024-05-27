import os
import sys
from pathlib import Path

# File for loading the model from an external datasource and saving it as a fast key json file
# The new saved model can be loaded with load_fast_key_json_model (create empty instance first)
# Loading models might take some level of massaging

# Fast key json files are something we introduce, they are simlar to normal json files but rather than store
# reactions in lists, dictionaries are used. This is primarily done to keep things simple when removing reactions
# as there is no need to hunt down and remove each spot where the reaction is found.
# When order is important lists are generated from dicts, this is always sorted alphabetically which ensures
# 2 equivalent models will yield identical results

# _____ Setting the CWD to be Flux_Code BEGIN _____
# The path to the Flux_Code dir can be specified directly below
# If not there is an attempt to find it, this is useful for running on the cluster where cwd does not start in the
# directory of the code being run
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

# define file paths
mat_file_location = Path("Data/Models/mat_models/zeromodel.mat")
save_file_location = Path("Data/Models/json_models/fast_key_format/recon_1_A549_4_30_2024/")

# The keys should not be altered, the values should agree with the loaded model
model_comp={"S":"S", "lb" : "lb", "ub" : "ub", "b" : "b","pinched_reactions" : None, "metabolite_names" : "mets", "reaction_names" : "rxns", "grRules" : "grRules", "genes" : "genes", "metabolite_comp" : "metFormulas"}
read_model_name = "model"
model_name_to_save = "recon_1b_raw"

# Load model from mat
recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_recon1_mat_model(mat_file_location,new_model_name = "Raw")



# Save model as fast key json
recon_flux_model.save_model_as_fast_key_json(save_file_location)

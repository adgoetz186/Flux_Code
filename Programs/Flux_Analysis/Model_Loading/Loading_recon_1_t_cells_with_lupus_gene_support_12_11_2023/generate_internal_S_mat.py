import copy
import sys
import os
import time
import scipy.io as sio
import pandas as pd
import numpy as np
from pathlib import Path


# This is to be used in conjuction with a matlab script containing the single line "Null_S_internal = null(internal_S,"rational")"
# Ideally this should be automated or replaced with pure python way to do this, but for now this is what works
# (Sympy was taking WAY too long and numerical approaches were generating dense matricies which added considerable
# computational costs)

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
import Programs.Flux_Analysis.Classes_And_Functions.FMC_Object_Functions as fof
# _____ Setting the CWD to be Flux_Code END _____


positive_min_model_location = Path("Data/Models/json_models/fast_key_format/recon_1_t_cells_12_11_23/edited_pos_min_model_Day_2_single_size")
null_space_S_path = Path("Data/null_space_S/recon_1_t_cells_12_11_23")

# Assumes all models share the same S. Order of S is normalized by generating function.

flux_model_list = []
for filename in os.listdir(positive_min_model_location):
	print(filename)
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_fast_key_json_model(positive_min_model_location / filename)
	flux_model_list.append(recon_flux_model)

flux_model = flux_model_list[0]
flux_model_copy = copy.deepcopy(flux_model)

flux_model_copy.convert_model_to_bidirectional_flux()
lb, ub, S, b, rxn_list, met_list = flux_model_copy.dicts_to_mats()
print(np.shape(S))
rxn_array = np.array(rxn_list)

print(np.shape(S))
print(len(flux_model_copy.get_exchanged_reaction_info()))
internal_lb, internal_ub, internal_S, b, internal_rxn, met_list = flux_model_copy.generate_mat_wo_exchange(
	prune_specific=["biomass_reaction"])
print(np.shape(internal_S))
sio.savemat(null_space_S_path/"internal_S.mat",{"internal_S":internal_S})


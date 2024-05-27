import copy
import os
import sys
from pathlib import Path

import numpy as np
import scipy.io as sio

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

# _____ Setting the CWD to be Flux_Code END _____


positive_min_model_location = Path("Data/Models/json_models/fast_key_format/recon_1_A549_4_30_2024/CurVar_pos_min_model.json")
null_space_S_path = Path("Data/null_space_S/recon_1_A549_4_30_2024")

# Assumes all models share the same S. Order of S is normalized by generating function.


recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_fast_key_json_model(positive_min_model_location)

flux_model_copy = copy.deepcopy(recon_flux_model)

lb, ub, S, b, rxn_list, met_list = flux_model_copy.dicts_to_mats()
print(np.shape(S))

flux_model_copy.convert_model_to_bidirectional_flux()
flux_model_copy.reaction_info("biomass_reaction")
print(flux_model_copy.test_feasibility())
input()

lb, ub, S, b, rxn_list, met_list = flux_model_copy.dicts_to_mats()
print(np.shape(S))
rxn_array = np.array(rxn_list)

print(np.shape(S))
print(len(flux_model_copy.get_exchanged_reaction_info()))
internal_lb, internal_ub, internal_S, b, internal_rxn, met_list = flux_model_copy.generate_mat_wo_exchange(
	prune_specific=["biomass_reaction"])
print(np.shape(internal_S))
print(flux_model_copy.test_feasibility())
input()
sio.savemat(null_space_S_path/"internal_S.mat",{"internal_S":internal_S})


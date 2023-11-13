import copy
import os
import time
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
		path_to_FC = Path.cwd().parents[
			[Path.cwd().parents[i].parts[-1] for i in range(len(Path.cwd().parts) - 1)].index(
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

positive_min_model_location = Path("Data/Models/pkl_models/recon_1b_t_cells/pos_min_models/")
# it should be possible to eventually just have HR ready models after exp_aligned
# This way breaks things up to save intermediates
null_space_S_path = Path("Data/null_space_S/recon_1b_t_cells")

flux_model_list = []
for filename in os.listdir(positive_min_model_location):
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_pkl_model(positive_min_model_location / filename)
	flux_model_list.append(recon_flux_model)

flux_model = flux_model_list[0]

lb, ub, pos_S, b, pos_rxn_list, met_list = flux_model.dicts_to_mats()

simp_neg_ind, comp_neg_ind, comp_perm, cutter, exchange_cutter = flux_model.precomp_positive_flux_point_to_bi_dir(pos_rxn_list,prune_specific=["biomass_reaction"])
test_bi_S = flux_model.positive_S_to_bi_dir(pos_S,simp_neg_ind,cutter,exchange_cutter=exchange_cutter)
print(np.shape(test_bi_S))
print("start")
start = time.time()
NS = flux_model.generate_null_space(test_bi_S)
print(time.time() - start)


np.save(null_space_S_path / "Nullspace", NS)


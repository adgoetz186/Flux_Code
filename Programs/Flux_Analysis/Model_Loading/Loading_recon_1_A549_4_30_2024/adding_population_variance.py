import os
import sys
import numpy as np
import copy as cp
import pandas as pd
from pathlib import Path

# Tweaks minimal model with biological knowledge

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

# removes reactions that cant carry flux
# removes curated reactions:
# ALR,
measured_rates_location = Path("Data/experimental_alignment_data/recon_1_A549/Measured_Rates.txt")
min_model_file_location = Path("Data/Models/json_models/fast_key_format/recon_1_A549_4_30_2024/Curated_pos_min_model.json")
output_file_location = Path("Data/Models/json_models/fast_key_format/recon_1_A549_4_30_2024")



recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_fast_key_json_model(min_model_file_location)

list_of_pinched_exchange = pd.read_csv(measured_rates_location)["reaction"].to_list()
print(list_of_pinched_exchange)
input()
lb, ub, S, b, rxn_list, met_list = recon_flux_model.dicts_to_mats()
recon_model_MK1 = cp.deepcopy(recon_flux_model)
for i in range(len(rxn_list)):
    print(rxn_list[i],lb[i],ub[i])
input()
for i in list_of_pinched_exchange:
    recon_model_MK1.update_reaction_bounds(i, min(lb[rxn_list.index(i)]*2,lb[rxn_list.index(i)]*.5), max(ub[rxn_list.index(i)]*2,ub[rxn_list.index(i)]*.5))
    recon_model_MK1.reaction_info(i)
biomass_rxn = (lb[rxn_list.index("biomass_reaction")] + ub[rxn_list.index("biomass_reaction")])/2
CV = np.sqrt(5.57) / 8.6
# https://www.spandidos-publications.com/10.3892/etm.2019.7873
pinch_dt = np.log(2) / biomass_rxn
sqrt_dt = pinch_dt / 8.6 * np.sqrt(5.57)
recon_model_MK1.update_reaction_bounds("biomass_reaction", np.log(2) / (pinch_dt + sqrt_dt),
										np.log(2) / (pinch_dt - sqrt_dt))
recon_model_MK1.reaction_info("biomass_reaction")
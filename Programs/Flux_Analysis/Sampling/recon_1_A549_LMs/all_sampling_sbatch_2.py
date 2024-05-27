import os
import sys
import pandas as pd
import numpy as np
import scipy.io as sio
from pathlib import Path
import time
import gurobipy as gp

# _____ Setting the CWD to be Flux_Code BEGIN _____
# Flux_Code path here
if not Path(os.getcwd()).parts[-1] == "Flux_Code":
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

from Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
# _____ Setting the CWD to be Flux_Code END _____

Internal_nullspace_S = sio.loadmat(Path("Data/null_space_S/recon_1_A549_4_30_2024/NS.mat"))['NS']
print(str(round(time.time()*1000)))


output_path = Path('Data/HR/HR_A549_therm_test_uniform')

recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_fast_key_json_model(Path('Data/Models/json_models/fast_key_format/recon_1_A549_4_30_2024/HR_ready_model.json'))


lb, ub, S, b, rxn_list, met_list = recon_flux_model.dicts_to_mats()
print(lb[rxn_list.index('biomass_reaction')],ub[rxn_list.index('biomass_reaction')])

RPS = np.zeros_like(lb)
RPS[rxn_list.index('biomass_reaction')] -= 1

# For the hypergator us this gurobi env key
#recon_flux_model.add_gp_key_env_to_model('grb-ts.ufhpc')

warmup_points = np.load(Path('Data/HR/HR_Warmup/Recon_1_A549_4_30_2024/warmup_HR_ready_model.npy'))



# might wanna add a test here to make sure names line up, right now just test for complete
# 10000, 5000
start = time.time()
HR_samples = np.array([])
#,thermo_const={"prune_specific": ["biomass_reaction"], "NS_internal": Internal_nullspace_S}
while np.size(HR_samples) == 0:
	HR_samples = recon_flux_model.HRSampler_gene_bias_lincomb_pinch(warmup_points, 2000, 1000,None,4,thermo_const={"prune_specific": ["biomass_reaction"], "NS_internal": Internal_nullspace_S})
#HR_samples = flux_model.HRSampler_lincom(warmup_dict[flux_model_name], 10000, 500)
print(HR_samples)

np.save(output_path / (f'2initial_{str(round(time.time()*1000))}'), HR_samples)
print(f"The time for this run was: {(time.time() - start)/60} minutes")
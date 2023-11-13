import os
import pandas as pd
import numpy as np
from pathlib import Path
import time
import scipy.io as sio
import matplotlib.pyplot as plt
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

HR_model_location = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/HR_ready_Day_2")
bi_directional_path = Path("Data/HR/HR_bi_directional_points/recon_1b_t_cells_gene_Day_2/")
mat_path = Path("Data/HR/HR_mat/recon_1b_t_cells_gene_Day_2/recon_1b_t_cells_gene_Day_2.mat")





# loads flux models
flux_model_list = []
for filename in os.listdir(HR_model_location):
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_fast_key_json_model(HR_model_location / filename)
	flux_model_list.append(recon_flux_model)

# loads flux models

recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_fast_key_json_model(HR_model_location / os.listdir(HR_model_location)[0])
recon_flux_model.convert_model_to_bidirectional_flux()

# loads reaction names
lb, ub, S, b, rxn_list, met_list = recon_flux_model.dicts_to_mats()
names = rxn_list
names.append("total_cost")
print(names)


mat_dict = {}
mat_dict["rxn_names"] = names
# generates dtf
for file_name in os.listdir(bi_directional_path):
	mat_dict[file_name.split(".")[0].replace("-","_")] = np.load(bi_directional_path / file_name)
	print(np.shape(np.load(bi_directional_path / file_name)))
sio.savemat(mat_path,mat_dict)
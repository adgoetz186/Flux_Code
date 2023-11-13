import copy

import numpy as np
import matplotlib.pyplot as plt
import time
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import os
import pickle
import pandas as pd
import scipy.io as sio
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

input_path = Path("Data/HR/HR_Points_analysis/recon_1b_t_cells/")
HR_model_location = Path("Data/Models/pkl_models/recon_1b_t_cells/HR_ready_models")
output_path = Path("Data/HR/HR_bi_directional_points/recon_1b_t_cells/")





# loads flux models
flux_model_list = []
for filename in os.listdir(HR_model_location):
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_pkl_model(HR_model_location / filename)
	flux_model_list.append(recon_flux_model)

# loads reaction names
lb, ub, S, b, rxn_list, met_list = flux_model_list[0].dicts_to_mats()
rxn_list.append("total_flux")
names = rxn_list


list_of_names = []

B6_dtf_list = []
TC_dtf_list = []
# generates dtf
for file_name in os.listdir(input_path):
	flux_samples = np.average(np.load(input_path / file_name),axis=0)
	flux_samples = np.reshape(flux_samples,(1,-1))
	print(np.shape(flux_samples))
	dtf = pd.DataFrame(flux_samples,columns=names)
	for rxn_name in dtf.columns:
		if "_inverted" in rxn_name:
			if rxn_name.replace("_inverted","") in dtf.columns:
				forward = rxn_name.replace("_inverted", "")
				dtf.loc[:,forward] = dtf.loc[:, forward].subtract(dtf.loc[:, rxn_name])
				dtf.pop(rxn_name)
			else:
				forward = rxn_name.replace("_inverted", "")
				dtf = dtf.rename(columns={rxn_name:forward})
				dtf.loc[:,forward] = -dtf.loc[:,forward]
	name_model = file_name.split("_")[0].split(".")[0]
	name_test = copy.deepcopy(dtf.columns).to_list()
	name_test.sort()
	if not np.alltrue(np.array([name_test[i] == dtf.columns[i] for i in range(len(name_test))])):
		print("Naming system seems corrupted, columns of resultant samples will not correspond to bi directional model")
	if "B6" in file_name:
		B6_dtf_list.append(dtf)
	elif "TC" in file_name:
		TC_dtf_list.append(dtf)
B6_dtf = pd.concat(B6_dtf_list)
TC_dtf = pd.concat(TC_dtf_list)
np.save(output_path/"B6_mice",B6_dtf.to_numpy())
np.save(output_path/"TC_mice",TC_dtf.to_numpy())



import numpy as np
import os
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
import json
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

raw_vert_sample_location = Path("Data/vert_samples/raw_vert_samples/recon1_t_cell")
positive_min_model_location = Path("Data/Models/pkl_models/recon1_t_cell/p_min_models")
processed_vert_sample_location = Path("Data/vert_samples/dataframe_verts/recon1_t_cell")

flux_model_list = []
for filename in os.listdir(positive_min_model_location):
	print(filename)
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_pkl_model(positive_min_model_location / filename)
	flux_model_list.append(recon_flux_model)

# Each model has the same reactions, inversion is what gives different counts
flux_model_rxns = [i.replace("[inverted]", "") for i in flux_model_list[0].reaction_names]
for fm in flux_model_list:
	flux_model_rxns_2 = [i.replace("[inverted]", "") for i in fm.reaction_names]
	for rxn in flux_model_rxns_2:
		if rxn not in flux_model_rxns:
			print(rxn)
print("done")

flux_model_rxn_names = [i.replace("[inverted]", "") for i in flux_model_list[0].reaction_names]
vert_list = []
for filename_ind in range(len(os.listdir(raw_vert_sample_location))):
	vert = np.load(raw_vert_sample_location / os.listdir(raw_vert_sample_location)[filename_ind])
	rxn_names = flux_model_list[filename_ind].reaction_names
	vert_dict = {}
	for i in flux_model_rxn_names:
		if i in rxn_names and (i + "[inverted]") not in rxn_names:
			vert_dict[i] = vert[:, rxn_names.index(i)]
		if i not in rxn_names and (i + "[inverted]") in rxn_names:
			vert_dict[i] = -vert[:, rxn_names.index((i + "[inverted]"))]
		if i in rxn_names and (i + "[inverted]") in rxn_names:
			vert_dict[i] = vert[:, rxn_names.index(i)] - vert[:, rxn_names.index((i + "[inverted]"))]
	# axes[filename_ind].hist(vert_dict['EX_o2(e)[inverted]'])
	with open(processed_vert_sample_location/(os.listdir(positive_min_model_location)[filename_ind].split("_")[0]+"_vert_dict"),"wb") as outfile:
		pickle.dump(vert_dict,outfile)

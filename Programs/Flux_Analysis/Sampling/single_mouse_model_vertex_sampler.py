import os
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


positive_min_model_location = Path("Data/Models/pkl_models/recon1_t_cell/min_models")
vert_sample_location = Path("Data/vert_samples/raw_vert_samples/recon1_t_cell")


flux_model_list = []
for filename in os.listdir(positive_min_model_location):
	print(filename)
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_pkl_model(positive_min_model_location / filename)
	print(np.shape(recon_flux_model.S))
	flux_model_list.append(recon_flux_model)

rxn_name_list = [i.reaction_names for i in flux_model_list]
for i in rxn_name_list:
	for j in rxn_name_list:
		for k in i:
			if k not in j:
				print(k)
	

for flux_model_ind in range(len(flux_model_list)):
	flux_model = flux_model_list[flux_model_ind]
	obj = flux_model.objective_function_generator(flux_model.reaction_names,[-1 for i in range(np.size(flux_model.reaction_names))])
	min_flux = -1*flux_model.find_objective_value(obj)
	max_allowed_flux = 1.1*min_flux
	flux_array = flux_model.generate_vertex_samples(100000,max_combined_flux=max_allowed_flux)
	np.save(vert_sample_location/(os.listdir(positive_min_model_location)[flux_model_ind].split("_")[0]+"_verts2"),flux_array)


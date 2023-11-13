import copy
import os
from functools import partial
import re
import pandas as pd
import numpy as np
from pathlib import Path
import scipy.io as sio
import time
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



pos_min_model = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/pos_min_model_Day_2")
warmup_sample_location = Path("Data/HR/HR_Warmup/recon_1b_t_cells")
RAS_npy_location = Path("Data/scRNA/Processed/RAS_npy")
Raw_gene_file_location = Path("Data/scRNA/Raw")
entrez_dict_location = Path("Data/Models/Gene_Data_For_Models/entrez_dict/RECON1.json")

with open(entrez_dict_location, ) as entrez:
	entrez_dict = json.load(entrez)
	

flux_model_dict = {}
for filename in os.listdir(pos_min_model):
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_fast_key_json_model(pos_min_model / filename)
	flux_model_dict[filename.split(".")[0]] = recon_flux_model

for folder in os.listdir(Raw_gene_file_location):
	features = np.loadtxt(Raw_gene_file_location / folder / "features.tsv.gz", delimiter="\t", dtype=str)
	gene_mat = np.array(sio.mmread(Raw_gene_file_location / folder / "matrix.mtx.gz").todense())
	print(folder, np.average(gene_mat))
	# gene_mat = np.load(scRNA_npy_location/(folder+".npy"))
	mdl_name = folder.split("-")[0][:2] + "-" + folder.split("-")[0][2:]
	model = flux_model_dict[mdl_name]
	gene_array = np.average(gene_mat, axis=1)
	RAS_dict = model.create_RAS_dict(gene_array,features)
	lb, ub, S, b, rxn_list, met_list = model.dicts_to_mats()
	RAS_array = np.array([RAS_dict[i] for i in rxn_list])
	
	np.save(RAS_npy_location / (folder + ".npy"), RAS_array)



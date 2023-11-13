import os
import pandas as pd
import json
import numpy as np
from pathlib import Path
import scipy.io as sio
import time
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

input_json_path = Path("Data/Models/json_models/RECON1.json")
HR_model_location = Path("Data/Models/pkl_models/recon_1b_t_cells/HR_ready_models")
gene_ortholog_location = Path("Data/scRNA/Ortholog_Data/human_mouse_hcop_six_column.txt.gz")
output_json_path = Path("Data/Models/Gene_Data_For_Models/entrez_dict/RECON1.json")

flux_model_dict = {}
for filename in os.listdir(HR_model_location):
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_pkl_model(HR_model_location / filename)
	flux_model_dict[filename.split("_")[0]] = recon_flux_model
first_model_ind = list(flux_model_dict.keys())[0]
model = flux_model_dict[first_model_ind]
used_gene_list = model.get_used_gene_list()

ortholog_matrix = np.loadtxt(gene_ortholog_location,dtype=str,delimiter="\t")

human_entrez_list = list(ortholog_matrix[1:,0])
mouse_ensembl_list = list(ortholog_matrix[1:,4])


human_entrez_dict = {}
for gene in used_gene_list:
	human_entrez_dict[gene] = {}
	count = 0
	for gene_2 in used_gene_list:
		if gene.split(".")[0] == gene_2.split(".")[0]:
			count += 1
	if gene.split(".")[0] in human_entrez_list:
		gene_convert_ind = human_entrez_list.index(gene.split(".")[0])
		human_entrez_dict[gene]["mouse_ensemble"] = mouse_ensembl_list[gene_convert_ind]
		human_entrez_dict[gene]["degen"] = count
		
	else:
		human_entrez_dict[gene]["degen"] = count

	
with open(output_json_path, "w") as mdl_file:
	json.dump(human_entrez_dict,mdl_file)
	

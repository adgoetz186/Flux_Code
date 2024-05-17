import os
from pathlib import Path
import numpy as np
import sys

# Adds reactions which are important for model behavior (GLUN, SERtm, FADRx)
# Updates some reaction bounds, clamps maintenance_ATP to be about 2
# finds the mouse ortholog id for the human genes given in recon1

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

# Defines paths
# path to raw model
raw_file_location = Path("Data/Models/json_models/fast_key_format/recon_1_t_cells_12_11_23/Raw.json")
# path to gene ortholog database
gene_ortholog_location = Path("Data/scRNA/Ortholog_Data/human_mouse_hcop_six_column.txt.gz")
# path to output model
output_file_location = Path("Data/Models/json_models/fast_key_format/recon_1_t_cells_12_11_23/")

# Loads model
recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_fast_key_json_model(raw_file_location)
recon_flux_model.model_dict["model_name"] = "Complete_prior"

# Updates reaction bounds
tol = 1e-7
recon_flux_model.update_reaction_bounds("maintenance_ATP",2-tol ,2+tol)
recon_flux_model.update_reaction_bounds("CLFORtex",0,0)
recon_flux_model.update_reaction_bounds("CITtbm","keep",0)
recon_flux_model.update_reaction_bounds("CK","keep",0)
recon_flux_model.update_reaction_bounds("biomass_reaction",0 ,200)

# Added to avoid feasible issues resulting from few reactions handling metabolite these metabolites
# Taken from RECON3D
recon_flux_model.model_dict["rxn_dict"]["FADRx"] = {
	'rxn_metabolites': {'h[c]': -1, 'fad[c]': -1, 'fadh2[c]': 1.0, 'nadh[c]': -1, "nad[c]": 1}, 'lb': 0, 'ub': 1000.0,
	'grRules': '', 'subSystems': 'Vitamin B2 metabolism',
	'rxnReferences': '',
	'confidenceScores': '', 'rxnECNumbers': '', 'rxnNotes': 'Added from Recon3D http://bigg.ucsd.edu/models/Recon3D/reactions/FADRx', 'rxnNames': 'FAD reductase'}

# Added to allow model reclaim observed excretion of L-Aspartate
recon_flux_model.update_reaction_bounds("ASPt6",-1000.0,'keep')

# updated gene rules for reactions, "HGNC:29570 or HGNC:4331"
recon_flux_model.model_dict["rxn_dict"]["GLUN"]["grRules"] = "(27165.1) or (2744.1)"
# "HGNC:16085"
recon_flux_model.model_dict["rxn_dict"]["SERtm"]["grRules"] = "(94081.1)"

# updates gene dict based on new gene rules
recon_flux_model.update_gene_dict()

# Recon1 uses human genes, our RNA genes are mouse based which requires conversion between the human entrez and
# stable Ensembl ID for mice, ENSMUSG###########. This is handled using an ortholog database provided by the HGNC:
# https://www.genenames.org/tools/hcop/
used_gene_list = recon_flux_model.model_dict["gene_dict"]
ortholog_matrix = np.loadtxt(gene_ortholog_location,dtype=str,delimiter="\t")
human_entrez_list = list(ortholog_matrix[1:,0])
mouse_ensembl_list = list(ortholog_matrix[1:,4])
for gene in used_gene_list.keys():
	if gene.split(".")[0] in human_entrez_list:
		if mouse_ensembl_list[human_entrez_list.index(gene.split(".")[0])][:7] == "ENSMUSG":
			used_gene_list[gene]["mouse_ortholog_id"] = mouse_ensembl_list[human_entrez_list.index(gene.split(".")[0])]
		else:
			used_gene_list[gene]["mouse_ortholog_id"] = "NO_MOUSE_ORTHOLOG_IN_HGNC"
	else:
		used_gene_list[gene]["mouse_ortholog_id"] = "GENE_NOT_FOUND_IN_HGNC_DATABASE"
recon_flux_model.save_model_as_fast_key_json(output_file_location)
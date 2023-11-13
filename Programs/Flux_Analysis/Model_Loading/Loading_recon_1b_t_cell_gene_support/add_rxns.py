from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import os
from pathlib import Path
import numpy as np

# Adds reactions which are important for model behavior (GLUN and SERtm)

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


raw_file_location = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/Raw.json")
gene_ortholog_location = Path("Data/scRNA/Ortholog_Data/human_mouse_hcop_six_column.txt.gz")
output_file_location = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/")


recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_fast_key_json_model(raw_file_location)
recon_flux_model.model_dict["model_name"] = "Complete_prior"

print(recon_flux_model.model_dict["rxn_dict"])

for rxn_name in recon_flux_model.model_dict["rxn_dict"].keys():
	if "atp" in rxn_name.lower():
		print(rxn_name)
for met_name in recon_flux_model.model_dict["met_dict"].keys():
	if "h2o" in met_name.lower():
		print(met_name)

recon_flux_model.reaction_info("DHCRD2")
recon_flux_model.reaction_info("FADH2tru")
recon_flux_model.reaction_info("FADH2tx")
recon_flux_model.reaction_info("ASPTAm")

recon_flux_model.reaction_info("EX_asp_L(e)")
recon_flux_model.metabolite_info("asp-L[e]")
recon_flux_model.metabolite_info("fad[c]")
recon_flux_model.metabolite_info("fadh2[c]")
print(recon_flux_model.model_dict["gene_dict"])

tol = 1e-7
recon_flux_model.update_reaction_bounds("maintenance_ATP",2-tol ,2+tol)
recon_flux_model.update_reaction_bounds("CLFORtex",0,0)
recon_flux_model.update_reaction_bounds("CITtbm","keep",0)
recon_flux_model.update_reaction_bounds("CK","keep",0)
recon_flux_model.update_reaction_bounds("biomass_reaction",0 ,200)

# Added to avoid feasible issues resulting from few reactions handling metabolite these metabolites
recon_flux_model.model_dict["rxn_dict"]["FADRx"] = {
	'rxn_metabolites': {'h[c]': -1, 'fad[c]': -1, 'fadh2[c]': 1.0, 'nadh[c]': -1, "nad[c]": 1}, 'lb': 0, 'ub': 1000.0,
	'grRules': '', 'subSystems': 'Vitamin B2 metabolism',
	'rxnReferences': '',
	'confidenceScores': '', 'rxnECNumbers': '', 'rxnNotes': 'Added from Recon3D http://bigg.ucsd.edu/models/Recon3D/reactions/FADRx', 'rxnNames': 'FAD reductase'}


#recon_flux_model.model_dict["rxn_dict"]["000"] = {'rxn_metabolites': {'h[c]': -1, 'fad[c]': -1, 'fadh2[c]': 1.0, 'nadh[c]': -1, "nad[c]": 1}, 'lb': 0, 'ub': 1000.0,'grRules': '(125061.1)', 'subSystems': 'Vitamin B2 metabolism','rxnReferences': '','confidenceScores': '', 'rxnECNumbers': '', 'rxnNotes': 'Added from Recon3D http://bigg.ucsd.edu/models/Recon3D/reactions/FADRx', 'rxnNames': 'FAD reductase'}

# Added to allow model reclaim observed excretion of L-Aspartate
recon_flux_model.update_reaction_bounds("ASPt6",-1000.0,'keep')

# "HGNC:29570 or HGNC:4331"
recon_flux_model.model_dict["rxn_dict"]["GLUN"]["grRules"] = "(27165.1) or (2744.1)"
# "HGNC:16085"
recon_flux_model.model_dict["rxn_dict"]["SERtm"]["grRules"] = "(94081.1)"


recon_flux_model.update_gene_dict()


#recon_flux_model.update_reaction_bounds("biomass_reaction",2-tol ,2+tol)
#recon_flux_model.reaction_info("biomass_reaction")
print(recon_flux_model.model_dict["gene_dict"])
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
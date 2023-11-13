from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import os
from pathlib import Path
import numpy as np

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
# File locations should end with "/"
# File names should not end with file specifier (.txt, .py, etc)
raw_file_location = Path("Data/Models/json_models/fast_key_format/Recon3D/Raw.json")

output_file_location = Path("Data/Models/json_models/fast_key_format/Recon3D/")


recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_fast_key_json_model(raw_file_location)
recon_flux_model.model_dict["model_name"] = "Rxn_added"


for met in recon_flux_model.model_dict["met_dict"].keys():
	
	if "gln" in met.lower():
		print(met)

for rxn in recon_flux_model.model_dict["rxn_dict"].keys():
	met = [i.lower() for i in recon_flux_model.model_dict["rxn_dict"][rxn]["rxn_metabolites"].keys()]
	if "gln__l_c" in met:
		print(rxn, met,recon_flux_model.model_dict["rxn_dict"][rxn])
	if "ser__l_c" in met and "ser__l_m" in met:
		print(rxn, met, recon_flux_model.model_dict["rxn_dict"][rxn])

# note SERtm is r1435 in recon3D

#Potentially remove HGNC:4331

for rxn_name in recon_flux_model.model_dict["rxn_dict"].keys():
	if rxn_name == "CLFORtex":
		print(rxn_name,recon_flux_model.model_dict["rxn_dict"][rxn_name])
# {'S': {'adp[c]': 1.0, 'atp[c]': -1.0, 'h2o[c]': -1.0, 'h[c]': 1.0, 'pi[c]': 1.0}, 'lb': 1.9999999, 'ub': 2.0000001, 'grRule': array([], dtype='<U1')}
# CK {'S': {'adp[m]': 1.0, 'atp[m]': -1.0, 'creat[m]': -1.0, 'pcreat[m]': 1.0}, 'lb': -1000.0, 'ub': 0, 'grRule': '(1159.1) or (1160.1)'}


tol = 1e-7
recon_flux_model.update_reaction_bounds("ATPM",2-tol ,2+tol)
recon_flux_model.update_reaction_bounds("BIOMASS_reaction",0.05-tol , 0.05+tol )
recon_flux_model.update_reaction_bounds("BIOMASS_maintenance",0 ,0 )
recon_flux_model.update_reaction_bounds("BIOMASS_maintenance_noTrTr",0 , 0 )
#recon_flux_model.update_reaction_bounds("CITtbm","keep",0)
#recon_flux_model.update_reaction_bounds("CLFORtex",0,0)

#recon_flux_model.update_reaction_bounds("CK","keep",0)


recon_flux_model.save_model_as_fast_key_json(output_file_location)
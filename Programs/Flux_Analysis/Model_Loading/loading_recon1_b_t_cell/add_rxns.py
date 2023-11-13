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


raw_file_location = Path("Data/Models/pkl_models/recon1_b_t_cell/raw.pkl")

output_file_location = Path("Data/Models/pkl_models/recon1_b_t_cell/rxn_added")


recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_pkl_model(raw_file_location)
tol = 1e-7
recon_flux_model.update_reaction_bounds("maintenance_ATP",2-tol ,2+tol)
recon_flux_model.reaction_info("biomass_reaction")

# "HGNC:29570 or HGNC:4331"
recon_flux_model.rxn_dict["GLUN"]["grRule"] = "(27165.1) or (2744.1)"
# "HGNC:16085"
recon_flux_model.rxn_dict["SERtm"]["grRule"] = "(94081.1)"


recon_flux_model.update_reaction_bounds("CITtbm","keep",0)
recon_flux_model.update_reaction_bounds("CK","keep",0)
recon_flux_model.update_reaction_bounds("CLFORtex",0,0)

recon_flux_model.save_model_as_pkl(output_file_location)
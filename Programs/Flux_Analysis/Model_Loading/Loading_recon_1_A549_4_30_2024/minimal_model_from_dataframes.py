import os
import sys
from pathlib import Path


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
import Programs.Flux_Analysis.Classes_And_Functions.FMC_Object_Functions as fof
# _____ Setting the CWD to be Flux_Code END _____

exp_align_location = Path("Data/Models/json_models/fast_key_format/recon_1_A549_4_30_2024/Experimentally_aligned.json")
min_model_matrix_location = Path("Data/Minimal_Model_Data/recon_1_A549_pos")

# min model location, make sure folder is created
min_model_location = Path("Data/Models/json_models/fast_key_format/recon_1_A549_4_30_2024")






recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_fast_key_json_model(exp_align_location)
recon_flux_model.model_dict['model_name'] = "Pos_min_model"
recon_flux_model.convert_model_to_positive_flux()

min_model = fof.minimal_flux_list_multimodel_model_differences([recon_flux_model],min_model_matrix_location)[0]
print(min_model.test_feasibility())
#remember to fun FVA after easing the bounds

min_model.save_model_as_fast_key_json(min_model_location)

import os
import sys
import time
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

model_file_path = Path("Data/Models/json_models/fast_key_format/recon_1_A549_4_30_2024/Experimentally_aligned.json")
output_model_folder = Path("Data/Minimal_Model_Data/recon_1_A549_pos")
recon1_exp_align = Flux_Balance_Model()
recon1_exp_align.load_fast_key_json_model(model_file_path)
recon1_exp_align.convert_model_to_positive_flux()
start = time.time()
recon1_exp_align.generate_essential_flux_dataframe(125,output_model_folder,"A549_pos_",print_progress = True,seed_name=True)
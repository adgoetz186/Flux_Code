import os
import time
import sys
from pathlib import Path
import numpy as np
import pandas as pd

# generates minimal model matrices. Binary matricies where each row describes a unique minimal model
# Columns correspond to reactions, entries are 1 when a reaction is essential and 0 otherwise.

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

# defines file locations
model_file_location = Path("Data/Models/json_models/fast_key_format/recon_1_t_cells_12_11_23/exp_aligned_Day_2_single_size")
output_model_folder = Path("Data/Minimal_Model_Data/recon_1_t_cells_12_11_23/")

# Finds minimal models for each condition, multiple should be run in parallel to generate enough models efficiently
for filename in os.listdir(model_file_location):
	model_file_path = model_file_location/filename
	recon1_exp_align = Flux_Balance_Model()
	recon1_exp_align.load_fast_key_json_model(model_file_path)
	start = time.time()
	recon1_exp_align.generate_essential_flux_dataframe(125,output_model_folder,recon1_exp_align.model_dict['model_name']+"_2_",print_progress = True,seed_name=True)
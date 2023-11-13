import os
import pandas as pd
import numpy as np
from pathlib import Path
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.FMC_Object_Functions as fof

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

min_model_location = Path("Data/Minimal_Model_Data/recon_1b_t_cells")


flux_model_list = []
name_list = []
for filename in os.listdir(min_model_location):
	if "header" in filename:
		names = np.load(min_model_location/filename)
	else:
		flux_model_list.append(np.load(min_model_location/filename))
total_essential = np.vstack(flux_model_list)
print(np.count_nonzero(np.where(np.average(total_essential,axis = 0)>=0.14525,np.zeros_like(np.average(total_essential,axis = 0)),np.ones_like(np.average(total_essential,axis = 0)))))
print(np.average(total_essential[:,np.where(names=="PGMT")]))
print(np.shape(total_essential))
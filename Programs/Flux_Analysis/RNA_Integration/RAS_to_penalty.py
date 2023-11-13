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

RAS_npy_location = Path("Data/scRNA/Processed/RAS_npy")
RPS_npy_location = Path("Data/scRNA/Processed/RPS_npy")


for RAS_file in os.listdir(RAS_npy_location):
	RAS_array = np.load(RAS_npy_location/RAS_file)
	min_RAS = np.nanmin(RAS_array[np.nonzero(RAS_array)[0]])
	max_RAS = np.nanmax(RAS_array[np.nonzero(RAS_array)[0]])
	median_RAS = np.nanmedian(RAS_array)
	edited_RAS = np.nan_to_num(RAS_array,nan=median_RAS)
	penalty_scores = 1/edited_RAS
	print(penalty_scores)
	np.save(RPS_npy_location/RAS_file,penalty_scores)
	
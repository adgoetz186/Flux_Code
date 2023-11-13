import os
import pandas as pd
import numpy as np
from pathlib import Path
import time
import matplotlib.pyplot as plt
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

HR_Points_location = Path("Data/HR/HR_Points_Distances/recon_1b_t_cells_gene_Day_2")

def distance(f):
	pass

samples = 1000
distance_per_saved = 500
for file_name in os.listdir(HR_Points_location):
	if file_name.split(".")[-1] == "npy":
		point_array = np.load(HR_Points_location/file_name)
		print(np.shape(point_array))
		dist_list = []
		for i in range(round(np.shape(point_array)[0]/samples-1)):
			avg_dist = 0
			for j in range(samples):
				avg_dist += np.linalg.norm(point_array[(i+1)*j]-point_array[(i+1)*(j+1)])
			dist_list.append(avg_dist/samples)
		print(dist_list)
		plt.plot(np.arange(round(np.shape(point_array)[0]/samples-1))*distance_per_saved,dist_list,label= file_name.split("_")[0])
plt.xlabel("steps between samples")
plt.ylabel("distance between points")
plt.legend()
plt.show()
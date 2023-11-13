import os
import pandas as pd
import numpy as np
import scipy.stats as st
from pathlib import Path
import matplotlib.pyplot as plt
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model

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

penalty_location = Path("Data/scRNA/Processed/RPS_npy")
HR_Points_location = Path("Data/HR/HR_Points_Distances/recon_1b_t_cells_gene_Day_2")


def sort_penalty(point,RPS):
	test = 1/RPS + np.random.random(np.size(RPS))/1e-100
	#return np.average(np.abs(st.rankdata(point)-st.rankdata(1/RPS)))
	return np.linalg.norm(st.rankdata(point)-st.rankdata(1/RPS))
	#return np.average(np.abs(st.rankdata(test) - st.rankdata(1 / RPS)))
	


print(st.rankdata([1,1,3]))

for file_name in os.listdir(HR_Points_location):
	print(file_name)
	if file_name.split(".")[-1] == "npy":
		penalty_array = np.zeros(1)
		for i in os.listdir(penalty_location):
			print(i)
			if file_name.replace("-", "").split("_")[0] in i and "D2" in i:
				penalty_array = np.load(penalty_location / i)
				point_array = np.load(HR_Points_location/file_name)
				error_terms = []
				for i in point_array:
					#error_terms.append(sort_penalty(i[:-1],penalty_array))
					error_terms.append(i[-1])
				plt.hist(error_terms)
				plt.xlim([min(error_terms),max(error_terms)])
				plt.show()
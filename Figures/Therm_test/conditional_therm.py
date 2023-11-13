import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
import scipy.io as sio
import time
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model

def con_step(feas_array,step):
	feas_ind = np.nonzero(feas_array)[0]
	#print(feas_ind)
	test_feas = feas_ind + step
	test_feas = test_feas[test_feas < (np.size(feas_array)-1)]
	infeas_ind = np.nonzero(1 - feas_array)[0]
	test_infeas = infeas_ind + step
	test_infeas = test_infeas[test_infeas < (np.size(feas_array) - 1)]
	return np.average(feas_array[test_feas]), 1-np.average(feas_array[test_infeas])

	
	

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

Therm_feas_mat_location = Path("Data/HR/Therm_test_mat")


for filename in os.listdir(Therm_feas_mat_location):
	feas_array = np.load(Therm_feas_mat_location / filename)
	prop_feas = []
	prop_infeas = []
	step_list = np.round(10**np.linspace(0,5,1000))
	for step in step_list:
		
		p_feas, p_infeas = con_step(feas_array,int(step))
		prop_feas.append(p_feas)
		prop_infeas.append(p_infeas)
	plt.plot(step_list,prop_feas)
	plt.plot(step_list,prop_infeas)
	plt.xlabel("n")
	plt.ylabel("prob $p_i$ is feasible\ngiven $p_{i-n}$ is feasible")
	plt.ylim([0,1])
	plt.show()

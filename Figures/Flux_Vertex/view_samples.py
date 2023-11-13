import numpy as np
import os
import pickle
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

vert_sample_location = Path("Data/vert_samples/dataframe_verts/recon1_t_cell")

vert_dicts = []
for filename in os.listdir(vert_sample_location):
	with open(vert_sample_location/filename,"rb") as readfile:
		vert_dicts.append(pickle.load(readfile))
print("done")

vert_list = []
for key in vert_dicts[0].keys():
	plt.hist(vert_dicts[0][key],bins = 100)
	plt.title(key)
	plt.show()
	#fig, axes = plt.subplots(8, 1, sharex=True)
	
	#for vert_ind in range(len(vert_dicts)):
		#axes[filename_ind].hist(vert_dict['EX_o2(e)[inverted]'])
		#axes[vert_ind].hist(vert_dicts[vert_ind][key],bins=100)
		#axes[filename_ind].hist(vert_dict['EX_nh4(e)'])
		#axes[filename_ind].hist(vert_dict['EX_co2(e)'])
		#axes[vert_ind].set_title(os.listdir(vert_sample_location)[vert_ind])
	#fig.suptitle(key)
	#plt.show()
	#plt.hist(vert_dict[1])
	
#possibly add zip function
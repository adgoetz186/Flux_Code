import os
import pandas as pd
import numpy as np
import scipy.sparse as sp
from pathlib import Path
from matplotlib import pyplot as plt
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

S_val = np.zeros((2,2))
#rxn
#S_val[0,0] = -1
#S_val[1,0] = -1
#S_val[2,0] = 2
#S_val[3,0] = 0

#S_val[2,1] = -1

#S_val[0,2] = 1
#S_val[1,2] = 1




ub = 2 * np.ones(2)
lb = 0 * np.ones(2)
flux_modelb = Flux_Balance_Model()
flux_modelb.create_model_from_components("Cube",S_val, lb,ub,np.zeros(2),[f"r_{i}" for i in range(2)],[f"m_{i}" for i in range(2)])


warmup = flux_modelb.generate_warmup_gb(50, [1,1.0001, 2],"","")
print(flux_modelb.model_dict["total_flux_limits"])
print("alright")

#samples_vs = flux_modelb.generate_vertex_samples(100)
print(warmup[0]-warmup[1])

samples_hr = flux_modelb.HRSampler_2_ss(warmup,10000,10000,)
print(samples_hr)
plt.hist(samples_hr[:,0])
plt.hist(samples_hr[:,1])
plt.hist(samples_hr[:,2])
plt.hist(samples_hr[:,3])
plt.show()
plt.hist(samples_hr[:,1])
plt.hist(samples_hr[:,0])
plt.show()
plt.hist(samples_hr[:,-1])
plt.show()
plt.hist(samples_hr[:,-2])
plt.show()
print(warmup)
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')


sequence_containing_x_vals_hr = list(samples_hr[:,0])
sequence_containing_y_vals_hr = list(samples_hr[:,1])
sequence_containing_z_vals_hr = list(samples_hr[:,2])
ax.scatter(sequence_containing_x_vals_hr, sequence_containing_y_vals_hr, sequence_containing_z_vals_hr)

sequence_containing_x_vals_vs = list(samples_vs[:,0])
sequence_containing_y_vals_vs = list(samples_vs[:,1])
sequence_containing_z_vals_vs = list(samples_vs[:,2])
ax.scatter(sequence_containing_x_vals_vs, sequence_containing_y_vals_vs, sequence_containing_z_vals_vs,s = 100,color = "red")
#ax.scatter(warmup[:,0], warmup[:,1], warmup[:,2],s = 50)
ax.set_xlabel('$X$', fontsize=20)
ax.set_ylabel('$Y$', fontsize=20)
ax.set_zlabel('$Z$', fontsize=20)
plt.show()
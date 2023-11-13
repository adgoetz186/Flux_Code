import scipy.io as sio
import numpy as np
from pathlib import Path
import os
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

uniform_path = Path("Data/HR/HR_Points_Distances/recon_1b_t_cells_gene_Day_2_single_B6_1/B6-1_HR_points.npy")
gene_bais_path = Path("Data/HR/HR_Points_Distances/recon_1b_t_cells_gene_Day_2_single_B6_1_gene/B6-1_HR_points.npy")
therm_and_gene_path = Path("Data/HR/HR_Points_Distances/recon_1b_t_cells_gene_Day_2_single_B6_1_gene_therm/B6-1_HR_points.npy")
therm_and_gene_path_B = Path("Data/HR/HR_Points_Distances/recon_1b_t_cells_gene_Day_2_single_B6_1_gene_therm_B/B6-1_HR_points.npy")

HR_model_location = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/HR_Ready_Day_2_single_size_B6_1_A")


flux_model_dict = {}
for filename in os.listdir(HR_model_location):
	print(filename)
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_fast_key_json_model(HR_model_location / filename)
	flux_model_dict[filename.split(".")[0]] = recon_flux_model
print(flux_model_dict["B6-1"])
file_dict = {"uniform":uniform_path,"gene_bias":gene_bais_path,"therm_and_gene_bias":therm_and_gene_path,"therm_and_gene_bias_B":therm_and_gene_path_B}

lb, ub, S, b, rxn_list, met_list = flux_model_dict['B6-1'].dicts_to_mats()
point_samples = {}
for key in file_dict.keys():
	points = np.load(file_dict[key])
	print(np.shape(points))
	#print(points)
	point_samples[key] = points
	
point_samples["rxn_names"] = np.array(rxn_list)
print(point_samples)
#sio.savemat("Point_Samples.mat", point_samples)
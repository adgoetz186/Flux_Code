import os
import time
import scipy.linalg as sla
import pandas as pd
import numpy as np
import copy
from pathlib import Path
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.FMC_Object_Functions as fof

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



positive_min_model_location = Path("Data/Models/pkl_models/recon_1b_t_cells/edited_pos_min_models")
# it should be possible to eventually just have HR ready models after exp_aligned
# This way breaks things up to save intermediates
null_space_S_path = Path("Data/null_space_S/recon_1b_t_cells/Nullspace.npy")

flux_model_list = []
for filename in os.listdir(positive_min_model_location):
	print(filename)
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_pkl_model(positive_min_model_location / filename)
	flux_model_list.append(recon_flux_model)


for model in flux_model_list:
	model_cpy = copy.deepcopy(model)
	print(0, model.rxn_dict["EX_gly(e)"])
	print(0, model.rxn_dict["EX_gly(e)_inverted"])
	model.convert_model_to_bidirectional_flux()
	print(1,model.rxn_dict["EX_gly(e)"])
	model.convert_model_to_positive_flux()
	print(2,model.rxn_dict["EX_gly(e)_inverted"])
	for rxn in model.rxn_dict.keys():
		if not fof.comp_rxn_dict(rxn,model,model_cpy):
			print(rxn, fof.comp_rxn_dict(rxn,model,model_cpy))
			print(model.rxn_dict[rxn],model_cpy.rxn_dict[rxn])
	print(0,model_cpy.rxn_dict["EX_gly(e)_inverted"])
	print(3,model.rxn_dict["EX_gly(e)_inverted"])
	input()
		

flux_model = flux_model_list[0]
print(flux_model.rxn_dict["maintenance_ATP"])
internal_lb, internal_ub, internal_S, b, internal_rxn, met_list = flux_model.generate_mat_wo_exchange(prune_specific=["biomass_reaction"])
print(np.shape(internal_S))

start = time.time()
test = np.arange(850)
for i in range(1000):
	
	ph = np.zeros_like(test)
	ph[0:-1] = test[1:]

print(ph)
print(time.time()-start)
input()
count = 0
NS = np.load(null_space_S_path,allow_pickle = True)
#NS = sla.null_space(internal_S)
NS = np.transpose(np.array(NS,dtype=float))
start = time.time()
for i in range(1000):
	point = (2*np.random.random(885)-1)*10
	feas = flux_model.generate_therm_model(NS,point)
	if feas == 2:
		count += 1
	else:
		print(feas)
print(time.time()-start)
print(count)
input()
print(np.shape(NS))
print(flux_model.generate_therm_model(NS,point,tol = 1e-9))
input()
print(NS)
print(list(NS))
print(np.transpose(NS))
print(np.sum(np.abs(np.matmul(internal_S,NS))))

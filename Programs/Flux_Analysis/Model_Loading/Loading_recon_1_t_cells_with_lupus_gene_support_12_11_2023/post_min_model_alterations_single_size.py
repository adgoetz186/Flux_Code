import os
import sys
import gurobipy as gp
from gurobipy import GRB
from pathlib import Path
import numpy as np

# Tweaks minimal model with biological knowledge

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

# removes reactions that cant carry flux
# removes curated reactions:
# ALR,

min_model_file_location = Path("Data/Models/json_models/fast_key_format/recon_1_t_cells_12_11_23/pos_min_model_Day_2_single_size/")
output_file_location = Path("Data/Models/json_models/fast_key_format/recon_1_t_cells_12_11_23/edited_pos_min_model_Day_2_single_size/")

flux_model_list = []
for filename in os.listdir(min_model_file_location):
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_fast_key_json_model(min_model_file_location / filename)
	flux_model_list.append(recon_flux_model)

for model in flux_model_list:
	print(model.test_feasibility())
	print(model.reaction_info("ALR"))
	del model.model_dict["rxn_dict"]["ALR"]
	
for model in flux_model_list:
	print(model.reaction_info("PGI_inverted"))
	del model.model_dict["rxn_dict"]["PGI_inverted"]
	
for model in flux_model_list:#
	print(model.reaction_info("SBTR"))
	del model.model_dict["rxn_dict"]["SBTR"]
	print(model.test_feasibility())
	print(model.model_dict.keys())
	print(model.dicts_to_mats())



#input()

#for model in flux_model_list:
#	print(model.reaction_info("G6PDH2r"))
#	del model.rxn_dict["G6PDH2r"]



#for flux_model in flux_model_list:
#	new_bounds = flux_model.fva()
#	for rxn_name in flux_model.rxn_dict.keys():
#		flux_model.update_reaction_bounds(rxn_name, new_bounds[rxn_name]["lb"], new_bounds[rxn_name]["ub"])

list_of_rxns_to_remove = []
for model in flux_model_list:
	list_of_rxns = list(model.model_dict["rxn_dict"])
	list_of_rxns_to_remove_model = []
	for rxn in list_of_rxns:
		if model.model_dict["rxn_dict"][rxn]["ub"] == model.model_dict["rxn_dict"][rxn]["lb"] and model.model_dict["rxn_dict"][rxn]["ub"] == 0:
			list_of_rxns_to_remove_model.append(rxn)
	list_of_rxns_to_remove.append(list_of_rxns_to_remove_model)

print(list_of_rxns_to_remove)


for i in list_of_rxns_to_remove[0]:
	keep = False
	for j in list_of_rxns_to_remove:
		if i not in j:
			keep = True
	if not keep:
		for model in flux_model_list:
			del model.model_dict["rxn_dict"][i]

for model in flux_model_list:
	rxn_bound_dict = model.fva()
	for rxn_name in model.model_dict["rxn_dict"]:
		model.update_reaction_bounds(rxn_name,rxn_bound_dict[rxn_name]['lb'],rxn_bound_dict[rxn_name]['ub'])

	



# write code to compare models

lb_0,ub_0,S_0,b_0,rxn_list_0,met_list_0 = flux_model_list[0].dicts_to_mats()
rxn_list_0 = np.array(rxn_list_0)
met_list_0 = np.array(met_list_0)
for model_ind in range(len(flux_model_list)):
	
	lb, ub, S, b, rxn_list, met_list = flux_model_list[model_ind].dicts_to_mats()
	

	#print(rxn_array[diff_ind])
	#print(point_to_use[diff_ind][-5:])
	#print(rxn_array[diff_ind][-5:])
	#print(lb[diff_ind][-5:])
	#print(ub[diff_ind][-5:])
	#print(diff)
	print(np.alltrue(rxn_list_0==rxn_list))
	print(np.alltrue(met_list_0 == met_list))
	print(flux_model_list[model_ind].test_feasibility())
	#print(np.alltrue(lb_0 == lb))
	#print(np.alltrue(ub_0 == ub))
	print(np.alltrue(S_0 == S))
	print(np.alltrue(b_0 == b))
	flux_model_list[model_ind].save_model_as_fast_key_json(output_file_location)






#recon_flux_model.update_reaction_bounds("CLFORtex",0,0)
#recon_flux_model.update_reaction_bounds("CITtbm","keep",0)
#recon_flux_model.update_reaction_bounds("CK","keep",0)
#recon_flux_model.update_reaction_bounds("biomass_reaction",0 ,200)
# "HGNC:29570 or HGNC:4331"
#recon_flux_model.rxn_dict["GLUN"]["grRule"] = "(27165.1) or (2744.1)"
# "HGNC:16085"
#recon_flux_model.rxn_dict["SERtm"]["grRule"] = "(94081.1)"
#recon_flux_model.update_reaction_bounds("biomass_reaction",2-tol ,2+tol)
#recon_flux_model.reaction_info("biomass_reaction")

#recon_flux_model.save_model_as_pkl(output_file_location)
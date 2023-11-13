import copy
import os
import time
import pandas as pd
import numpy as np
from pathlib import Path
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


positive_min_model_location = Path("Data/Models/pkl_models/recon_1b_t_cells/HR_ready_models")
# it should be possible to eventually just have HR ready models after exp_aligned
# This way breaks things up to save intermediates
null_space_S_path = Path("Data/null_space_S/recon_1b_t_cells")
positive_points_path = Path("Data/HR/HR_Points_analysis/recon_1b_t_cells/")
bidirectional_points_path = Path("Data/HR/HR_bi_directional_points/recon_1b_t_cells/")

positive_points = np.load(positive_points_path / os.listdir(positive_points_path)[0])[:, :-1]
print(np.shape(positive_points))
test_positive_point = copy.deepcopy(positive_points[0])
print(test_positive_point[4])
print(np.shape(test_positive_point))


# print(test_positive_point)
# input()


print(test_positive_point)

flux_model_list = []
for filename in os.listdir(positive_min_model_location):
	print(filename)
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_pkl_model(positive_min_model_location / filename)
	flux_model_list.append(recon_flux_model)

flux_model = flux_model_list[0]
flux_model_copy = copy.deepcopy(flux_model)
flux_model_copy.convert_model_to_bidirectional_flux()
lb, ub, S, b, rxn_list, met_list = flux_model_copy.dicts_to_mats()
print(np.shape(S))
rxn_array = np.array(rxn_list)
print(rxn_list)

print(np.size(rxn_array))

lb, ub, pos_S, b, pos_rxn_list, met_list = flux_model.dicts_to_mats()

print(np.shape(pos_S))
input()

dtf = pd.DataFrame(copy.deepcopy(positive_points[:1, :]), columns=pos_rxn_list)
print(np.shape(positive_points[:1, :]))
timestart = time.time()
for rxn_name in dtf.columns:
	if "_inverted" in rxn_name:
		if rxn_name.replace("_inverted", "") in dtf.columns:
			forward = rxn_name.replace("_inverted", "")
			dtf.loc[:, forward] = dtf.loc[:, forward].subtract(dtf.loc[:, rxn_name])
			dtf.pop(rxn_name)
		else:
			forward = rxn_name.replace("_inverted", "")
			dtf = dtf.rename(columns={rxn_name: forward})
			dtf.loc[:, forward] = -dtf.loc[:, forward]
print(timestart - time.time())
# print(dtf.to_numpy()[0])
# print(np.shape(dtf.to_numpy()[0]))
# print(test_positive_point[4])
simp_neg_ind, comp_neg_ind, comp_perm, cutter, exchange_cutter = flux_model.precomp_positive_flux_point_to_bi_dir(pos_rxn_list,prune_specific=["biomass_reaction"])
test_positive_point_to_test = flux_model.positive_flux_point_to_bi_dir(test_positive_point,simp_neg_ind, comp_neg_ind, comp_perm, cutter,exchange_cutter = exchange_cutter)
test_bi_S = flux_model.positive_S_to_bi_dir(pos_S,simp_neg_ind,cutter,exchange_cutter=exchange_cutter)
input()
pos_S

print("start")
print(timestart - time.time())
print("done")
print(list(test_positive_point_to_test))
print(list(dtf.to_numpy()[0]))
print(np.alltrue(dtf.to_numpy()[0] == test_positive_point_to_test))
input()
print(comp_neg_ind)
print(np.shape(test_positive_point))
input()

input()
comp_neg_ind = np.where(
	np.array([("_inverted" in i and i.replace("_inverted", "") in pos_rxn_array) for i in pos_rxn_array]), -1 * ones,
	zeros)
comp_perm_array = np.nonzero(
	np.where(np.array([("_inverted" not in i and (i + "_inverted") in pos_rxn_array) for i in pos_rxn_array]), zeros,
	         ones))[0]
print(cutter)
print(pos_rxn_list)
place_h = pos_rxn_array[cutter]
for i in range(np.size(place_h)):
	place_h[i] = place_h[i].replace("_inverted", "")
print(np.alltrue(place_h == rxn_array))
input()
print(rxn_array[cutter])
print(np.shape(S))
print(np.shape(internal_S))
start = time.time()
print("start")
input()
NS = flux_model.generate_null_space(internal_S)
print(time.time() - start)
input()
print(flux_model.rxn_dict["maintenance_ATP"])
flux_model.pos_S()
internal_lb, internal_ub, internal_S, b, internal_rxn, met_list = flux_model.generate_mat_wo_exchange(
	prune_specific=["biomass_reaction"])
print(np.shape(internal_S))

np.save(null_space_S_path / "Nullspace", NS)


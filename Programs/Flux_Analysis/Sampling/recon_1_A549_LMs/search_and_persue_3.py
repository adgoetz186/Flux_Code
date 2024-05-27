import os
import sys
import scipy.optimize as so
import pickle as pkl
import pandas as pd
import numpy as np
import scipy.io as sio
from pathlib import Path
import time
import gurobipy as gp

# _____ Setting the CWD to be Flux_Code BEGIN _____
# Flux_Code path here
if not Path(os.getcwd()).parts[-1] == "Flux_Code":
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

from Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
# _____ Setting the CWD to be Flux_Code END _____


def sampling_w_error(LMs,RPS_names,reaction_objectives,model,warmup_points,Internal_nullspace_S,alpha = 0.05):
	lb, ub, S, b, rxn_list, met_list = recon_flux_model.dicts_to_mats()
	print(lb[rxn_list.index('biomass_reaction')], ub[rxn_list.index('biomass_reaction')])
	RPS = np.zeros_like(lb)
	for name_ind in range(len(RPS_names)):
		RPS[rxn_list.index(RPS_names[name_ind])] += LMs[name_ind]
	# might wanna add a test here to make sure names line up, right now just test for complete
	# 10000, 5000
	start = time.time()
	HR_samples = np.array([])
	# ,thermo_const={"prune_specific": ["biomass_reaction"], "NS_internal": Internal_nullspace_S}
	while np.size(HR_samples) == 0:
		HR_samples = model.HRSampler_gene_bias_lincomb_pinch(warmup_points, 2000, 1000, RPS, 4,
																		thermo_const={
																			"prune_specific": ["biomass_reaction"],
																			"NS_internal": Internal_nullspace_S})
	objectives = []
	new_LMs = []
	for name_ind in range(len(RPS_names)):
		objectives.append(np.abs(reaction_objectives[name_ind] - np.average(HR_samples[499::100,rxn_list.index(RPS_names[name_ind])]))/reaction_objectives[name_ind])
		LM_error = reaction_objectives[name_ind] - np.average(HR_samples[499::100,rxn_list.index(RPS_names[name_ind])])
		new_LM = LMs[name_ind] - alpha*LM_error
		new_LMs.append(new_LM)

	return new_LMs,np.average(np.array(objectives)), HR_samples
	# HR_samples = flux_model.HRSampler_lincom(warmup_dict[flux_model_name], 10000, 500)

Internal_nullspace_S = sio.loadmat(Path("Data/null_space_S/recon_1_A549_4_30_2024/NS.mat"))['NS']
print(str(round(time.time()*1000)))
print(np.arange(2000)[499::100])

output_path = Path('Data/HR/HR_A549_therm_test_uniform')
op_LM = Path('Data/HR/optimal_LMs')
recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_fast_key_json_model(Path('Data/Models/json_models/fast_key_format/recon_1_A549_4_30_2024/HR_ready_model.json'))

# Used to get model pinched values
expt_recon_flux_model = Flux_Balance_Model()
expt_recon_flux_model.load_fast_key_json_model(Path('Data/Models/json_models/fast_key_format/recon_1_A549_4_30_2024/Experimentally_aligned.json'))










points = np.load('Data/HR/HR_A549_therm_test_uniform/tg7_initial_1715329936069.npy')
measured_exchange_rates_file = Path(f"Data/experimental_alignment_data/recon_1_A549/Measured_Rates.txt")
measured_exchange_rates = pd.read_csv(measured_exchange_rates_file)
exchange_names = measured_exchange_rates['reaction'].to_list()
print(exchange_names)
print('done')
measured_biomass_rate_file = Path(f"Data/experimental_alignment_data/recon_1_A549/Biomass_Data.csv")
measured_biomass_rates = pd.read_csv(measured_biomass_rate_file)
exchange_names.append(measured_biomass_rates['reaction'].to_list()[0])



lb, ub, S, b, rxn_list, met_list = expt_recon_flux_model.dicts_to_mats()
#list_of_reactions_to_opt = ['EX_lac_L(e)','EX_arg_L(e)_inverted']
list_of_reactions_to_opt = ['EX_lac_L(e)','EX_arg_L(e)_inverted','EX_asn_L(e)_inverted','EX_cys_L(e)_inverted','EX_glu_L(e)','EX_pro_L(e)']
dict_opt = {}
for rxn_name in list_of_reactions_to_opt:
	if "_inverted" not in rxn_name:
		dict_opt[rxn_name] = (lb[rxn_list.index(rxn_name)] + ub[rxn_list.index(rxn_name)])/2
	else:
		dict_opt[rxn_name] = -1*(lb[rxn_list.index(rxn_name.replace('_inverted',''))] + ub[rxn_list.index(rxn_name.replace('_inverted',''))]) / 2
print(dict_opt)



# For the hypergator us this gurobi env key
#recon_flux_model.add_gp_key_env_to_model('grb-ts.ufhpc')

warmup_points = np.load(Path('Data/HR/HR_Warmup/Recon_1_A549_4_30_2024/warmup_HR_ready_model.npy'))

reaction_objectives = list(dict_opt.values())
RPS_names = list(dict_opt.keys())
#sampling_w_error(RPS_names,)
all_LMs = np.array([0,0,0,0,0,0])
new_LMs,error_term,samples = sampling_w_error(all_LMs,RPS_names,reaction_objectives,recon_flux_model,warmup_points,Internal_nullspace_S,alpha = 50)
all_LMs = np.vstack((all_LMs,new_LMs))
error_terms = [error_term]
while error_term > 0.05 and len(error_terms) <= 25:
	new_LMs, error_term,samples = sampling_w_error(new_LMs, RPS_names, reaction_objectives, recon_flux_model,
										   warmup_points, Internal_nullspace_S, alpha=(50/(len(error_terms)+1)**1.1))
	all_LMs = np.vstack((all_LMs, new_LMs))
	print(all_LMs)
	print(error_term)
	error_terms.append(error_term)
	time.sleep(20)
np.save(output_path / (f'opt3_{str(round(time.time()*1000))}'), samples)
LM_ID = f'opt3_{str(round(time.time()*1000))}'
np.save(op_LM / LM_ID, all_LMs)
np.save(op_LM / ("Perf3_"+LM_ID), np.array(error_terms))

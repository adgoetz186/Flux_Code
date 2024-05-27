import os
import sys
import pandas as pd
import numpy as np
import scipy.io as sio
from pathlib import Path
import matplotlib.pyplot as plt
import time
import scipy.stats as st
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

recon_flux_model = Flux_Balance_Model()
recon_flux_model.load_fast_key_json_model(Path('Data/Models/json_models/fast_key_format/recon_1_A549_4_30_2024/HR_ready_model.json'))

expt_recon_flux_model = Flux_Balance_Model()
expt_recon_flux_model.load_fast_key_json_model(Path('Data/Models/json_models/fast_key_format/recon_1_A549_4_30_2024/Experimentally_aligned.json'))

#points = np.load('Data/HR/HR_A549_therm_test_uniform/tg3_sparse_initial_1715685471016.npy')
points = np.zeros(0)
point_folder = Path('Data/HR/HR_A549_therm_test_lactate_biased')
for point_array_loc in os.listdir(point_folder):
	if np.size(points) == 0:
		points = np.load(point_folder/point_array_loc)[500::100,:]
		print(np.shape(points))
	else:
		#continue
		points = np.vstack((points,np.load(point_folder/point_array_loc)[::100,:]))
		print(np.shape(points))
measured_exchange_rates_file = Path(f"Data/experimental_alignment_data/recon_1_A549/Measured_Rates.txt")
measured_exchange_rates = pd.read_csv(measured_exchange_rates_file)
exchange_names = measured_exchange_rates['reaction'].to_list()
exchange_values = measured_exchange_rates[' value'].to_list()
print(exchange_names)
print('done')
measured_biomass_rate_file = Path(f"Data/experimental_alignment_data/recon_1_A549/Biomass_Data.csv")
measured_biomass_rates = pd.read_csv(measured_biomass_rate_file)
exchange_names.append(measured_biomass_rates['reaction'].to_list()[0])
exchange_values.append(measured_biomass_rates['value'].to_list()[0])
print(measured_exchange_rates)
print(exchange_names)
print(exchange_values)

lb, ub, S, b, rxn_list, met_list = recon_flux_model.dicts_to_mats()

rxn_list.append('total_flux')
point_dtf = pd.DataFrame(data=points,columns=rxn_list)
for col_name in point_dtf.columns:
	if '_inverted' in col_name:
		if col_name.replace('_inverted','') in point_dtf.columns:

			point_dtf.loc[:,col_name.replace('_inverted','')] -= point_dtf.loc[:,col_name]
			point_dtf.drop(col_name,axis=1,inplace=True)

		else:

			point_dtf.loc[:, col_name] *= -1
			point_dtf.rename({col_name:col_name.replace('_inverted','')},inplace=True,axis=1)

print(point_dtf.columns.to_list())

print(len(rxn_list))
print(point_dtf)

print(rxn_list)
print(np.shape(points))
recon_flux_model.convert_model_to_bidirectional_flux()
lb, ub, S, b, rxn_list, met_list = recon_flux_model.dicts_to_mats()
e_lb, e_ub, e_S, e_b, e_rxn_list, e_met_list = expt_recon_flux_model.dicts_to_mats()
for rxn_ind in range(len(exchange_names)):
	points_to_use = point_dtf.loc[:,exchange_names[rxn_ind]].to_numpy()
	print(points_to_use)
	plt.scatter(np.arange(np.size(points_to_use)),points_to_use)
	#plt.title(exchange_names[rxn_ind])

	plt.hlines(lb[rxn_list.index(exchange_names[rxn_ind])],0,np.size(points_to_use))
	plt.hlines(e_ub[e_rxn_list.index(exchange_names[rxn_ind])], 0, np.size(points_to_use),color = 'red')
	print(st.percentileofscore(points_to_use,e_ub[e_rxn_list.index(exchange_names[rxn_ind])]))
	plt.hlines(ub[rxn_list.index(exchange_names[rxn_ind])], 0, np.size(points_to_use))
	plt.title(f"{exchange_names[rxn_ind]}\nmean is {np.round(st.percentileofscore(points_to_use,e_ub[e_rxn_list.index(exchange_names[rxn_ind])]),2)} percentile")
	plt.show()
	print(points)


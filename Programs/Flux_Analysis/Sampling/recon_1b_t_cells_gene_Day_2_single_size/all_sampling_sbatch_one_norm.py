import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.io as sio
from pathlib import Path
import time
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

model_location_B6_3 = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/HR_Ready_Day_2_single_size_all/B6-3.json")
warmup_sample_location_B6_3 = Path("Data/HR/HR_Warmup/recon_1b_t_cells_gene_Day_2_single_size_all/warmup_B6-3.npy")
gene_location_B6_3 = Path("Data/scRNA/Bulk/B63-D2-GEX-B7")

#model_location_TC_5 = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/HR_Ready_Day_2_single_size/TC-5.json")
#warmup_sample_location_TC_5 = Path("Data/HR/HR_Warmup/recon_1b_t_cells_gene_Day_2_single_size/warmup_TC-5.npy")
#gene_location_TC_5 = Path("Data/scRNA/Bulk/TC5-D2-GEX-B1")

Internal_nullspace_S = sio.loadmat(Path("Data/null_space_S/recon_1b_t_cells_single_size/Null_S_internal.mat"))['Null_S_internal']
trans_NS = np.transpose(Internal_nullspace_S)

run_type_arg_dict = {}
#run_type_arg_dict["B6_1_uniform"] = {"repeat":100,"model_location":model_location_B6_1,"gene_RNA_script":gene_location_B6_1,"warmup_sample_location":warmup_sample_location_B6_1,"thermo_const":None,"gene_penalty_mod":1e-12,"points_to_save":1000,"points_to_skip":100,"output_folder_path":Path("Data/HR/recon_1b_t_cells_gene_Day_2_single_B6_1")}
run_type_arg_dict["B6_3_gene_bias"] = {"repeat":5,"model_location":model_location_B6_3,"gene_RNA_script":gene_location_B6_3,"warmup_sample_location":warmup_sample_location_B6_3,"thermo_const":None,"gene_penalty_mod":1e-12,"points_to_save":10000,"points_to_skip":5000,"output_folder_path":Path("Data/HR/recon_1b_t_cells_gene_Day_2_single_B6_3_gene")}

#run_type_arg_dict["TC_5_uniform"] = {"repeat":100,"model_location":model_location_TC_5,"gene_RNA_script":gene_location_TC_5,"warmup_sample_location":warmup_sample_location_TC_5,"thermo_const":None,"gene_penalty_mod":1e-12,"points_to_save":1000,"points_to_skip":100,"output_folder_path":Path("Data/HR/recon_1b_t_cells_gene_Day_2_single_TC_5")}
#run_type_arg_dict["TC_5_gene_bias"] = {"repeat":5,"model_location":model_location_TC_5,"gene_RNA_script":gene_location_TC_5,"warmup_sample_location":warmup_sample_location_TC_5,"thermo_const":{"prune_specific": ["biomass_reaction"], "NS_internal": Internal_nullspace_S},"gene_penalty_mod":1e-12,"points_to_save":100,"points_to_skip":10,"output_folder_path":Path("Data/HR/recon_1b_t_cells_gene_Day_2_single_TC_5_gene")}
#run_type_arg_dict["TC_5_gene_bias"] = {"repeat":5,"model_location":model_location_TC_5,"gene_RNA_script":gene_location_TC_5,"warmup_sample_location":warmup_sample_location_TC_5,"thermo_const":None,"gene_penalty_mod":1e-12,"points_to_save":1000,"points_to_skip":1000,"output_folder_path":Path("Data/HR/recon_1b_t_cells_gene_Day_2_single_TC_5_gene")}


run_list = []
for run_name in run_type_arg_dict.keys():
	run_list += [run_name for i in range(run_type_arg_dict[run_name]["repeat"])]
print(run_list)
print(len(run_list))

for value in range(1):
	#run_for_core = run_list[]
	#value = int(sys.argv[1])
	#value = 3
	run_for_core = run_list[value]
	
	
	recon_flux_model = Flux_Balance_Model()
	recon_flux_model.load_fast_key_json_model(run_type_arg_dict[run_for_core]["model_location"])
	print(recon_flux_model.model_dict)
	print(recon_flux_model.model_dict.keys())
	#print(recon_flux_model.model_dict['gurobi_token'])
	#print(recon_flux_model.model_dict)
	#input()
	
	warmup_points = np.load(run_type_arg_dict[run_for_core]["warmup_sample_location"])
	lb_full, ub_full, S_full, b, rxn_list, met_list = recon_flux_model.dicts_to_mats()
	simp_neg_ind, comp_neg_ind, comp_perm, cutter, exchange_cutter = recon_flux_model.precomp_positive_flux_point_to_bi_dir(
		rxn_list, prune_specific=["biomass_reaction"])
	count = 0
	for i in warmup_points:
		ent_flux_point = recon_flux_model.positive_flux_point_to_bi_dir(i[:-1], simp_neg_ind,
		                                                    comp_neg_ind, comp_perm,
		                                                    cutter,
		                                                    exchange_cutter=exchange_cutter)
		therm_feas = recon_flux_model.generate_therm_model_new(trans_NS, ent_flux_point)
		if therm_feas == 2:
			print("good")
			count+=1
			print(count)
		elif therm_feas != 3:
			print(recon_flux_model.generate_therm_model_new(trans_NS, ent_flux_point))
	print(":(")
	print(warmup_points[:,-1])
	print(np.min(warmup_points[:,-1]),np.max(warmup_points[:,-1]))
	#plt.hist(warmup_points[:,-1])
	#plt.show()
	time_list = []
	
	gene_array = np.load(run_type_arg_dict[run_for_core]["gene_RNA_script"] / "matrix.npy")
	feature_array = np.load(run_type_arg_dict[run_for_core]["gene_RNA_script"] /"features.npy")
	gene_array = gene_array/np.average(gene_array)
	start = time.time()
	RAS = np.array(recon_flux_model.create_RAS_values(gene_array,feature_array,False))
	RPS = recon_flux_model.RAS_to_RPS_mat(RAS)
	# might wanna add a test here to make sure names line up, right now just test for complete
	# 100000, 500
	#print((np.max(warmup_points[:,-1]) - np.min(warmup_points[:,-1]))/np.min(warmup_points[:,-1]))
	#input()
	HR_samples = recon_flux_model.HRSampler_gene_bias_lincomb_pinch(warmup_points, run_type_arg_dict[run_for_core]["points_to_save"], run_type_arg_dict[run_for_core]["points_to_skip"],RPS,1.01,thermo_const=run_type_arg_dict[run_for_core]["thermo_const"],gene_penalty_mod=run_type_arg_dict[run_for_core]["gene_penalty_mod"])
	print(np.shape(HR_samples))
	print(HR_samples)
	lb_full, ub_full, S_full, b, rxn_list, met_list = recon_flux_model.dicts_to_mats()
	simp_neg_ind, comp_neg_ind, comp_perm, cutter, exchange_cutter = recon_flux_model.precomp_positive_flux_point_to_bi_dir(
		rxn_list, prune_specific=["biomass_reaction"])
	count = 0
	for i in HR_samples:
		ent_flux_point = recon_flux_model.positive_flux_point_to_bi_dir(i[:-1], simp_neg_ind,
		                                                    comp_neg_ind, comp_perm,
		                                                    cutter,
		                                                    exchange_cutter=exchange_cutter)
		therm_feas = recon_flux_model.generate_therm_model_new(trans_NS, ent_flux_point)
		if therm_feas == 2:
			print("good")
			count+=1
			print(count)
		elif therm_feas != 3:
			print(recon_flux_model.generate_therm_model_new(trans_NS, ent_flux_point))
	print(":(")
	input()
	#HR_samples = flux_model.HRSampler_lincom(warmup_dict[flux_model_name], 10000, 500)
	#np.save(run_type_arg_dict[run_for_core]["output_folder_path"] / (run_for_core + f"_HR_points_{run_list[:value].count(run_list[value])}"), HR_samples)
	print(f"The time for this run was: {(time.time() - start)/60} minutes")
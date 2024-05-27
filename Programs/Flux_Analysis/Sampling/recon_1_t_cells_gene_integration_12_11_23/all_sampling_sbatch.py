import os
import sys
import pandas as pd
import numpy as np
import scipy.io as sio
from pathlib import Path
import time
import gurobipy as gp

# Performs polytope sampling, see readme.txt for suggestions for running

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

model_location_B6_1 = Path("Data/Models/json_models/fast_key_format/recon_1_t_cells_12_11_23/HR_Ready_Day_2_single_size/B6-1.json")
warmup_sample_location_B6_1 = Path("Data/HR/HR_Warmup/recon_1_t_cells_gene_Day_2_single_size_12_11_23/warmup_B6-1.npy")
gene_location_B6_1 = Path("Data/scRNA/Bulk/B61-D2-GEX-B5")

model_location_B6_2 = Path("Data/Models/json_models/fast_key_format/recon_1_t_cells_12_11_23/HR_Ready_Day_2_single_size/B6-2.json")
warmup_sample_location_B6_2 = Path("Data/HR/HR_Warmup/recon_1_t_cells_gene_Day_2_single_size_12_11_23/warmup_B6-2.npy")
gene_location_B6_2 = Path("Data/scRNA/Bulk/B62-D2-GEX-B6")

model_location_B6_3 = Path("Data/Models/json_models/fast_key_format/recon_1_t_cells_12_11_23/HR_Ready_Day_2_single_size/B6-3.json")
warmup_sample_location_B6_3 = Path("Data/HR/HR_Warmup/recon_1_t_cells_gene_Day_2_single_size_12_11_23/warmup_B6-3.npy")
gene_location_B6_3 = Path("Data/scRNA/Bulk/B63-D2-GEX-B7")

model_location_B6_4 = Path("Data/Models/json_models/fast_key_format/recon_1_t_cells_12_11_23/HR_Ready_Day_2_single_size/B6-4.json")
warmup_sample_location_B6_4 = Path("Data/HR/HR_Warmup/recon_1_t_cells_gene_Day_2_single_size_12_11_23/warmup_B6-4.npy")
gene_location_B6_4 = Path("Data/scRNA/Bulk/B64-D2-GEX-B8")

model_location_TC_5 = Path("Data/Models/json_models/fast_key_format/recon_1_t_cells_12_11_23/HR_Ready_Day_2_single_size/TC-5.json")
warmup_sample_location_TC_5 = Path("Data/HR/HR_Warmup/recon_1_t_cells_gene_Day_2_single_size_12_11_23/warmup_TC-5.npy")
gene_location_TC_5 = Path("Data/scRNA/Bulk/TC5-D2-GEX-B1")

model_location_TC_6 = Path("Data/Models/json_models/fast_key_format/recon_1_t_cells_12_11_23/HR_Ready_Day_2_single_size/TC-6.json")
warmup_sample_location_TC_6 = Path("Data/HR/HR_Warmup/recon_1_t_cells_gene_Day_2_single_size_12_11_23/warmup_TC-6.npy")
gene_location_TC_6 = Path("Data/scRNA/Bulk/TC6-D2-GEX-B2")

model_location_TC_7 = Path("Data/Models/json_models/fast_key_format/recon_1_t_cells_12_11_23/HR_Ready_Day_2_single_size/TC-7.json")
warmup_sample_location_TC_7 = Path("Data/HR/HR_Warmup/recon_1_t_cells_gene_Day_2_single_size_12_11_23/warmup_TC-7.npy")
gene_location_TC_7 = Path("Data/scRNA/Bulk/TC7-D2-GEX-B3")

model_location_TC_8 = Path("Data/Models/json_models/fast_key_format/recon_1_t_cells_12_11_23/HR_Ready_Day_2_single_size/TC-8.json")
warmup_sample_location_TC_8 = Path("Data/HR/HR_Warmup/recon_1_t_cells_gene_Day_2_single_size_12_11_23/warmup_TC-8.npy")
gene_location_TC_8 = Path("Data/scRNA/Bulk/TC8-D2-GEX-B4")

Internal_nullspace_S = sio.loadmat(Path("Data/null_space_S/recon_1b_t_cells_single_size/Null_S_internal.mat"))['Null_S_internal']

skip_uniform = 5000
keep_uniform = 10000

skip_gene_therm = 500
keep_gene_therm = 4100
run_type_arg_dict = {}
run_type_arg_dict["B6_1_uniform"] = {"repeat":1,"model_location":model_location_B6_1,"gene_RNA_script":gene_location_B6_1,"warmup_sample_location":warmup_sample_location_B6_1,"thermo_const":None,"gene_penalty_mod":1e-12,"points_to_save":keep_uniform,"points_to_skip":skip_uniform,"output_folder_path":Path("Data/HR/HR_fin_2/B6_1")}
run_type_arg_dict["B6_1_gene_bias_therm"] = {"repeat":40,"model_location":model_location_B6_1,"gene_RNA_script":gene_location_B6_1,"warmup_sample_location":warmup_sample_location_B6_1,"thermo_const":{"prune_specific": ["biomass_reaction"], "NS_internal": Internal_nullspace_S},"gene_penalty_mod":1e-0,"points_to_save":keep_gene_therm,"points_to_skip":skip_gene_therm,"output_folder_path":Path("Data/HR/HR_fin_2/B6_1_gene_therm")}
run_type_arg_dict["B6_2_uniform"] = {"repeat":1,"model_location":model_location_B6_2,"gene_RNA_script":gene_location_B6_2,"warmup_sample_location":warmup_sample_location_B6_2,"thermo_const":None,"gene_penalty_mod":1e-12,"points_to_save":keep_uniform,"points_to_skip":skip_uniform,"output_folder_path":Path("Data/HR/HR_fin_2/B6_2")}
run_type_arg_dict["B6_2_gene_bias_therm"] = {"repeat":40,"model_location":model_location_B6_2,"gene_RNA_script":gene_location_B6_2,"warmup_sample_location":warmup_sample_location_B6_2,"thermo_const":{"prune_specific": ["biomass_reaction"], "NS_internal": Internal_nullspace_S},"gene_penalty_mod":1e-0,"points_to_save":keep_gene_therm,"points_to_skip":skip_gene_therm,"output_folder_path":Path("Data/HR/HR_fin_2/B6_2_gene_therm")}
run_type_arg_dict["B6_3_uniform"] = {"repeat":1,"model_location":model_location_B6_3,"gene_RNA_script":gene_location_B6_3,"warmup_sample_location":warmup_sample_location_B6_3,"thermo_const":None,"gene_penalty_mod":1e-12,"points_to_save":keep_uniform,"points_to_skip":skip_uniform,"output_folder_path":Path("Data/HR/HR_fin_2/B6_3")}
run_type_arg_dict["B6_3_gene_bias_therm"] = {"repeat":40,"model_location":model_location_B6_3,"gene_RNA_script":gene_location_B6_3,"warmup_sample_location":warmup_sample_location_B6_3,"thermo_const":{"prune_specific": ["biomass_reaction"], "NS_internal": Internal_nullspace_S},"gene_penalty_mod":1e-0,"points_to_save":keep_gene_therm,"points_to_skip":skip_gene_therm,"output_folder_path":Path("Data/HR/HR_fin_2/B6_3_gene_therm")}
run_type_arg_dict["B6_4_uniform"] = {"repeat":1,"model_location":model_location_B6_4,"gene_RNA_script":gene_location_B6_4,"warmup_sample_location":warmup_sample_location_B6_4,"thermo_const":None,"gene_penalty_mod":1e-12,"points_to_save":keep_uniform,"points_to_skip":skip_uniform,"output_folder_path":Path("Data/HR/HR_fin_2/B6_4")}
run_type_arg_dict["B6_4_gene_bias_therm"] = {"repeat":40,"model_location":model_location_B6_4,"gene_RNA_script":gene_location_B6_4,"warmup_sample_location":warmup_sample_location_B6_4,"thermo_const":{"prune_specific": ["biomass_reaction"], "NS_internal": Internal_nullspace_S},"gene_penalty_mod":1e-0,"points_to_save":keep_gene_therm,"points_to_skip":skip_gene_therm,"output_folder_path":Path("Data/HR/HR_fin_2/B6_4_gene_therm")}
run_type_arg_dict["TC_5_uniform"] = {"repeat":1,"model_location":model_location_TC_5,"gene_RNA_script":gene_location_TC_5,"warmup_sample_location":warmup_sample_location_TC_5,"thermo_const":None,"gene_penalty_mod":1e-12,"points_to_save":keep_uniform,"points_to_skip":skip_uniform,"output_folder_path":Path("Data/HR/HR_fin_2/TC_5")}
run_type_arg_dict["TC_5_gene_bias_therm"] = {"repeat":40,"model_location":model_location_TC_5,"gene_RNA_script":gene_location_TC_5,"warmup_sample_location":warmup_sample_location_TC_5,"thermo_const":{"prune_specific": ["biomass_reaction"], "NS_internal": Internal_nullspace_S},"gene_penalty_mod":1e-0,"points_to_save":keep_gene_therm,"points_to_skip":skip_gene_therm,"output_folder_path":Path("Data/HR/HR_fin_2/TC_5_gene_therm")}
run_type_arg_dict["TC_6_uniform"] = {"repeat":1,"model_location":model_location_TC_6,"gene_RNA_script":gene_location_TC_6,"warmup_sample_location":warmup_sample_location_TC_6,"thermo_const":None,"gene_penalty_mod":1e-12,"points_to_save":keep_uniform,"points_to_skip":skip_uniform,"output_folder_path":Path("Data/HR/HR_fin_2/TC_6")}
run_type_arg_dict["TC_6_gene_bias_therm"] = {"repeat":40,"model_location":model_location_TC_6,"gene_RNA_script":gene_location_TC_6,"warmup_sample_location":warmup_sample_location_TC_6,"thermo_const":{"prune_specific": ["biomass_reaction"], "NS_internal": Internal_nullspace_S},"gene_penalty_mod":1e-0,"points_to_save":keep_gene_therm,"points_to_skip":skip_gene_therm,"output_folder_path":Path("Data/HR/HR_fin_2/TC_6_gene_therm")}
run_type_arg_dict["TC_7_uniform"] = {"repeat":1,"model_location":model_location_TC_7,"gene_RNA_script":gene_location_TC_7,"warmup_sample_location":warmup_sample_location_TC_7,"thermo_const":None,"gene_penalty_mod":1e-12,"points_to_save":keep_uniform,"points_to_skip":skip_uniform,"output_folder_path":Path("Data/HR/HR_fin_2/TC_7")}
run_type_arg_dict["TC_7_gene_bias_therm"] = {"repeat":40,"model_location":model_location_TC_7,"gene_RNA_script":gene_location_TC_7,"warmup_sample_location":warmup_sample_location_TC_7,"thermo_const":{"prune_specific": ["biomass_reaction"], "NS_internal": Internal_nullspace_S},"gene_penalty_mod":1e-0,"points_to_save":keep_gene_therm,"points_to_skip":skip_gene_therm,"output_folder_path":Path("Data/HR/HR_fin_2/TC_7_gene_therm")}
run_type_arg_dict["TC_8_uniform"] = {"repeat":1,"model_location":model_location_TC_8,"gene_RNA_script":gene_location_TC_8,"warmup_sample_location":warmup_sample_location_TC_8,"thermo_const":None,"gene_penalty_mod":1e-12,"points_to_save":keep_uniform,"points_to_skip":skip_uniform,"output_folder_path":Path("Data/HR/HR_fin_2/TC_8")}
run_type_arg_dict["TC_8_gene_bias_therm"] = {"repeat":40,"model_location":model_location_TC_8,"gene_RNA_script":gene_location_TC_8,"warmup_sample_location":warmup_sample_location_TC_8,"thermo_const":{"prune_specific": ["biomass_reaction"], "NS_internal": Internal_nullspace_S},"gene_penalty_mod":1e-0,"points_to_save":keep_gene_therm,"points_to_skip":skip_gene_therm,"output_folder_path":Path("Data/HR/HR_fin_2/TC_8_gene_therm")}

run_list = []
for run_name in run_type_arg_dict.keys():
	run_list += [run_name for i in range(run_type_arg_dict[run_name]["repeat"])]


value = int(1)
print(run_list[value])

run_for_core = run_list[value]

recon_flux_model = Flux_Balance_Model()

recon_flux_model.load_fast_key_json_model(run_type_arg_dict[run_for_core]["model_location"])

# For the hypergator us this gurobi env key
#recon_flux_model.add_gp_key_env_to_model('grb-ts.ufhpc')

warmup_points = np.load(run_type_arg_dict[run_for_core]["warmup_sample_location"])

time_list = []

gene_array = np.load(run_type_arg_dict[run_for_core]["gene_RNA_script"] / "matrix.npy")

feature_array = np.load(run_type_arg_dict[run_for_core]["gene_RNA_script"] /"features.npy")

gene_array = gene_array/np.average(gene_array)
start = time.time()
RAS = np.array(recon_flux_model.create_RAS_values(gene_array,feature_array,False))
RPS = recon_flux_model.RAS_to_RPS_mat(RAS)

# might wanna add a test here to make sure names line up, right now just test for complete
# 100000, 500
HR_samples = np.array([])
while np.size(HR_samples) == 0:
	HR_samples = recon_flux_model.HRSampler_gene_bias_lincomb_pinch(warmup_points, run_type_arg_dict[run_for_core]["points_to_save"], run_type_arg_dict[run_for_core]["points_to_skip"],RPS,1.01,thermo_const=run_type_arg_dict[run_for_core]["thermo_const"],gene_penalty_mod=run_type_arg_dict[run_for_core]["gene_penalty_mod"])
#HR_samples = flux_model.HRSampler_lincom(warmup_dict[flux_model_name], 10000, 500)
print(HR_samples)
if not os.path.isdir(run_type_arg_dict[run_for_core]["output_folder_path"]):
	os.mkdir(run_type_arg_dict[run_for_core]["output_folder_path"])
np.save(run_type_arg_dict[run_for_core]["output_folder_path"] / (run_for_core + f"_HR_points_{run_list[:value].count(run_list[value])}"), HR_samples)
print(f"The time for this run, {run_for_core}, was: {(time.time() - start)/60} minutes")
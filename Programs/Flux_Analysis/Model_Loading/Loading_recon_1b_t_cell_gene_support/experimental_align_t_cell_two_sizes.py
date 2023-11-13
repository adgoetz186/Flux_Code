import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import scipy.io as sio
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

model_file_location = Path("Data/Models/json_models/fast_key_format/recon_1b_t_cells/Complete_prior.json")
gene_data_file_location = Path("Data/scRNA/Bulk/")
cfp = Path("Data/experimental_alignment_data/recon_1b_t_cells_gene_support/Individual_Mice")
# Sets the day to read the scRNA file
Day = 2
ofp = Path(f"Data/Models/json_models/fast_key_format/recon_1b_t_cells/exp_aligned_Day_{Day}_two_sizes")

# Dictionary of mice names
mouse_number_filename_dict = {1: "B6-1", 2: "B6-2", 3: "B6-3", 4: "B6-4", 5: "TC-5", 6: "TC-6", 7: "TC-7", 8: "TC-8"}

# specifies paths for:
# allowed exchange reactions, measured rates, biomass datasets, experimental alignment results, experimental
# alignment results, and cell sizes for each biological replicate using the keys given above
expt_allow_ex_files_dict = {i: cfp / mouse_number_filename_dict[i] / "Allowed_Exchange_Reactions" for i in
                            mouse_number_filename_dict.keys()}
expt_measured_rates_dict = {i: cfp / mouse_number_filename_dict[i] / "Measured_Rates.csv" for i in
                            mouse_number_filename_dict.keys()}
biomass_rates_dict = {i: cfp / mouse_number_filename_dict[i] / "Biomass_Data.csv" for i in
                      mouse_number_filename_dict.keys()}
output_opt_data_file_location_a_dict = {
	i: cfp / mouse_number_filename_dict[i] / f"experimental_alignment_result_data_a_Day_{Day}" for i in
	mouse_number_filename_dict.keys()}
output_opt_data_file_location_b_dict = {
	i: cfp / mouse_number_filename_dict[i] / f"experimental_alignment_result_data_b_Day_{Day}_two_sizes" for i in
	mouse_number_filename_dict.keys()}
cell_size_file_dict = {i: cfp / mouse_number_filename_dict[i] / f"cell_size.txt" for i in
                       mouse_number_filename_dict.keys()}

list_of_mrna = []
list_of_median_cell_mrna = []
list_of_cell_count = []
weighted_average_denom = 0
# for folder in os.listdir(gene_data_file_location):
#    if "D" + str(Day) in folder:
#        print(folder, Day)
# Loads list of gene names in ensemble Stable IDs (ENSMUSG###########) format
#        features = np.loadtxt(gene_data_file_location / folder / "features.tsv.gz", delimiter="\t", dtype=str)
#        feature_list = list(features[:, 0])

# Loads gene matrix
#        gene_mat = np.array(sio.mmread(gene_data_file_location / folder / "matrix.mtx.gz").todense())

# counts the number of cells in the scRNAseq matrix for calculating the bulk scRNAseq matrix
#        weighted_average_denom += np.shape(gene_mat)[1]

# saves the average of the scRNAseq matrix
#        list_of_mrna.append(np.average(gene_mat))

# finds the expression levels of each cell in the scRNAseq matrix then obtains the median expression level
#        list_of_median_cell_mrna.append(np.median(np.average(gene_mat,axis=0)))

#        list_of_cell_count.append(np.shape(gene_mat)[1])
#        print(np.shape(gene_mat))
#        print(gene_mat)
#        print(list_of_mrna)
#        print(weighted_average_denom)
# plt.scatter(list_of_mrna,list_of_cell_count)
# plt.show()
# plt.scatter(list_of_median_cell_mrna,list_of_cell_count)
# plt.show()

# print(list_of_mrna)
# input()
# bulk_gene_mat = np.sum(np.hstack(list_of_mrna), axis=1) / weighted_average_denom

error_array = ""
for mouse_number in mouse_number_filename_dict.keys():
	output_opt_data_file_location = output_opt_data_file_location_a_dict[mouse_number]
	
	
	
	with open(output_opt_data_file_location.parent / (output_opt_data_file_location.stem + ".pkl"), "rb") as readfile:
		model_dict = pickle.load(readfile)
	X_array = model_dict["model_alignment"]
	exp = model_dict["experimental_fluxes"]
	ind = model_dict["experimental_flux_index"]
	alpha_list = model_dict["alpha_list"]

	if error_array == "":
		error_array = np.zeros((len(mouse_number_filename_dict), len(alpha_list)))
	
	labels = [model_dict["saved_rxn_name_order"][i] for i in ind]
	
	# remove valine flux from consideration in error term
	del ind[labels.index('EX_val_L(e)')]
	
	for alpha_index in range(np.size(alpha_list)):
		alpha = alpha_list[alpha_index]
		expt = np.array([exp[i] for i in ind])
		mdl = np.array([X_array[alpha_index][i] for i in ind])
		error_array[mouse_number - 1, alpha_index] = np.sum((mdl - expt * alpha) ** 2 / (expt * alpha) ** 2)
print(error_array)
error_array_B6 = error_array[:4]
error_vec_B6 = np.average(error_array_B6, axis=0)
min_ind_B6 = np.argmin(error_vec_B6)
error_array_TC = error_array[4:]
error_vec_TC = np.average(error_array_TC, axis=0)
min_ind_TC = np.argmin(error_vec_TC)
print(error_array_B6)
print(error_array_TC)
print(error_vec_B6)
print(error_vec_TC)

mass_of_cell  = 10 ** (-12) * (1 / np.array(alpha_list))
plt.scatter(mass_of_cell[min_ind_B6] * 10 ** 12, error_vec_B6[min_ind_B6])
plt.scatter(mass_of_cell[min_ind_TC] * 10 ** 12, error_vec_TC[min_ind_TC])
plt.plot(mass_of_cell * 10 ** 12, error_vec_B6)
plt.plot(mass_of_cell * 10 ** 12, error_vec_TC)
plt.legend()
plt.show()
# for i in error_array:
for i in mouse_number_filename_dict.keys():
	plt.plot(mass_of_cell * 10 ** 12, error_array[i - 1], label=mouse_number_filename_dict[i])
plt.legend()
plt.scatter(mass_of_cell[min_ind_B6] * 10 ** 12, error_vec_B6[min_ind_B6])
plt.scatter(mass_of_cell[min_ind_TC] * 10 ** 12, error_vec_TC[min_ind_TC])
plt.show()
print("done")
mass_of_cell_B6 = mass_of_cell[min_ind_B6]
mass_of_cell_TC = mass_of_cell[min_ind_TC]
print(mass_of_cell_B6, mass_of_cell_TC)
for mouse_number in mouse_number_filename_dict.keys():
	experimental_Allowed_Exchange_file_location = expt_allow_ex_files_dict[mouse_number]
	experimental_measured_rates = expt_measured_rates_dict[mouse_number]
	biomass_rates = biomass_rates_dict[mouse_number]
	output_opt_data_file_location = output_opt_data_file_location_b_dict[mouse_number]
	cell_size_file = cell_size_file_dict[mouse_number]
	if mouse_number <= 4:
		mass_of_cell = mass_of_cell_B6
	else:
		mass_of_cell = mass_of_cell_TC
	print(f"mass of cell = {mass_of_cell}")

	alpha_list = [10 ** (-12) * (1 / mass_of_cell)]
	if error_array == "":
		error_array = np.zeros((len(mouse_number_filename_dict), len(alpha_list)))
	recon1_rxn_added = Flux_Balance_Model()
	recon1_rxn_added.load_fast_key_json_model(model_file_location)
	recon1_rxn_added.model_dict["model_name"] = mouse_number_filename_dict[mouse_number]
	recon1_rxn_added.pinch_restricted_exchange_reactions(experimental_Allowed_Exchange_file_location,
	                                                     restore_essential=False, true_exchange_signifier="(e)",
	                                                     restore_previously_pinched=True)
	tol = 1e-7
	br = pd.read_csv(biomass_rates)["value"].to_numpy()[0]
	recon1_rxn_added.update_reaction_bounds("biomass_reaction", br - tol, br + tol)
	# There did seem to be some risk of not allowed fluxes being just very close to 0
	# Pinching wont fix this
	# Write option to pinch to zero if zero falls between lb and ub
	# flux_modelb.pinch_reactions(10**(-5))
	# alpha_list = 10**np.linspace(20,22,30)
	
	# A549 has 244 pg of protein
	# the fraction of biomass protein is 0.706
	
	for folder in os.listdir(gene_data_file_location):
		if "D" + str(Day) == folder.split("-")[1] and mouse_number_filename_dict[mouse_number].replace("-", "") == \
				folder.split("-")[0]:
			features = np.load(gene_data_file_location / folder / "features.npy")
			bulk_gene_mat = np.load(gene_data_file_location / folder / "matrix.npy")
	bulk_gene_mat = bulk_gene_mat/np.average(bulk_gene_mat)

	RAS_dict = recon1_rxn_added.create_RAS_dict(bulk_gene_mat, features)

	for rxn_name in RAS_dict.keys():
		if RAS_dict[rxn_name] == 0:
			del recon1_rxn_added.model_dict["rxn_dict"][rxn_name]
	
	recon1_rxn_added.fit_to_experimental_data(experimental_measured_rates, alpha_list,
	                                          search_save_path=output_opt_data_file_location)
	
	if recon1_rxn_added.test_feasibility():
		recon1_rxn_added.save_model_as_fast_key_json(ofp)
	else:
		print("Model Infeasible")


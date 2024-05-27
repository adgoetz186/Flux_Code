import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Performs initial alignment of t-cell sizes with experimental data
# Each conditions error is found uniquely here allowing for one or more conditions to be grouped to find more than
# one cell size. At the time of writting this we only allow for one cell size.

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
print(os.getcwd())
model_file_location = Path("Data/Models/json_models/fast_key_format/recon_1_A549_4_30_2024/Complete_prior.json")

gene_features_location = Path("Data/scRNA/A549_Bulk/entrez_symbol_features.npy")
gene_matrix_location = Path("Data/scRNA/A549_Bulk/matrix.npy")
# Sets the day to read the scRNA file



# specifies paths for:
# allowed exchange reactions, measured rates, biomass datasets, experimental alignment results, experimental
# alignment results, and cell sizes for each biological replicate using the keys given above
experimental_Allowed_Exchange_file_location = Path("Data/experimental_alignment_data/recon_1_A549/Allowed_Exchange_Reactions")
experimental_measured_rates = Path("Data/experimental_alignment_data/recon_1_A549/Measured_Rates.txt")
biomass_rates = Path("Data/experimental_alignment_data/recon_1_A549/Biomass_Data.csv")
output_opt_data_file_location = Path("Data/experimental_alignment_data/recon_1_A549/align_result_data/experimental_alignment_result_data.pkl")
ofp = Path("Data/Models/json_models/fast_key_format/recon_1_A549_4_30_2024")
cell_size_file = Path("Data/experimental_alignment_data/recon_1_A549/cell_size.txt")


list_of_mrna = []
list_of_median_cell_mrna = []
list_of_cell_count = []
weighted_average_denom = 0
error_array = ""



cell_size_array = np.loadtxt(cell_size_file,delimiter=",")
mass_of_cell = 10**np.linspace(np.log10(cell_size_array[0]),np.log10(cell_size_array[1]),int(cell_size_array[2]))
print(f"mass of cell = {mass_of_cell}")
mass_of_cell *= 10**(-12)
alpha_list = list(10**(-12)*(1/mass_of_cell))
#if isinstance(error_array, str):
#    error_array = np.zeros((len(mouse_number_filename_dict),len(alpha_list)))
recon1_rxn_added = Flux_Balance_Model()
recon1_rxn_added.load_fast_key_json_model(model_file_location)
recon1_rxn_added.model_dict["model_name"] = "Experimentally_aligned"


print(recon1_rxn_added.test_feasibility())
recon1_rxn_added.pinch_restricted_exchange_reactions(experimental_Allowed_Exchange_file_location,
                                                     restore_essential=False, true_exchange_signifier="(e)",
                                                     restore_previously_pinched=True)


tol = 1e-7
br = pd.read_csv(biomass_rates)["value"].to_numpy()[0]
recon1_rxn_added.update_reaction_bounds("biomass_reaction", br - tol, br + tol)



entrez_features = np.load(gene_features_location)[:,0]
symbol_features = np.load(gene_features_location)[:,1]
#print(list(symbol_features).index('ADSS'))
#print(list(entrez_features)[3163])

#for i in range(np.shape(features)[0]):
#    print(features[i])


#features = np.reshape(features,(-1,1))
print(entrez_features)

print(np.size(entrez_features.flatten()) - np.size(np.unique(entrez_features.flatten())))
print(np.size(np.char.upper(entrez_features.flatten())) - np.size(np.unique(np.char.upper(entrez_features.flatten()))))
print(np.size(np.char.upper(entrez_features.flatten())))
print(np.size(np.unique(np.char.upper(entrez_features.flatten()))))

bulk_gene_mat = np.load(gene_matrix_location)

ulist = entrez_features.flatten()
upper_count = 0
upper_gene = 0
lower_count = 0
lower_gene = 0
for i in range(len(ulist)):
    if (ulist[i].isupper()):
        upper_count += 1
        upper_gene += bulk_gene_mat[i]
    else:
        lower_count +=1
        lower_gene += bulk_gene_mat[i]
print(upper_count, lower_count)
print(upper_gene, lower_gene)

print(entrez_features)


print(np.load("Data/scRNA/Bulk/B61-D0-GEX-A5/features.npy"))
print(np.load("Data/scRNA/Bulk/B61-D0-GEX-A5/matrix.npy"))
print(bulk_gene_mat)
print(bulk_gene_mat)
#bulk_gene_mat+=1

#print(bulk_gene_mat[list(symbol_features).index("ASPCTr")])
#print(recon1_rxn_added.model_dict['rxn_dixt']["ASPCTr"])
recon1_rxn_added.reaction_info('ADSS')
#print(entrez_features)

#RAS_dict = recon1_rxn_added.create_RAS_values(bulk_gene_mat,np.reshape(entrez_features,(-1,1)))

#count = 0
#for rxn_name in RAS_dict.keys():
#    count+=1
#    if RAS_dict[rxn_name] == 0 and rxn_name not in ['']:
#        del recon1_rxn_added.model_dict["rxn_dict"][rxn_name]
#    elif RAS_dict[rxn_name] == 0:
#        print(rxn_name)


# A549 has 244 pg of protein
# the fraction of biomass protein is 0.706
recon1_rxn_added.fit_to_experimental_data(experimental_measured_rates,alpha_list,search_save_path = output_opt_data_file_location)
print(recon1_rxn_added.test_feasibility())
if recon1_rxn_added.test_feasibility():
    recon1_rxn_added.save_model_as_fast_key_json(ofp)
else:
    print("Model Infeasible")


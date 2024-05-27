import os
from pathlib import Path
import numpy as np
import sys
import pandas as pd
import scipy.io as sio

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

gene_folder_loc = Path('Data/scRNA/Raw')
gene_expression_dict = {}
gene_expression_bin_dict = {}
metnames = np.array([i.replace("'","").replace(" ","") for i in np.loadtxt(Path('Data/RNA_extra_data/gene_names_metabolism.csv'),dtype=str,delimiter=",")])
for foldername in os.listdir(gene_folder_loc):
    names = np.loadtxt(gene_folder_loc / foldername / 'features.tsv.gz',dtype=str,delimiter='\t')[:,0]
    index_list = []
    for i in metnames:
        index_list.append(np.where(names==i))

    sorter = np.argsort(names)

    sorted = sorter[np.searchsorted(names, metnames, sorter=sorter)]

    seq_mat = sio.mmread(gene_folder_loc/foldername/'matrix.mtx.gz').todense()
    #seq_mat = np.random.randint(0,3,(32285, 6011))
    met_mat = seq_mat[sorted,:]


    binarymat = np.ceil(met_mat/np.max(met_mat))
    gene_expression_bin_dict[foldername] = np.average(binarymat)
    gene_expression_dict[foldername] = np.average(met_mat)
    print(gene_expression_dict)
    print(gene_expression_bin_dict)

for name in gene_expression_dict.keys():
    day = name.split("-")[1]
    mouse = name.split("-")[0]
    if day == "D0":
        for name_2 in gene_expression_dict.keys():
            day_2 = name_2.split("-")[1]
            mouse_2 = name_2.split("-")[0]
            if mouse == mouse_2 and day_2 == "D2":
                print(f"Mouse: {mouse}, day 2 - day 0 exp: {gene_expression_dict[name_2]-gene_expression_dict[name] }, day 0 avg: {gene_expression_dict[name]}, day 2 avg {gene_expression_dict[name_2]}")
                print(
                    f"Mouse: {mouse}, day 2 - day 0 binary exp: {gene_expression_bin_dict[name_2] - gene_expression_bin_dict[name]}, day 0 binary avg: {gene_expression_bin_dict[name]}, day 2 binary avg {gene_expression_bin_dict[name_2]}")


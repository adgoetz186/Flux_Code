import os

import pickle
import numpy as np

os.chdir("../../..")

model_scRNAseq = np.loadtxt("Data/Intermediate/scRNAseq_for_model/A549_Recon_2_2_hgnc.txt",delimiter=",",dtype=str)
gene_names = model_scRNAseq[:,0]
scRNAseq_mat = model_scRNAseq[:,1:].astype(float)

flux_model_2 = pickle.load(open("Data/Intermediate/Essential_Flux_Model/A549_recon2_2d.mdl", "rb"))

fva_results = flux_model_2.fva()
for i in range(len(flux_model_2.reaction_names)):
	print(flux_model_2.reaction_names[i],fva_results[i],flux_model_2.lb[i],flux_model_2.ub[i])
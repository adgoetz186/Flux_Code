import os
import pickle
import numpy as np

os.chdir("../../..")

model_scRNAseq = np.loadtxt("Data/Intermediate/scRNAseq_for_model/A549_Recon_2_2_hgnc.txt",delimiter=",",dtype=str)
print(model_scRNAseq)
gene_names = model_scRNAseq[:,0]
print(len(gene_names))
input(0)
scRNAseq_mat = model_scRNAseq[:,1:].astype(float)

flux_model_2 = pickle.load(open("Data/Intermediate/Essential_Flux_Model/A549_recon2_2d.mdl", "rb"))

print(gene_names)
RAS_matrix = flux_model_2.convert_scRNAseq_to_RAS_matrix(scRNAseq_mat,gene_names)
reactions = np.array(flux_model_2.reaction_names,dtype=str)
reactions = np.reshape(np.array(reactions),(-1,1))
array_to_save = np.hstack((reactions,RAS_matrix))
np.savetxt("Data/Intermediate/RAS/A549_Recon_2_2_hgnc.txt",array_to_save,delimiter=",",fmt="%s")


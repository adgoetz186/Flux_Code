import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.io as sio

os.chdir("../../")
print(np.loadtxt("Data/Input/Experimental_Data/sc_RNA_seq/B61-D0_GEX_A5/filtered_feature_bc_matrix/barcodes.tsv.gz",delimiter=" ",dtype=str))
print(np.loadtxt("Data/Input/Experimental_Data/sc_RNA_seq/B61-D0_GEX_A5/filtered_feature_bc_matrix/features.tsv.gz",delimiter="\t",dtype=str))

print(np.shape(np.loadtxt("Data/Input/Experimental_Data/sc_RNA_seq/B61-D0_GEX_A5/filtered_feature_bc_matrix/barcodes.tsv.gz",delimiter=" ",dtype=str)))
print(np.shape(np.loadtxt("Data/Input/Experimental_Data/sc_RNA_seq/B61-D0_GEX_A5/filtered_feature_bc_matrix/features.tsv.gz",delimiter="\t",dtype=str)))

print(np.shape(sio.mmread("Data/Input/Experimental_Data/sc_RNA_seq/B61-D0_GEX_A5/filtered_feature_bc_matrix/matrix.mtx.gz")))

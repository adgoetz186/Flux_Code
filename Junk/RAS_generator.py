import os
import pickle
import Flux_Code.Programs.Flux_Analysis.Model_Loading.Functions as lf

os.chdir("../")
scRNAseq_gene_list = list(lf.RNA_Seq_Load("Data/Input/scRNA-seq-data/A549_sciCAR_data.mat")["RNA"]["Features"])
scRNA_mat = lf.RNA_Seq_Load("Data/Input/scRNA-seq-data/A549_sciCAR_data.mat")["RNA"]["data"].todense()
flux_model_2 = pickle.load(open("Data/Intermediate/Essential_Flux_Model/A549_recon2_2d.mdl", "rb"))

print(scRNAseq_gene_list)
print(flux_model_2.convert_scRNAseq_to_RAS_matrix(scRNA_mat,scRNAseq_gene_list))
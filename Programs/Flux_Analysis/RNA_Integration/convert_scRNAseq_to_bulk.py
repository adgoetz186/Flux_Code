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

gene_data_file_location = Path("Data/scRNA/Raw/")
bulk_gene_data_file_location = Path("Data/scRNA/Bulk/")



for folder in os.listdir(gene_data_file_location):
	print(f"Starting: {folder}")
	features = np.loadtxt(gene_data_file_location / folder / "features.tsv.gz", delimiter="\t", dtype=str)
	gene_mat = np.array(sio.mmread(gene_data_file_location / folder / "matrix.mtx.gz").todense())
	bulk_gene_mat = np.average(gene_mat, axis=1)
	if not os.path.exists(bulk_gene_data_file_location / folder):
		os.mkdir(bulk_gene_data_file_location / folder)
	np.save(bulk_gene_data_file_location / folder / "matrix",bulk_gene_mat)
	np.save(bulk_gene_data_file_location / folder / "features",features)



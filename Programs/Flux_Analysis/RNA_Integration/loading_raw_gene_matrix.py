import os
import pandas as pd
import numpy as np
from pathlib import Path
import scipy.io as sio
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

Raw_gene_file_location = Path("Data/scRNA/Raw")
scRNA_npy_location = Path("Data/scRNA/Processed/npy_matrix")

barcodes = np.loadtxt(Raw_gene_file_location/"barcodes.tsv.gz",delimiter="\t",dtype = str)
features = np.loadtxt(Raw_gene_file_location/"features.tsv.gz",delimiter="\t",dtype = str)
matrix = sio.mmread(Raw_gene_file_location/"matrix.mtx.gz").todense()
np.save(scRNA_npy_location/"trial",matrix)



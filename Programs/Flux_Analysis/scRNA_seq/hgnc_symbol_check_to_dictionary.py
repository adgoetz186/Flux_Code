import os
import json
import Flux_Code.Programs.Flux_Analysis.Model_Loading.Functions as lf

# Takes in a list of symbol gene names and a HGNC checker .csv file and gives a dictionary to convert
# The list to the HGNC nomenclature.

# Navigates to main folder
os.chdir("../../..")

# The below csv file should be obtained by creating a list of all genes and running it in the hgnc checker
# (https://www.genenames.org/tools/multi-symbol-checker/). A case sensitive search should be done,
# as case helps differentiate between species
with open("Data/Input/nomenclature_files/HGNC_multi-symbol_checker_results/A549_HGNC_check_results.csv","r") as HGNC_check:
	file_list = HGNC_check.readlines()
	
# The list of gene names in the scRNAseq data
scRNAseq_gene_list = list(lf.RNA_Seq_Load("Data/Input/scRNA-seq-data/A549_sciCAR_data.mat")["RNA"]["Features"])

# breaks the HGNC check file into column header and a list of search results
header = file_list[1].replace('"','').split(",")
search_results = [file_list[i].replace('"','').split(",") for i in range(2,len(file_list))]

# specifies the priority order of match types from left to right
# So if a gene has multiple hits the one with the most preferred match will be used
# Currently ambiguity will occur if a gene has two or more hits for its most preferred match type
match_type_priority_list = ["Approved symbol","Alias symbol","Previous symbol"]

# Negative index is used because the 4th column of the csv obtained from the HGNC gives gene names
# Something which can contain commas
HGNC_negative_index = header.index("HGNC ID")-len(header)
# Generates a dictionary to go from gene symbols in the scRNAseq list to HGNC nomenclature
to_HGNC_dict = {}
for gene_name in scRNAseq_gene_list:
	# Creates list of hits, each hit occurs when a gene from the scRNAseq list matches the input section of the
	# checker. This is needed as single gene searching often yield multiple matches
	hit_list = []
	for gene_search_entry in search_results:
		if gene_name == gene_search_entry[0]:
			hit_list.append([gene_search_entry[1],search_results.index(gene_search_entry)])
	# This happens if the HGNC checker does not like the name of the gene
	if len(hit_list) == 0:
		print(f"There were no hits for {gene_name}")
	# This is the typical outcome, one hit for one gene
	elif len(hit_list) == 1:
		to_HGNC_dict[gene_name] = search_results[hit_list[0][1]][HGNC_negative_index]
	# For multiple hits, runs through the match type priority list and assigns HGNC symbol with best match
	else:
		for match_type in match_type_priority_list:
			if gene_name not in to_HGNC_dict.keys():
				for hit in hit_list:
					if search_results[hit[1]][1] == match_type:
						to_HGNC_dict[gene_name] = search_results[hit[1]][HGNC_negative_index]
		if gene_name not in to_HGNC_dict.keys():
			print("Multiple hits found, but no desired match type:")
			print(hit_list)
			to_HGNC_dict[gene_name] = search_results[hit_list[0][1]][HGNC_negative_index]
# Writes dictionary to file
with open("Data/Input/nomenclature_files/to_hgnc_dictionary/to_hgnc_A549_genes.json","w") as to_HGNC:
	json.dump(to_HGNC_dict,to_HGNC)
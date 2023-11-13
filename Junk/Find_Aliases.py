import pickle
import os
import Flux_Code.Programs.Flux_Analysis.Model_Loading.Functions as lf
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import numpy as np
import mygene
import sys
import json
from Bio import Entrez

#os.chdir("../../..")

def retrieve_annotation(id_list):
    Entrez.email = "agoetz@ufl.edu"
    """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information.
    Returns a list of dictionaries with the annotations."""

    request = Entrez.epost("gene", id=",".join(id_list))
    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        # FIXME: How generate NAs instead of causing an error with invalid IDs?
        print("An error occurred while retrieving the annotations.")
        print("The error returned was %s" % e)
        sys.exit(-1)

    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key=queryKey)
    annotations = Entrez.read(data)

    print("Retrieved %d annotations for %d genes" % (len(annotations), len(id_list)))

    return annotations



#with open("Data/Input/Symbol_Handling/Model_scRNA_mismatch_symbols/recon2_2","r") as mismatch_file:
#    symbol_list = eval(mismatch_file.readlines()[0])

#Converts symbols to entrez id
mg = mygene.MyGeneInfo()
mygene_search_results = mg.querymany(['ENSMUSG00000051951','ENSMUSG00000089699'])
#entrez_list = [i["entrezgene"] for i in mygene_search_results]
print(mygene_search_results)
input()
# *Always* tell NCBI who you are
Entrez.email = "agoetz@ufl.edu"
print(retrieve_annotation(['ENSMUSG00000051951','ENSMUSG00000089699']))
#annotations = json.dumps(retrieve_annotation(entrez_list))
#with open("Data/Input/Symbol_Handling/Model_scRNA_mismatch_database_dictionary/recon2_2.json","w") as writefile:
#    writefile.write(annotations)
from Flux_Code.Programs.Flux_Analysis.Classes_And_Functions.Flux_Model_Class import Flux_Balance_Model
import scipy.sparse as sp
import numpy as np
import scipy.io as spio

def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

def create_model_from_file(input_filename,model_name,output_filename):
    # filename should contain full path to get from cwd to file
    # model_name should be the key for the model in the dictionary of the mat file after it is opened by scipy.io
    # This does assume there is a convention to the name of files and that might not hold
    # New conventions should just be added in if possible
    mat = loadmat(input_filename)
    try:
        S = sp.csr_matrix(mat[model_name]["S"]).astype(np.float)
        lb = mat[model_name]["lb"].astype(np.float)
        ub = mat[model_name]["ub"].astype(np.float)
        c = mat[model_name]["c"].astype(np.float)
        b = mat[model_name]["b"].astype(np.float)
        print(mat[model_name].keys())
        #print("rules")
        #print(mat[model_name]["rules"])
        
        print(mat[model_name].keys())

        print(mat[model_name]["mets"])
        for i in mat[model_name]["mets"]:
            if "lac_L" in i:
                print(i)
        print(mat[model_name]["metNames"])
        input()
        print(mat[model_name]["description"])
        input()
        rxnGeneMat = sp.csr_matrix(mat[model_name]["rxnGeneMat"]).astype(np.float)
        grRules = list(mat[model_name]["grRules"])
        print(len(grRules))
        for i in range(len(grRules)):
            grRules[i] = str(grRules[i])
        for i in grRules:
            if i != "[]":
                print(i)
        gene_list = list(mat[model_name]["genes"])

        # I use rxns rather than reaction names to conserve model space. Currently this is somewhat important as multithreading involves handing several models at once
        # It should be beyond possible to add a small section to allow for retriving descriptions in the analysis section
        reaction_names = list(mat[model_name]["rxns"])
        species_names = list(mat[model_name]["mets"])
        species_comp = list(mat[model_name]["metFormulas"])
        input("stop")
        FBA_model = Flux_Balance_Model(model_name, S, lb, ub, b, reaction_names=reaction_names,species_names=species_names, species_comp=species_comp, grRules = grRules,genes = gene_list)
    except KeyError:
        print(f"It is likely model name was not valid, valid model name will likely be one of the following:")
        for key in mat.keys():
            print(key)
        print("If the model is valid, then check to see if your model uses standard names for model components ('S' for stochiometric matrix, 'lb' for lower bound list, etc.)")



def RNA_Seq_Load(input_filename):
    # filename should contain full path to get from cwd to file
    # model_name should be the key for the model in the dictionary of the mat file after it is opened by scipy.io
    # This does assume there is a convention to the name of files and that might not hold
    # New conventions should just be added in if possible
    mat = loadmat(input_filename)
    return mat
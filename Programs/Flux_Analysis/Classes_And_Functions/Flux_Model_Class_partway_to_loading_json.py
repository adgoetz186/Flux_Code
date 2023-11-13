import pickle
import sys
import sympy as sy
from pympler import asizeof
import gurobipy as gp
import copy
import gzip
import re
from pathlib import Path
import pandas as pd
from gurobipy import GRB
import scipy as sc
import scipy.linalg
import time
import json
import scipy.sparse as sp
import numpy as np
import igraph as ig
print(ig.__version__)
import scipy.io as spio
import scipy.linalg as sla
import random
from matplotlib import pyplot as plt


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


def size_order_key(e):
	return len(e)


class Flux_Balance_Model:
	# Im not adding an add species option, species should all be defined during initialization
	
	model = "none"
	
	def __init__(self, model_name=None, rxn_dict=None, met_dict=None, gene_dict = None):
		# model name should specify the exact model, for example the 3rd B6 mouse's model should be "B6-3" and
		# not "recon1"
		self.model_name = model_name
		
		# dictionary describing reactions, indexed by reaction name
		# should contain:
		# "S" key which specifies which metabolites are created/consumed in the reaction
		# "lb" and "ub" keys which specify bounds
		# optionally "grRules" key which gives the genes which facilitate the reaction
		self.rxn_dict = rxn_dict
		
		# dictionary describing metabolites, indexed by metabolite name
		# should contain:
		# "comp" key which specifies metabolite composition
		# "b" key which gives corresponding value of b in Sv = b
		# optionally "grRules" which gives the genes which facilitate the reaction
		self.met_dict = met_dict
		self.env = None
		self.gene_dict = None
		self.total_flux_limits = None
	
	def add_gp_key_env_to_model(self, TokenServer_name):
		# for hpg, TokenServer_name should be 'grb-ts.ufhpc'
		env = gp.Env(empty=True)
		env.setParam('TokenServer', TokenServer_name)
		env.start()
		self.env = env
	
	def save_model_as_pkl(self, filepath):
		dict_of_model = {}
		dict_of_model["model_name"] = self.model_name
		dict_of_model["total_flux_limits"] = self.total_flux_limits
		dict_of_model["rxn_dict"] = self.rxn_dict
		dict_of_model["met_dict"] = self.met_dict
		with open(filepath.parent/(filepath.stem+".pkl"), "wb") as outfile:
			pickle.dump(dict_of_model, outfile)
	
	def load_pkl_model(self, filepath):
		# Loads a pickled model which was saved using the "save_model_as_pkl" function from this class
		with open(filepath, "rb") as readfile:
			model_dict = pickle.load(readfile)
		self.model_name = model_dict["model_name"]
		# This should not be needed but allows models created before total flux limits were added to be loaded
		if "total_flux_limits" in model_dict.keys():
			self.total_flux_limits = model_dict["total_flux_limits"]
		self.rxn_dict = model_dict["rxn_dict"]
		self.met_dict = model_dict["met_dict"]
	
	def load_mat_model(self, filename, read_model_name="", model_comp=None,gz = False):
		# if model_name == "" the program guesses what the model key is by looking at the size of all dict entries
		# the largest entry is likely to be the one that holds the model
		# note sys.getsizeof() is insufficient as the model is often a dictionary which references large objects
		# and sys.getsizeof() will omit those referenced values
		# model_comp: dict which contains the names for each model component with the key being the "standard name"
		#   standard names: Stochiometric matrix = "S", lower bounds = "lb", upper bounds = "ub",
		#   rhs of Sv = b constraint = "b", dictionary of pinched reactions = "pinched_reactions",
		#   list of reaction names = "reaction_names", gene reaction rule list = "grRules",
		#   list of gene names = "genes", list of metabolite names = "metabolite_names"
		#   list of metabolite compositions = "metabolite_comp"
		#   In the case where an optional model component is not expected to be found in the model being loaded
		#   the key should remain unchanged, but the corresponding value should be None
		#   This will be assumed to be the case for pinched_reactions
		if gz:
			with gzip.open(filename, 'rb') as f:
				mat = loadmat(f)
		else:
			mat = loadmat(filename)
		if read_model_name == "":
			size_list = [asizeof.asizeof(mat[key]) for key in mat.keys()]
			read_model_name = list(mat.keys())[size_list.index(max([asizeof.asizeof(mat[key]) for key in mat.keys()]))]
		if model_comp == None:
			model_comp = {}
			model_comp["S"] = "S"
			model_comp["lb"] = "lb"
			model_comp["ub"] = "ub"
			model_comp["b"] = "b"
			model_comp["pinched_reactions"] = None
			model_comp["metabolite_names"] = "metabolite_names"
			model_comp["reaction_names"] = "reaction_names"
			model_comp["grRules"] = "grRules"
			model_comp["genes"] = "genes"
			model_comp["metabolite_comp"] = "metabolite_comp"
		try:
			test_key = mat[read_model_name]
		except KeyError:
			print(f"It is likely model name was not valid, valid model name will likely be one of the following:")
			for key in mat.keys():
				print(key)
		req_keys = {"S", "lb", "ub", "b", "pinched_reactions", "reaction_names", "grRules", "genes", "metabolite_names",
		            "metabolite_comp"}
		curr_keys = set(model_comp.keys())
		if len(req_keys - curr_keys) > 0:
			print(f"There are {len(req_keys - curr_keys)} required key(s) missing from model_comp")
			for i in (req_keys - curr_keys):
				print(i)
		
		try:
			model_dict = {}
			tried_key = model_comp["S"]
			S = sp.csr_matrix(mat[read_model_name][tried_key]).astype(np.float)
			S = S.toarray()
			model_dict["S"] = S
			tried_key = model_comp["lb"]
			lb = mat[read_model_name][tried_key].astype(np.float)
			model_dict["lb"] = lb
			tried_key = model_comp["ub"]
			ub = mat[read_model_name][tried_key].astype(np.float)
			model_dict["ub"] = ub
			tried_key = model_comp["b"]
			b = mat[read_model_name][tried_key].astype(np.float)
			model_dict["b"] = b
			
			tried_key = model_comp["reaction_names"]
			if tried_key == None:
				if S != None:
					reaction_names = [f"Reaction_{i}" for i in range(np.shape(self.S)[1])]
				else:
					reaction_names = None
			else:
				reaction_names = list(mat[read_model_name][tried_key])
			model_dict["reaction_names"] = reaction_names
			
			tried_key = model_comp["grRules"]
			if tried_key == None:
				grRules = None
			else:
				grRules = list(mat[read_model_name][tried_key])
			model_dict["grRules"] = grRules
			
			tried_key = model_comp["metabolite_names"]
			if tried_key == None:
				if S != None:
					metabolite_names = [f"metabolite_{i}" for i in range(np.shape(self.S)[0])]
				else:
					metabolite_names = None
			else:
				metabolite_names = list(mat[read_model_name][tried_key])
			model_dict["metabolite_names"] = metabolite_names
			
			tried_key = model_comp["metabolite_comp"]
			if tried_key == None:
				metabolite_comp = None
			else:
				metabolite_comp = list(mat[read_model_name][tried_key])
			model_dict["metabolite_comp"] = metabolite_comp
		
		except KeyError:
			print(f"The model uses a different key for {tried_key}")
			print(f"possible keys: {mat[read_model_name].keys()}")
		return model_dict
	
	def load_json_model(self, filename, rxn_comp=None,met_comp = None, gene_comp = None, keep_extra = True, gz=False):
		# if model_name == "" the program guesses what the model key is by looking at the size of all dict entries
		# the largest entry is likely to be the one that holds the model
		# note sys.getsizeof() is insufficient as the model is often a dictionary which references large objects
		# and sys.getsizeof() will omit those referenced values
		# model_comp: dict which contains the names for each model component with the key being the "standard name"
		#   standard names: Stochiometric matrix = "S", lower bounds = "lb", upper bounds = "ub",
		#   rhs of Sv = b constraint = "b", dictionary of pinched reactions = "pinched_reactions",
		#   list of reaction names = "reaction_names", gene reaction rule list = "grRules",
		#   list of gene names = "genes", list of metabolite names = "metabolite_names"
		#   list of metabolite compositions = "metabolite_comp"
		#   In the case where an optional model component is not expected to be found in the model being loaded
		#   the key should remain unchanged, but the corresponding value should be None
		#   This will be assumed to be the case for pinched_reactions
		if rxn_comp is None:
			rxn_comp = {"rxn_name" : "id", "rxn_metabolites":"metabolites", "lb" : "lower_bound", "ub" : "upper_bound", "pinched_reactions" : None,"subsystem":"subsystem", "grRules" : "gene_reaction_rule"}
			rxn_comp={"rxn_metabolites":"metabolites", "lb" : "lower_bound", "ub" : "upper_bound", "b" : "b","pinched_reactions" : None,"subsystem":"subsystem", "met_name" : "id", "rxn_name" : "id", "grRules" : "gene_reaction_rule", "met_composition" : "formula","met_compartment": "compartment"}
		if met_comp is None:
			met_comp = {"met_name" : "id","b" : "b", "met_composition" : "formula","met_compartment": "compartment"}
		if gene_comp is None:
			gene_comp = {"gene_name": "id", "symbol": "name"}
		with open(filename,) as json_in:
			model = json.load(json_in)
		print(model.keys())
		print(model["genes"])
		for i in model["genes"]:
			print(i)
			input()
		input()
		self.rxn_dict = {}
		for i in model["reactions"]:
			self.rxn_dict[i[model_comp["rxn_name"]]] = {}
			for key in model_comp.keys():
				if model_comp[key] in i and key != "rxn_name":
					self.rxn_dict[i[model_comp["rxn_name"]]][key] = i[model_comp[key]]
			extra = {}
			for key_orig in i.keys():
				if key_orig not in [model_comp[key_new] for key_new in self.rxn_dict[i[model_comp["rxn_name"]]].keys()] and key_orig != model_comp["rxn_name"]:
					extra[key_orig] = i[key_orig]
			if keep_extra:
				self.rxn_dict[i[model_comp["rxn_name"]]]["extra"] = extra
		self.met_dict = {}
		for i in model["metabolites"]:
			self.met_dict[i[model_comp["met_name"]]] = {}
			for key in model_comp.keys():
				if model_comp[key] in i and key != "met_name":
					self.met_dict[i[model_comp["met_name"]]][key] = i[model_comp[key]]
			extra = {}
			for key_orig in i.keys():
				if key_orig not in [model_comp[key_new] for key_new in
				                    self.met_dict[i[model_comp["met_name"]]].keys()] and key_orig != model_comp[
					"met_name"]:
					extra[key_orig] = i[key_orig]
			if keep_extra:
				self.met_dict[i[model_comp["met_name"]]]["extra"] = extra
			print(i)
			print(self.met_dict[i[model_comp["met_name"]]])
			input()
		input()
		if read_model_name == "":
			size_list = [asizeof.asizeof(mat[key]) for key in mat.keys()]
			read_model_name = list(mat.keys())[size_list.index(max([asizeof.asizeof(mat[key]) for key in mat.keys()]))]
		if model_comp == None:
			model_comp = {}
			model_comp["S"] = "S"
			model_comp["lb"] = "lb"
			model_comp["ub"] = "ub"
			model_comp["b"] = "b"
			model_comp["pinched_reactions"] = None
			model_comp["metabolite_names"] = "metabolite_names"
			model_comp["reaction_names"] = "reaction_names"
			model_comp["grRules"] = "grRules"
			model_comp["genes"] = "genes"
			model_comp["metabolite_comp"] = "metabolite_comp"
		try:
			test_key = mat[read_model_name]
		except KeyError:
			print(f"It is likely model name was not valid, valid model name will likely be one of the following:")
			for key in mat.keys():
				print(key)
		req_keys = {"S", "lb", "ub", "b", "pinched_reactions", "reaction_names", "grRules", "genes", "metabolite_names",
		            "metabolite_comp"}
		curr_keys = set(model_comp.keys())
		if len(req_keys - curr_keys) > 0:
			print(f"There are {len(req_keys - curr_keys)} required key(s) missing from model_comp")
			for i in (req_keys - curr_keys):
				print(i)
		
		try:
			model_dict = {}
			tried_key = model_comp["S"]
			S = sp.csr_matrix(mat[read_model_name][tried_key]).astype(np.float)
			S = S.toarray()
			model_dict["S"] = S
			tried_key = model_comp["lb"]
			lb = mat[read_model_name][tried_key].astype(np.float)
			model_dict["lb"] = lb
			tried_key = model_comp["ub"]
			ub = mat[read_model_name][tried_key].astype(np.float)
			model_dict["ub"] = ub
			tried_key = model_comp["b"]
			b = mat[read_model_name][tried_key].astype(np.float)
			model_dict["b"] = b
			
			tried_key = model_comp["reaction_names"]
			if tried_key == None:
				if S != None:
					reaction_names = [f"Reaction_{i}" for i in range(np.shape(self.S)[1])]
				else:
					reaction_names = None
			else:
				reaction_names = list(mat[read_model_name][tried_key])
			model_dict["reaction_names"] = reaction_names
			
			tried_key = model_comp["grRules"]
			if tried_key == None:
				grRules = None
			else:
				grRules = list(mat[read_model_name][tried_key])
			model_dict["grRules"] = grRules
			
			tried_key = model_comp["metabolite_names"]
			if tried_key == None:
				if S != None:
					metabolite_names = [f"metabolite_{i}" for i in range(np.shape(self.S)[0])]
				else:
					metabolite_names = None
			else:
				metabolite_names = list(mat[read_model_name][tried_key])
			model_dict["metabolite_names"] = metabolite_names
			
			tried_key = model_comp["metabolite_comp"]
			if tried_key == None:
				metabolite_comp = None
			else:
				metabolite_comp = list(mat[read_model_name][tried_key])
			model_dict["metabolite_comp"] = metabolite_comp
		
		except KeyError:
			print(f"The model uses a different key for {tried_key}")
			print(f"possible keys: {mat[read_model_name].keys()}")
		return model_dict
	
	def create_model_from_components(self,model_name,S,lb,ub,b,rxn_names,met_names,met_comp = None,grRules = None):
		self.model_name = model_name
		rxn_dict = {}
		for i in range(len(rxn_names)):
			rxn_sub_dict = {}
			rxn_met = [met_names[j] for j in np.nonzero(S[:, i])[0]]
			stochio_values_met = [S[:, i][j] for j in np.nonzero(S[:, i])[0]]
			rxn_sub_dict["S"] = dict(zip(rxn_met, stochio_values_met))
			rxn_sub_dict["lb"] = lb[i]
			rxn_sub_dict["ub"] = ub[i]
			if grRules != None:
				rxn_sub_dict["grRule"] = grRules[i]
			rxn_dict[rxn_names[i]] = rxn_sub_dict
		self.rxn_dict = rxn_dict
		met_dict = {}
		for i in range(len(met_names)):
			met_sub_dict = {}
			if met_comp != None:
				met_sub_dict["comp"] = met_comp[i]
			met_sub_dict["b"] = b[i]
			met_dict[met_names[i]] = met_sub_dict
		self.met_dict = met_dict
	
	def update_reaction_bounds(self, name, new_lb, new_ub):
		if not new_lb == "keep":
			self.rxn_dict[name]["lb"] = new_lb
		if not new_ub == "keep":
			self.rxn_dict[name]["ub"] = new_ub
		if self.rxn_dict[name]["ub"] < self.rxn_dict[name]["lb"]:
			print(new_ub, new_lb)
			print(f"WARNING! Upper bound smaller than lower bound at {name}")
			raise ValueError
	
	def dicts_to_mats(self,alphabetize = True):
		# if omit zero flux is not positive all fluxes are used, if it is positive then if a lb and ub are within
		# that distance the reaction will not be included
		# lists hold teh order of the final metabolites and reactions
		# sometimes this will matter, other times it will not
		# it is better to make it unimportant when possible
		rxn_list = []
		for i in self.rxn_dict.keys():
			rxn_list.append(i)
		met_list = list(self.met_dict.keys())
		if alphabetize:
			rxn_list.sort()
			met_list.sort()
		lb = np.empty((len(rxn_list)))
		ub = np.empty((len(rxn_list)))
		S = np.zeros((len(met_list),len(rxn_list)))
		for rxn_ind in range(len(rxn_list)):
			lb[rxn_ind] = self.rxn_dict[rxn_list[rxn_ind]]['lb']
			ub[rxn_ind] = self.rxn_dict[rxn_list[rxn_ind]]['ub']
			for met in self.rxn_dict[rxn_list[rxn_ind]]['S'].keys():
				S[met_list.index(met),rxn_ind] = self.rxn_dict[rxn_list[rxn_ind]]['S'][met]
		b = np.empty((len(met_list)))
		for met_ind in range(len(met_list)):
			if "b" in self.met_dict[met_list[met_ind]]:
				b[met_ind] = self.met_dict[met_list[met_ind]]['b']
			else:
				b[met_ind] = 0
		return lb,ub,S,b,rxn_list,met_list
		
	def get_rxn_comp_flattened(self,comp_name,alphabetize_by_rxn_name = True):
		rxn_list = list(self.rxn_dict.keys())
		if alphabetize_by_rxn_name:
			rxn_list.sort()
		com_list = []
		for rxn_name in rxn_list:
			com_list.append(self.rxn_dict[rxn_name][comp_name])
		return com_list
	
	def limit_flux_bounds(self, absolute_max_value):
		for rxn in self.rxn_dict.keys():
			self.update_reaction_bounds(rxn, -1 * absolute_max_value, absolute_max_value)
	
	def add_reaction(self,rxn_name, S_dict, lb, ub, grRule = None):
		# S_dict must be a dict with the metabolites appearing in the reaction as keys and their coefficients as values
		# Assumes all metabolites of reaction are in model
		if grRule == None:
			print("Gene Reaction Rule should be included, currently assuming no gene association with this reaction")
			self.rxn_dict[rxn_name] = {"S":S_dict,"lb":lb,"ub":ub,"grRule":"[]"}
		else:
			self.rxn_dict[rxn_name] = {"S": S_dict, "lb": lb, "ub": ub, "grRule": grRule}
	
	def convert_scRNAseq_to_RAS_matrix(self, scRNAseq_matrix, gene_name_list):
		# gene_name_list should correspond to the row/column in the scRNAseq_matrix
		RAS_matrix = np.empty((len(self.reaction_names), np.shape(scRNAseq_matrix)[1]))
		for ci in range(np.shape(scRNAseq_matrix)[1]):
			gene_value_dict = dict(zip(gene_name_list, scRNAseq_matrix[:, ci]))
			for ri in range(len(self.reaction_names)):
				# Fixes notation error in recon2.2 model where HGNC:987 and HGNC:2898 appear as HGNC:HGNC:987 and
				# HGNC:HGNC:2898 (respectively)
				print(self.grRules[ri])
				print(self.reaction_names[ri])
				grRule_fixed = self.grRules[ri].replace("HGNC:HGNC:", "HGNC:")
				if self.reaction_names[ri] == "ACCOAC":
					print(grRule_fixed)
					print(self.get_reaction_activity_scores(grRule_fixed, gene_value_dict))
					input()
				RAS_matrix[ri, ci] = self.get_reaction_activity_scores(grRule_fixed, gene_value_dict)
				print(ri, ci)
		return RAS_matrix
	
	def get_reaction_activity_scores(self, string_to_handle, gene_value_dictionary):
		# additional arguements will be needed for different notation systems
		# gene value dictionary should be in the form "gene name" and "abundance"
		if string_to_handle == "[]":
			return 0
		string_to_evaluate = string_to_handle.replace("(", "{").replace(")", "}")
		vart = True
		while vart:
			m = re.search('{([^{}]+)}', string_to_evaluate)
			if m is not None:
				if "or" in m.group(0):
					new_str = m.group(0).replace("or", "+").replace("{", "(").replace("}", ")")
				if "and" in m.group(0):
					new_str = m.group(0).replace("and", ",").replace("{", "min(").replace("}", ")")
				string_to_evaluate = string_to_evaluate.replace(m.group(0), new_str)
			else:
				if "or" in string_to_evaluate:
					string_to_evaluate = "(" + string_to_evaluate.replace("or", "+") + ")"
				if "and" in string_to_evaluate:
					string_to_evaluate = "min(" + string_to_evaluate.replace("and", ",") + ")"
				vart = False
		gene_list = list(gene_value_dictionary.keys())
		
		gene_list.sort(reverse=True, key=size_order_key)
		for i in gene_list:
			string_to_evaluate = string_to_evaluate.replace(i, str(gene_value_dictionary[i]))
		return eval(string_to_evaluate)
	
	def flux_balance_test(self, point, cons_tol):
		# remove false
		if np.sum(np.matmul(self.S.toarray(), np.transpose(point)) - self.b) > -cons_tol and np.sum(
				np.matmul(self.S.toarray(), np.transpose(point)) - self.b) < cons_tol and False:
			print(f"Error of point: {np.abs(np.sum(np.matmul(self.S.toarray(), np.transpose(point))) - self.b)}")
	
	def normalize(self, vec):
		norm = np.linalg.norm(vec)
		if norm == 0:
			return vec
		else:
			return vec / norm
	
	def get_grRule_list(self,alphabetize = True):
		lb,ub,S,b,rxn_list,met_list = self.dicts_to_mats(alphabetize=alphabetize)
		gr_Rule_list = []
		for i in rxn_list:
			if len(self.rxn_dict[i]["grRule"]) == 0:
				gr_Rule_list.append("")
			else:
				gr_Rule_list.append(self.rxn_dict[i]["grRule"])
		return gr_Rule_list
	
	def get_used_gene_list(self):
		gene_list = []
		gr_Rule_list = self.get_grRule_list()
		for i in gr_Rule_list:
			for id in re.findall("\d+\.\d+",i):
				if id not in gene_list:
					gene_list.append(id)
		return gene_list
	
	def grRule_to_RAS(self,grRule,RNA_seq_dict):
		pass
		
	def find_objective_value(self, c):
		env = gp.Env(empty=True)
		#env.setParam('TokenServer', 'grb-ts.ufhpc')
		env.start()
		model = gp.Model("opt_model_" + self.model_name, env=env)
		model.Params.LogToConsole = 0
		react_flux = model.addMVar(shape=self.S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=self.lb,
		                           ub=self.ub)
		obj = np.reshape(c, (1, -1))
		#print(obj.shape)
		#print(obj @ react_flux)
		model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		rhs = np.transpose(self.b)
		model.addConstr(self.S @ react_flux == rhs, name="c")
		model.optimize()
		if model.Status == 2:
			return model.objVal
		else:
			print(model.Status)
			return np.nan
	
	def fva(self):
		# Returns list of [min,max] elements
		# If reactions are edited this fva no longer applies
		# This is valid for removed reactions, but could lead to issues if 2 reactions are simply switched around
		lb,ub,S,b,rxn_list,met_list = self.dicts_to_mats()
		if self.env != None:
			model = gp.Model("fva_model_" + self.model_name, env=self.env)
		else:
			model = gp.Model("fva_model_" + self.model_name)
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb,ub=ub)
		model.addConstr(S @ react_flux == b, name="c")
		model.Params.LogToConsole = 0
		obj = np.zeros((1, S.shape[1]))
		variability_list = [{"lb":0,"ub":0} for i in range(len(lb))]
		for i in range(len(lb)):
			print(len(lb)-i)
			min_max = []
			obj = np.zeros(obj.shape)
			obj[0, i] = 1
			model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
			model.optimize()
			min_max.append(model.objVal)
			model.reset()
			obj[0, i] = -1
			model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
			model.optimize()
			min_max.append(model.objVal*-1)
			variability_list[i]["ub"] = max(min_max)
			variability_list[i]["lb"] = min(min_max)
		return dict(zip(rxn_list, variability_list))
	
	def test_feasibility(self, lb = None, ub = None, S = None, b = None,**kwargs):
		if lb is None:
			lb,ub,S,b,rxn_list,met_list = self.dicts_to_mats()
		if self.env != None:
			model = gp.Model("fva_model_" + self.model_name, env=self.env)
		else:
			model = gp.Model("fva_model_" + self.model_name)
		if "unique_bounds" in kwargs.keys():
			tflb = kwargs["unique_bounds"][0]
			tfub = kwargs["unique_bounds"][1]
		else:
			tflb = lb
			tfub = ub
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=tflb, ub=tfub)
		model.addConstr(S @ react_flux == b, name="c")
		obj = np.zeros((1, S.shape[1]))
		model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.Params.LogToConsole = 0
		#model.Params.NumericFocus = 3
		#model.Params.ScaleFlag = 0
		#model.Params.FeasibilityTol = 1e-9
		model.optimize()
		return model.Status == 2
	
	def test_feasibility_pass_model(self, lb, ub, model, react_flux):
		react_flux.lb = lb
		react_flux.ub = ub
		obj = np.zeros((1, self.S.shape[1]))
		model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.Params.LogToConsole = 0
		model.optimize()
		return model.Status == 2
	
	def initalize_model(self, **kwargs):
		# Performs FVA and updates lb and ub
		varibility_list = self.fva(**kwargs)
		
		for i in range(len(varibility_list)):
			self.update_reaction_bounds(i, varibility_list[i][0], varibility_list[i][1])
		
		if "shrink_S" in kwargs.keys():
			self.pinch_reactions(kwargs["shrink_S"])
	
	def pinch_reactions(self, distance_between_bounds, zero_out=False):
		# This is unstable and should only be used as one of the final steps of model refinement
		# Improvements could be made such as using a dictionary for the metabolite constant storage
		constant_var_index = np.where((self.ub - self.lb) <= distance_between_bounds)[0]
		variable_var_index = np.where((self.ub - self.lb) > distance_between_bounds)[0]
		
		contant_var_val = (self.ub[constant_var_index] + self.lb[constant_var_index]) / 2
		
		cS = self.S.toarray()[:, constant_var_index]
		self.b -= np.matmul(cS, contant_var_val)
		constant_var_index = constant_var_index.tolist()
		constant_var_index.sort(reverse=True)
		for i in constant_var_index:
			metabolite_stochio_storage = {}
			row_values = self.S[:, i].nonzero()[0]
			for j in row_values:
				metabolite_stochio_storage[self.metabolite_names[j]] = self.S[j, i]
			self.pinched_reactions[self.reaction_names[i]] = [(self.lb[i] + self.ub[i]) / 2, metabolite_stochio_storage]
			self.del_reaction(i)
		print(f"Restricted {len(constant_var_index)} of {len(constant_var_index) + len(variable_var_index)} fluxes")
	
	def get_exchanged_reaction_info(self,true_exchange_signifier = ""):
		exch_rxn = []
		exch_met = []
		for i in self.rxn_dict.keys():
			if len(self.rxn_dict[i]["S"]) == 1 and true_exchange_signifier in i:
				exch_rxn.append(i)
				exch_met.append(list(self.rxn_dict[i]["S"].keys())[0])
		return dict(zip(exch_rxn,exch_met))
	
	def metabolite_info(self, name, flux_dtf = None):
		print(f"Information on {name}:")
		print(f"Molar mass = {np.round(self.get_met_mol_weight(name),2)} g/mol")
		print(self.met_dict[name])
		print(f"The reactions in which {name} is produced or consumed are:")
		list_of_con = []
		for i in self.rxn_dict.keys():
			if name in self.rxn_dict[i]["S"].keys():
				print(i)
				if flux_dtf is not None and i in flux_dtf.columns:
					print(f"Rxn name: {i}, average met produced: {self.rxn_dict[i]['S'][name]*np.average(flux_dtf.loc[:,i].to_numpy())}")
					list_of_con.append(self.rxn_dict[i]['S'][name]*np.average(flux_dtf.loc[:,i].to_numpy()))
	
	def fast_make_positive_HR(self,point):
		pass
	
	def generate_mat_wo_exchange(self,prune_specific = []):
		# convert this back to dict for check to make sure you didnt mess with any reactions
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		rxn_array = np.array(rxn_list)
		iti = [not(np.size(np.nonzero(S[:, rxn_ind])) == 1 or rxn_list[rxn_ind] in prune_specific) for rxn_ind in range(len(rxn_list))]
		internal_rxn = rxn_array[iti]
		internal_lb = lb[iti]
		internal_ub = ub[iti]
		internal_S = S[:,iti]
		return internal_lb, internal_ub, internal_S, b, internal_rxn, met_list

	def pos_S(self,inverse_tag = "_inverted"):
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		print(rxn_list)
		for i in range(len(rxn_list)):
			if "_inverted" in rxn_list[i]:
				positive_ent = copy.deepcopy(rxn_list[i]).replace("_inverted","")
				if positive_ent in rxn_list and rxn_list[i-1] != positive_ent:
					print(positive_ent)
					print(rxn_list[i])
					print(rxn_list[i-1])
					print(rxn_list[rxn_list.index(positive_ent)])
					input()
		print(np.array(rxn_list))
		input()

	
	def generate_null_space(self,internal_S):
		for i in range(9):
			rounded = np.sum(np.round(internal_S*10**i))
			if rounded == np.sum(internal_S*10**i):
				internal_Sy = sy.matrices.Matrix.zeros(np.shape(internal_S)[0], np.shape(internal_S)[1])
				for row_ind in range(np.shape(internal_S)[0]):
					for col_ind in range(np.shape(internal_S)[1]):
						if internal_S[row_ind, col_ind] != 0:
							internal_Sy[row_ind, col_ind] = sy.S(int(np.round(internal_S[row_ind, col_ind]*10**i)))/(sy.S(10)**sy.S(i))
				break
		Ny = np.hstack([np.array(i) for i in internal_Sy.nullspace()])
		print(Ny)
		return Ny
	
	def generate_null_space_fast(self,internal_S):
		Ny = sla.null_space(internal_S)
		#Ny = np.where(np.abs(Ny) > 1e-10,Ny,np.zeros_like(Ny))
		#print(NzNy)
		flat_array = np.log10(np.abs(Ny.flatten())+1e-16)
		return Ny
	
	def generate_therm_model(self,trans_N, point,u_bound_tol = 1e-7):
		if self.env != None:
			model = gp.Model("therm_model_" + self.model_name, env=self.env)
		else:
			model = gp.Model("therm_model_" + self.model_name)
		
		lb_mu = np.where(point > 0, -np.inf, u_bound_tol)
		lb_mu = np.where(point == 0, -np.inf, lb_mu)
		ub_mu = np.where(point > 0, -1*u_bound_tol, np.inf)
		ub_mu = np.where(point == 0, np.inf, ub_mu)
		#lb_negative_mu = np.where(point >= 0, -1000, tol)
		#ub_positive_mu = np.where(point <= 0, 1000, -1*tol)
		ts1 = trans_N.shape[1]
		mu = model.addMVar(shape=ts1, vtype=GRB.CONTINUOUS, name="mu",lb = lb_mu,ub=ub_mu)
		model.Params.LogToConsole = 0
		model.addConstr(trans_N @ mu == 0, name="c")
		
		obj = np.zeros((1, ts1))
		model.setObjective(obj @ mu, GRB.MAXIMIZE)
		model.optimize()
		return model.Status
	
	def test_thermo_feasibilty(self,N, therm_model):
		pass
		
	
	def element_edge_vert_generator(self,rxn,rxn_dtf = None,element = "compound",break_element_list = []):
		edge_list = []
		weight_list = []
		for met in self.rxn_dict[rxn]["S"].keys():
			if element == "compound":
				element_per_stochio = 1
			else:
				element_per_stochio = self.get_met_comp_dict(met)[element]
			if rxn_dtf is not None:
				avg_flux = np.average(rxn_dtf.loc[:,rxn].to_numpy())
			else:
				avg_flux = 1
			print(rxn,avg_flux)
			#input()
			# obtains the amount of element going to a specific product in the reaction
			if element_per_stochio * self.rxn_dict[rxn]["S"][met]*avg_flux < 0:
				if avg_flux < 0:
					if met in break_element_list:
						edge_list.append((met+"_for_"+rxn, "rev-" + rxn + "_RXN"))
					else:
						edge_list.append((met, "rev-" + rxn + "_RXN"))
				else:
					if met in break_element_list:
						edge_list.append((met+"_for_"+rxn,rxn+"_RXN"))
					else:
						edge_list.append((met,rxn+"_RXN"))
				weight_list.append(-1*element_per_stochio* self.rxn_dict[rxn]["S"][met]*avg_flux)
			elif element_per_stochio*self.rxn_dict[rxn]["S"][met]*avg_flux > 0:
				if avg_flux < 0:
					if met in break_element_list:
						edge_list.append(("rev-" + rxn + "_RXN", met+"_for_"+rxn))
					else:
						edge_list.append(("rev-" + rxn + "_RXN", met))
				else:
					if met in break_element_list:
						edge_list.append((rxn+"_RXN", met+"_for_"+rxn))
					else:
						edge_list.append((rxn + "_RXN", met))
				weight_list.append(element_per_stochio* self.rxn_dict[rxn]["S"][met]*avg_flux)
		return edge_list,weight_list
			#if self.rxn_dict[rxn]["S"][met] < 0:
			#	reactants[met] = np.average(rxn_dtf.loc[:,rxn].to_numpy())*-1*self.met_dict[met]['comp']
		
	def generate_rxn_met_graph(self,edge_list,weight_list,color_edge_of_vert = None,default_edge_color = "grey",vert_leak_tol = None,max_vec_width = 5,weight_display="raw"):
		g = ig.Graph(directed=True)
		edge_color_list = []
		vert_list = []
		for i in edge_list:
			if i[0] not in vert_list:
				vert_list.append(i[0])
			if i[1] not in vert_list:
				vert_list.append(i[1])
			if color_edge_of_vert is not None:
				if i[0].replace("_RXN","").replace("rev-","") in color_edge_of_vert.keys():
					edge_color_list.append(color_edge_of_vert[i[0].replace("_RXN","").replace("rev-","")])
				elif i[1].replace("_RXN", "").replace("rev-", "") in color_edge_of_vert.keys():
					edge_color_list.append(color_edge_of_vert[i[1].replace("_RXN", "").replace("rev-", "")])
				else:
					edge_color_list.append(default_edge_color)
			else:
				edge_color_list.append(default_edge_color)

		vert_leak_list_color = []
		for i in vert_list:
			vert_leak = 0
			if vert_leak_tol is not None:
				for j in range(len(edge_list)):
					if i == edge_list[j][0]:
						vert_leak-=weight_list[j]
					elif i == edge_list[j][1]:
						vert_leak += weight_list[j]
				if abs(vert_leak) > vert_leak_tol:
					if vert_leak < 0:
						vert_leak_list_color.append("red")
					else:
						vert_leak_list_color.append("blue")
				else:
					vert_leak_list_color.append("grey")
			else:
				vert_leak_list_color.append("grey")

		
				
		
		g.add_vertices(vert_list)
		g.vs["is_RXN"] = ["_RXN" not in i for i in vert_list]
		g.add_edges(edge_list)
		g.es["weight"] = weight_list
		g.es["color"] = edge_color_list
		g.vs["color"] = vert_leak_list_color
		
		
		

		print(edge_color_list)
		
		visual_style = {}
		visual_style["vertex_size"] = [float(i)/5+.01 for i in g.vs["is_RXN"]]
		visual_style["vertex_color"] = g.vs["color"]
		visual_style["edge_color"] = g.es["color"]
		visual_style["vertex_label"] = [i.replace("_RXN","") for i in g.vs["name"]]
		visual_style["edge_width"] = [i/max(list(g.es["weight"]))*max_vec_width for i in g.es["weight"]]
		if weight_display == "raw":
			visual_style["edge_label"] = [f"{np.round(i,3)}" for i in g.es["weight"]]
		elif weight_display == "rel":
			visual_style["edge_label"] = [f"{np.round(i/max(list(g.es['weight'])), 3)}" for i in g.es["weight"]]
		#visual_style["edge_color"] = [color_list[int(i)] for i in g.es["met_created"]]
		visual_style["layout"] = g.layout("kk")
		visual_style["bbox"] = (100000, 100000)
		visual_style["margin"] = [0,0,0,0]
		fig, ax = plt.subplots()
		ig.plot(g, target=ax, **visual_style)
		plt.show()
		
	def reaction_info(self, name):
		# Reaction name can be the name of the reaction or the index location of the reaction
		print(f"Information on {name}:")
		print(f"Lower Bound: {self.rxn_dict[name]['lb']}, Upper Bound {self.rxn_dict[name]['ub']}")
		print(f"The relevant chemical species participating in the reaction and their coefficients are:")
		print(self.rxn_dict[name]["S"])
		print(f"The reaction gene rule is:\n{self.rxn_dict[name]['grRule']}")
	
	def comparable_reaction_report(self, name):
		# Reaction name can be the name of the reaction or the index location of the reaction
		# Returns 2d list, the ith entry being ["metabolite_name_i",reaction_coefficient_for metabolite_i]
		if isinstance(name, str):
			name = self.reaction_names.index(name)
		metabolite_names = []
		row_values = self.S[:, name].nonzero()[0]
		for i in row_values:
			metabolite_names.append([self.metabolite_names[i], self.S[i, name]])
		return metabolite_names
	
	def objective_function_generator(self, list_of_reaction_names, list_of_weights):
		c = np.zeros(np.shape(self.S)[1])
		for i in range(len(list_of_reaction_names)):
			if isinstance(list_of_reaction_names[i], str):
				c[self.reaction_names.index(list_of_reaction_names[i])] = list_of_weights[i]
			else:
				c[list_of_reaction_names[i]] = list_of_weights[i]
		return c
		
		
	
	def lin_dep_np_ctp(self, old_on_basis, origin,accepted_points, new_point,old_rank):
		# Description:
		#   Takes in a set of basis vectors and a new proposed basis vector and either accepts or rejects
		#       the proposed vector based on whether the addition of the new vector increases the dimensionality of the
		#       basis set
		# Input:
		#   old_basis: 2d numpy array containing points which define a basis wrt the origin point
		#   origin: vectors are defined as the difference between points of interest and the origin
		#       Value should be within feasible space, either the first generated vertex or the center of the point
		#       matrix are valid approaches, we will use the first generated vertex
		#   new_point: proposed point which creates a vector wrt the origin
		# Output:
		#   either the old basis and old rank are returned or, if the proposed vector added to the basis creates
		#       a new larger basis, that basis and rank are returned
		
		# find the dimensionality of the original basis set
		#true_old_rank = np.linalg.matrix_rank(old_basis - origin)

		# make new candidate basis vector
		diff = (new_point - origin)
		potential_basic_vec = diff/np.sqrt(diff.dot(diff))

		# make new basis set
		if old_rank == 0:
			new_on_basis = np.zeros((1, np.size(origin)))
			new_on_basis[0] = potential_basic_vec
		else:
			new_on_basis = np.vstack((old_on_basis, potential_basic_vec))

		# find the dimensionality of the new basis set
		new_rank = np.linalg.matrix_rank(new_on_basis)

		#if len(np.shape(old_on_basis)) > 1:
		#	print(np.linalg.svd(old_on_basis)[1][-5:])
		#	print(np.linalg.svd(new_on_basis)[1][-5:])

			
		# if the new basis has a higher dimensionality, the added vector is indeed spanning and is kept
		if new_rank <= old_rank:
			return accepted_points ,old_on_basis, old_rank
		else:
			if np.size(accepted_points) == 0:
				accepted_points = np.reshape(new_point,(1,-1))
			else:
				accepted_points = np.vstack([accepted_points,new_point])
			return accepted_points,new_on_basis, new_rank
	
	
	def min_flux(self):
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		# Defines Model
		env = gp.Env(empty=True)
		# env.setParam('TokenServer', 'grb-ts.ufhpc')
		env.start()
		model = gp.Model("warmup_model_" + self.model_name, env=env)
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
		
		model.addConstr(S @ react_flux == b, name="c")
		
		obj = -1 * np.ones((1, np.size(ub)))
		model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.optimize()
		min_flux = model.objVal * -1
		
		return min_flux
	
	def generate_warmup_gb(self, max_search_number, sampling_levels, points_save_path, model_save_path,
	                          required_rank=np.inf, rxn_penalty_vec=None, **kwargs):
		# rel_distance_to_min_flux
		# Description:
		#   Generates warmup vectors for coordinate direction hit and run
		# Input:
		#   max_search_number: integer value which determines how many proposed vectors are subsequently rejected before
		#       the algorithm terminates
		# Output:
		#   2d array with each row giving a point in the feasible space the top row provides an origin point
		#       for the space while the following points gives a spanning set of vectors wrt the origin point
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		
		# Defines Model
		env = gp.Env(empty=True)
		# env.setParam('TokenServer', 'grb-ts.ufhpc')
		env.start()
		
		if rxn_penalty_vec is None:
			rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
		else:
			rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
		S = np.vstack((S, rxn_penalty_vec))
		S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
		lb = np.hstack((lb, 0))
		ub = np.hstack((ub, np.inf))
		b = np.hstack((b, 0))
		S[-1, -1] = -1
		rxn_list.append("total_flux")
		
		model = gp.Model("warmup_model_" + self.model_name, env=env)
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
		
		model.addConstr(S @ react_flux == b, name="c")
		
		obj = np.zeros((1, np.size(ub)))
		obj[0,-1] = -1
		print(obj)
		model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.Params.NumericFocus=3
		model.optimize()
		min_flux = model.objVal * -1
		print(min_flux)
		#input()
		N = sc.linalg.null_space(S)
		
		total_points = np.array([])
		origin = np.array([])
		
		# only 1 origin point should be made
		origin_created = False
		
		for sampling_level in sampling_levels:
			# Defines origin and points as empty
			
			accepted_points = np.array([])
			ortho_norm_basis = np.array([])
			
			ub[-1] = sampling_level * min_flux
			# if you want purely slices use
			# lb[-1] = sampling_level * min_flux
			model = gp.Model("warmup_model_" + self.model_name, env=env)
			react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
			# print(react_flux)
			model.addConstr(S @ react_flux == b, name="c")
			
			# Suppress printing to console
			model.Params.LogToConsole = 0
			
			# Initializes number of subsequent vectors which have been rejected
			search_number = 0
			
			#
			rank = 0
			
			print(f"SAMPLING LEVEL {sampling_level}")
			
			while (search_number < max_search_number and rank < required_rank):
				# print(search_number,rank)
				# Finds random vertex
				while model.Status != 2:
					model.reset()
					print(model.Status)
					obj = (2 * np.random.random_sample(S.shape[1]) - 1)
					# for variable objective total flux
					# obj[-1] = (sampling_level * min_flux*(np.heaviside(obj[-1],0)*2-1))
					model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
					model.optimize()
				
				# The below ensures the points fall inside the bounds
				# Does seem to potentially break constraints
				# entries of points are just being moved, the points should be scaled
				# if this becomes a problem, define scaling function
				# point_to_use = react_flux.X
				point_to_use = np.matmul(N, np.matmul(np.transpose(N), react_flux.X))
				point_to_use = np.minimum(point_to_use, ub)
				point_to_use = np.maximum(point_to_use, lb)
				
				#print(point_to_use)
				#print(np.max(point_to_use[:-1]*rxn_penalty_vec[0]))
				#print(rxn_list[np.argmax(point_to_use[:-1]*rxn_penalty_vec[0])])
				#print(rxn_penalty_vec[0,np.argmax(point_to_use[:-1]*rxn_penalty_vec[0])])
				#input()
				# print(point_to_use[-1])
				# Saves rank to later check for rank change
				old_rank = rank
				print(search_number, rank, point_to_use[-1], sampling_level * min_flux)
				# The first point becomes and origin point
				# The rank is 0 with the origin point since no vectors have been defined
				if origin_created:
					accepted_points, ortho_norm_basis, rank = self.lin_dep_np_ctp(ortho_norm_basis, origin,
					                                                              accepted_points, point_to_use, rank)
				if not origin_created:
					origin = point_to_use
					origin_created = True
				# Reset the model to discard the current solution
				model.reset()
				print(rank)
				# Reset the search_number if a vector is accepted as basis
				if rank != old_rank or np.size(accepted_points) == 0:
					search_number = 0
				else:
					search_number += 1
			print(rank)
			if np.size(total_points) == 0:
				total_points = np.vstack((origin, accepted_points))
			else:
				total_points = np.vstack((total_points, accepted_points))
		
		# returns set of points with top point being basis
		self.total_flux_limits = np.array([min_flux, max(sampling_levels) * min_flux])
		np.save(points_save_path / ("warmup_" + self.model_name), total_points)
		self.save_model_as_pkl(model_save_path / (self.model_name + "_HR_ready"))
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
	
	def generate_warmup_gb_cp(self, max_search_number, sampling_levels, points_save_path, model_save_path,
	                       required_rank=np.inf, rxn_penalty_vec=None, **kwargs):
		# rel_distance_to_min_flux
		# Description:
		#   Generates warmup vectors for coordinate direction hit and run
		# Input:
		#   max_search_number: integer value which determines how many proposed vectors are subsequently rejected before
		#       the algorithm terminates
		# Output:
		#   2d array with each row giving a point in the feasible space the top row provides an origin point
		#       for the space while the following points gives a spanning set of vectors wrt the origin point
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		
		# Defines Model
		env = gp.Env(empty=True)
		# env.setParam('TokenServer', 'grb-ts.ufhpc')
		env.start()
		model = gp.Model("warmup_model_" + self.model_name, env=env)
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
		
		model.addConstr(S @ react_flux == b, name="c")
		
		obj = -1 * np.ones((1, np.size(ub)))
		model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.optimize()
		min_flux = model.objVal * -1
		if rxn_penalty_vec is None:
			rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
		else:
			rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
		S = np.vstack((S, rxn_penalty_vec))
		S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
		lb = np.hstack((lb, 0))
		ub = np.hstack((ub, min_flux))
		b = np.hstack((b, 0))
		S[-1, -1] = -1
		rxn_list.append("total_flux")
		N = sc.linalg.null_space(S)
		
		total_points = np.array([])
		origin = np.array([])
		
		# only 1 origin point should be made
		origin_created = False
		
		for sampling_level in sampling_levels:
			# Defines origin and points as empty
			
			accepted_points = np.array([])
			ortho_norm_basis = np.array([])
			
			ub[-1] = sampling_level * min_flux
			# if you want purely slices use
			# lb[-1] = sampling_level * min_flux
			model = gp.Model("warmup_model_" + self.model_name, env=env)
			react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
			# print(react_flux)
			model.addConstr(S @ react_flux == b, name="c")
			
			# Suppress printing to console
			model.Params.LogToConsole = 0
			
			# Initializes number of subsequent vectors which have been rejected
			search_number = 0
			
			#
			rank = 0
			
			print(f"SAMPLING LEVEL {sampling_level}")
			
			while (search_number < max_search_number and rank < required_rank):
				# print(search_number,rank)
				# Finds random vertex
				while model.Status != 2:
					obj = (2 * np.random.random_sample(S.shape[1]) - 1)
					# for variable objective total flux
					# obj[-1] = (sampling_level * min_flux*(np.heaviside(obj[-1],0)*2-1))
					model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
					model.optimize()
				
				# The below ensures the points fall inside the bounds
				# Does seem to potentially break constraints
				# entries of points are just being moved, the points should be scaled
				# if this becomes a problem, define scaling function
				# point_to_use = react_flux.X
				point_to_use = np.minimum(react_flux.X, ub)
				point_to_use = np.maximum(point_to_use, lb)
				point_to_use = np.matmul(N, np.matmul(np.transpose(N), point_to_use))
				
				# print(point_to_use[-1])
				# Saves rank to later check for rank change
				old_rank = rank
				print(search_number, rank, point_to_use[-1], sampling_level * min_flux)
				# The first point becomes and origin point
				# The rank is 0 with the origin point since no vectors have been defined
				if origin_created:
					accepted_points, ortho_norm_basis, rank = self.lin_dep_np_ctp(ortho_norm_basis, origin,
					                                                              accepted_points, point_to_use, rank)
				if not origin_created:
					origin = point_to_use
					origin_created = True
				# Reset the model to discard the current solution
				model.reset()
				print(rank)
				# Reset the search_number if a vector is accepted as basis
				if rank != old_rank or np.size(accepted_points) == 0:
					search_number = 0
				else:
					search_number += 1
			print(rank)
			if np.size(total_points) == 0:
				total_points = np.vstack((origin, accepted_points))
			else:
				total_points = np.vstack((total_points, accepted_points))
		
		# returns set of points with top point being basis
		self.total_flux_limits = np.array([min_flux, max(sampling_levels) * min_flux])
		np.save(points_save_path / ("warmup_" + self.model_name), total_points)
		self.save_model_as_pkl(model_save_path / (self.model_name + "_HR_ready"))
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
	
	def test_warmup_gb(self, sample_size, save_path, required_rank=np.inf, rel_distance_to_min_flux=-1,
	                       **kwargs):
		# rel_distance_to_min_flux
		# Description:
		#   Generates warmup vectors for coordinate direction hit and run
		# Input:
		#   max_search_number: integer value which determines how many proposed vectors are subsequently rejected before
		#       the algorithm terminates
		# Output:
		#   2d array with each row giving a point in the feasible space the top row provides an origin point
		#       for the space while the following points gives a spanning set of vectors wrt the origin point
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		# Defines Model
		env = gp.Env(empty=True)
		# env.setParam('TokenServer', 'grb-ts.ufhpc')
		env.start()
		model = gp.Model("warmup_model_" + self.model_name, env=env)
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
		
		model.addConstr(S @ react_flux == b, name="c")
		
		if rel_distance_to_min_flux > 0:
			obj = -1 * np.ones((1, np.size(ub)))
			model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
			model.optimize()
			min_flux = model.objVal * -1
			max_allowed_flux = rel_distance_to_min_flux * min_flux
			print(max_allowed_flux)
			
			S = np.vstack((S, np.ones((1, np.shape(S)[1]))))
			S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
			lb = np.hstack((lb, 0))
			ub = np.hstack((ub, max_allowed_flux))
			b = np.hstack((b, 0))
			S[-1, -1] = -1
			rxn_list.append("total_flux")
			model = gp.Model("warmup_model_" + self.model_name, env=env)
			react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
			# print(react_flux)
			model.addConstr(S @ react_flux == b, name="c")
		
		# Suppress printing to console
		model.Params.LogToConsole = 0
		rxn_count = np.size(ub)
		rxn_ind = np.arange(rxn_count)
		# Initializes rank
		rank = -1
		average_of_average_flux = np.zeros(rxn_count)
		for i in range(rxn_count):
			print(i)
			# Finds random vertex
			
			af = 0
			for j in range(sample_size):
				
				while model.Status != 2:
					np.random.shuffle(rxn_ind)
					obj = (np.random.random_sample(S.shape[1]) - 1)
					obj[rxn_ind[:i]] *= -1
					model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
					model.optimize()
				point_to_use = np.minimum(react_flux.X, ub)
				point_to_use = np.maximum(point_to_use, lb)
				af += point_to_use[-1]
				
				# Reset the model to discard the current solution
				model.reset()
			average_of_average_flux[i] = af/sample_size
			print(af/sample_size)
			
			
			
		print(rank)
		# returns set of points with top point being basis
		np.save(save_path / (self.model_name), average_of_average_flux)
	
	def generate_vertex_samples(self, n,save_path,rel_distance_to_min_flux = -1, **kwargs):
		# Description:
		#   Generates warmup vectors for coordinate direction hit and run
		# Input:
		#   max_search_number: integer value which determines how many proposed vectors are subsequently rejected before
		#       the algorithm terminates
		# Output:
		#   2d array with each row giving a point in the feasible space the top row provides an origin point
		#       for the space while the following points gives a spanning set of vectors wrt the origin point
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		# Defines Model
		env = gp.Env(empty=True)
		#env.setParam('TokenServer', 'grb-ts.ufhpc')
		env.start()
		model = gp.Model("warmup_model_" + self.model_name, env=env)
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb,
		                           ub=ub)
		model.addConstr(S @ react_flux == b, name="c")
		if max_combined_flux > 0:
			obj = -1 * np.ones((1, np.size(ub)))
			model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
			model.optimize()
			min_flux = model.objVal * -1
			max_allowed_flux = rel_distance_to_min_flux * min_flux
			S = np.vstack((S, np.ones((1, np.shape(S)[1]))))
			S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
			lb = np.hstack((lb, 0))
			ub = np.hstack((ub, max_allowed_flux))
			b = np.hstack((b, 0))
			S[-1, -1] = -1
			rxn_list.append("total_flux")
			model = gp.Model("warmup_model_" + self.model_name, env=env)
			react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb, ub=ub)
			# print(react_flux)
			model.addConstr(S @ react_flux == b, name="c")
		
		# Suppress printing to console
		model.Params.LogToConsole = 0
		
		# Initializes number of subsequent vectors which have been rejected
		search_number = 0
		
		# Initializes rank
		point_matrix = np.zeros((n,len(lb)))
		for i in range(n):
			if i % 100 == 0:
				print(i)
			# Finds random vertex
			while model.Status == 1:
				obj = (2 * np.random.random_sample(S.shape[1]) - 1)
				model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
				model.optimize()
			point_matrix[i] = react_flux.X
			model.reset()
		np.save(save_path / ("warmup_" + self.model_name), point_matrix)
		np.save(save_path / ("warmup_header_" + self.model_name), rxn_list)
		return point_matrix
		
	
	def HRSampler(self, origin_and_warmup_points, number_of_points, stepsPerPoint, rxn_penalty_vec=None, **kwargs):
		# Description:
		#   Performs the coordinate direction hit and run algorithm
		# Inputs:
		#   origin_and_warmup_points: 2d numpy array which defines the spanning basis for the space, each row provides a
		#       point in the space. With the first row giving the origin point and the following rows giving points
		#       which become a spanning set of vectors wrt the point
		#   number_of_points: number of points the algorithm should try to generate
		#   steps_perPoint: number of points to skip between saved points
		#   artificial_centering: If True using artificial centering hit and run algorithm
		#   kwargs:
		#       objective_constraint: Dictionary to add objective constraint
		# Outputs:
		#   returns set of randomly selected points
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		maxMinTol = 1e-9
		uTol = 1e-9
		dTol = 1e-9
		Sv_cont_tol = 1e-9
		prev = 0
		
		# find center
		centerPoint = np.mean(origin_and_warmup_points, axis=0)
		print("centerPoint")
		# flux_balance_test(centerPoint,S,b,lb,ub,c)
		
		warmup_points = origin_and_warmup_points[1:]
		# origin = origin_and_warmup_points[0]
		
		# Making the origin the center of the basis avoids starting it on a boundary
		origin = centerPoint
		
		if "maxTime" in kwargs.keys():
			maxTime = kwargs["maxTime"]
		else:
			maxTime = 1000 * 3600
		
		if "initial_point" in kwargs.keys():
			prev_point = kwargs["initial_point"]
		else:
			prev_point = origin

		if self.total_flux_limits is not None:
			if rxn_penalty_vec is None:
				rxn_penalty_vec = np.ones((1, np.shape(S)[1]))
			else:
				rxn_penalty_vec = np.reshape(rxn_penalty_vec, (1, -1))
			S = np.vstack((S, rxn_penalty_vec))
			S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
			lb = np.hstack((lb, min(self.total_flux_limits)))
			ub = np.hstack((ub, max(self.total_flux_limits)))
			b = np.hstack((b, 0))
			S[-1, -1] = -1
			rxn_list.append("total_flux")
			
		
		N = sc.linalg.null_space(S)
		
		total_steps = 0
		total_steps_needed = number_of_points * stepsPerPoint
		
		# random_vector = np.random.choice(total_steps_needed)
		t0 = time.time()
		
		points = np.zeros((number_of_points, len(origin)))
		point_count = 1
		
		while point_count <= number_of_points:
			# draws several random numbers at a time for improved performance
			randVector = np.random.random_sample(stepsPerPoint)
			
			step_count = 1
			while step_count <= stepsPerPoint:
				#if step_count %50 ==0 :
				#	print(point_count,step_count)
				# print(step_count)
				# Pick a random warmup point
				random_row_index = np.random.choice(np.shape(warmup_points)[0])
				rand_point = warmup_points[random_row_index]
				
				# Get a direction from the basis point to the warmup point
				u = self.normalize(rand_point - origin)
				
				# Figure out the distances to upper and lower bounds
				distUb = (ub - prev_point)
				distLb = (prev_point - lb)
				
				# Figure out if we are too close to a boundary
				#validDir = ((distUb > dTol) & (distLb > dTol))
				
				# Finds list of entries of unit vector which are either positive or negative
				# and which fall within tolerance values
				posValidDir = ((distUb > dTol) & (distLb > dTol)) & (u > uTol)
				negValidDir = ((distUb > dTol) & (distLb > dTol)) & (u < (-1 * uTol))
				
				# Find the largest allowed step size in both the positive and negative directions
				# How many lengths of u can you move in any valid direction before you hit a wall
				pos_u_rel_dis_ub = np.compress(posValidDir, distUb, axis=0) / np.compress(posValidDir, u, axis=0)
				neg_u_rel_dis_ub = np.compress(negValidDir, distUb, axis=0) / np.compress(negValidDir, u, axis=0)
				
				# additional upperbound constraint to prevent exceeding imposed flux limit
				# Cases:
				# sum u is negative
				# sum u is positive
				# sum u is small
				#if max_combined_flux != -1:
				#	current_flux = np.sum(prev_point)
				#	max_add_from_u = max_combined_flux - current_flux

					# This might look weird, but we are just seeing how many unit vectors we can move before
					# we hit the flux limit, if that number is positive and if it is smaller than
					# the first positive distance to an upper bound it replaces it as it is more restrictive.
					# later on the code will check to see if it is more restrictive than the other ub distances
					
					# If we have to move a negative amount in the direction of the unit vector the code
					# also just checks to see if that is more restrictive than the first negative distance to an
					# upper bound
				#	rel_u_to_add = max_add_from_u/np.sum(u)
				#	if rel_u_to_add > 0 and rel_u_to_add < pos_u_rel_dis_ub[0]:
				#		pos_u_rel_dis_ub[0] = rel_u_to_add
				#	elif rel_u_to_add < 0 and rel_u_to_add > neg_u_rel_dis_ub[0]:
				#		neg_u_rel_dis_ub[0] = rel_u_to_add
				
				# Find the smallest allowed step size in both the positive and negative directions
				neg_u_rel_dis_lb  = -1 * np.compress(negValidDir, distLb, axis=0) / np.compress(negValidDir, u, axis=0)
				pos_u_rel_dis_lb = -1 * np.compress(posValidDir, distLb, axis=0) / np.compress(posValidDir, u, axis=0)
				
				# find the current value of the point
			
				
				# Find the smallest and largest allowed step sizes
				maxStep = np.min(np.hstack((pos_u_rel_dis_ub, neg_u_rel_dis_lb)))
				minStep = np.max(np.hstack((pos_u_rel_dis_lb, neg_u_rel_dis_ub)))
				
				
				if ((abs(minStep) < maxMinTol) & (abs(maxStep) < maxMinTol)) or (minStep > maxStep):
					print(f"Warning: \nMin step = {minStep}\n Max step = {maxStep}")
					continue
				
				# grab random value from pre-generated list for the step distance
				stepDist = randVector[step_count - 1] * (maxStep - minStep) + minStep
				curPoint = prev_point + stepDist * u
				
				

				if total_steps % 100 == 0:
					curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				
				if total_steps % 200000 == 0 and False:
					print(f"Error {max(curPoint - ub)}, {max(lb - curPoint)}")
				# Both values should be negative
				
				time_elapsed = time.time() - t0
				
				if step_count % 500000 == 0:
					time_per_step = time_elapsed / step_count
					print(
						f"{point_count}, {step_count}, {time_elapsed / 60}, {(total_steps_needed - step_count) * time_per_step / 60}")
				
				if np.any((ub - curPoint) < 0) or np.any((curPoint - lb) < 0):
					indexer = np.arange(0, np.size(curPoint))
					overInd = indexer[curPoint > ub]
					underInd = indexer[curPoint < lb]
					#print(overInd)
					#print(underInd)
					# scaling might be useful here
					for j in overInd:
						curPoint[j] = ub[j]
					for j in underInd:
						curPoint[j] = lb[j]
				if total_steps % 2000000 == 0:
					print(f"fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
				prev_point = curPoint
				step_count += 1
				
				total_steps += 1
				if np.floor(total_steps / total_steps_needed * 100) > prev:
					prev = total_steps / total_steps_needed * 100
					print(total_steps / total_steps_needed)
					print(curPoint[-1])
					#print(np.sum(points))
					#print(np.shape(points))
					#print(time_elapsed/(total_steps / total_steps_needed)/60)
				
				if time_elapsed > maxTime:
					points[point_count] = curPoint
					print(f"Time constraint reached after {point_count} points")
					return points
			points[point_count - 1] = curPoint
			
			point_count += 1
		
		return points
	
	def HRSampler_therm(self, origin_and_warmup_points, number_of_points, stepsPerPoint,skip_therm, NS,simp_neg_ind,comp_neg_ind,comp_perm,cutter,exchange_cutter=None, **kwargs):
		# Description:
		#   Performs the coordinate direction hit and run algorithm
		# Inputs:
		#   origin_and_warmup_points: 2d numpy array which defines the spanning basis for the space, each row provides a
		#       point in the space. With the first row giving the origin point and the following rows giving points
		#       which become a spanning set of vectors wrt the point
		#   number_of_points: number of points the algorithm should try to generate
		#   steps_perPoint: number of points to skip between saved points
		#   artificial_centering: If True using artificial centering hit and run algorithm
		#   kwargs:
		#       objective_constraint: Dictionary to add objective constraint
		# Outputs:
		#   returns set of randomly selected points
		lb, ub, S, b, rxn_list, met_list = self.dicts_to_mats()
		maxMinTol = 1e-9
		uTol = 1e-9
		dTol = 1e-9
		Sv_cont_tol = 1e-9
		prev = 0
		
		# find center
		centerPoint = np.mean(origin_and_warmup_points, axis=0)
		print("centerPoint")
		# flux_balance_test(centerPoint,S,b,lb,ub,c)
		
		warmup_points = origin_and_warmup_points[1:]
		# origin = origin_and_warmup_points[0]
		
		# Making the origin the center of the basis avoids starting it on a boundary
		origin = centerPoint
		
		if "maxTime" in kwargs.keys():
			maxTime = kwargs["maxTime"]
		else:
			maxTime = 1000 * 3600
		
		if "initial_point" in kwargs.keys():
			prev_point = kwargs["initial_point"]
		else:
			prev_point = origin
		
		if self.total_flux_limits is not None:
			S = np.vstack((S, np.ones((1, np.shape(S)[1]))))
			S = np.hstack((S, np.zeros((np.shape(S)[0], 1))))
			lb = np.hstack((lb, min(self.total_flux_limits)))
			ub = np.hstack((ub, max(self.total_flux_limits)))
			b = np.hstack((b, 0))
			S[-1, -1] = -1
			rxn_list.append("total_flux")
		
		N = sc.linalg.null_space(S)
		
		total_steps = 0
		total_steps_needed = number_of_points * stepsPerPoint
		
		# random_vector = np.random.choice(total_steps_needed)
		t0 = time.time()
		
		points = np.zeros((number_of_points, len(origin)))
		point_count = 1
		stepsPerPoint_cp = copy.deepcopy(stepsPerPoint)
		
		first_pass = True
		
		while point_count <= number_of_points:
			if first_pass:
				stepsPerPoint = stepsPerPoint_cp
			else:
				stepsPerPoint = skip_therm
			# draws several random numbers at a time for improved performance
			randVector = np.random.random_sample(stepsPerPoint)
			
			step_count = 1
			while step_count <= stepsPerPoint:
				# if step_count %50 ==0 :
				#	print(point_count,step_count)
				# print(step_count)
				# Pick a random warmup point
				random_row_index = np.random.choice(np.shape(warmup_points)[0])
				rand_point = warmup_points[random_row_index]
				
				# Get a direction from the basis point to the warmup point
				u = self.normalize(rand_point - origin)
				
				# Figure out the distances to upper and lower bounds
				distUb = (ub - prev_point)
				distLb = (prev_point - lb)
				
				# Figure out if we are too close to a boundary
				# validDir = ((distUb > dTol) & (distLb > dTol))
				
				# Finds list of entries of unit vector which are either positive or negative
				# and which fall within tolerance values
				posValidDir = ((distUb > dTol) & (distLb > dTol)) & (u > uTol)
				negValidDir = ((distUb > dTol) & (distLb > dTol)) & (u < (-1 * uTol))
				
				# Find the largest allowed step size in both the positive and negative directions
				# How many lengths of u can you move in any valid direction before you hit a wall
				pos_u_rel_dis_ub = np.compress(posValidDir, distUb, axis=0) / np.compress(posValidDir, u, axis=0)
				neg_u_rel_dis_ub = np.compress(negValidDir, distUb, axis=0) / np.compress(negValidDir, u, axis=0)
				
				# additional upperbound constraint to prevent exceeding imposed flux limit
				# Cases:
				# sum u is negative
				# sum u is positive
				# sum u is small
				# if max_combined_flux != -1:
				#	current_flux = np.sum(prev_point)
				#	max_add_from_u = max_combined_flux - current_flux
				
				# This might look weird, but we are just seeing how many unit vectors we can move before
				# we hit the flux limit, if that number is positive and if it is smaller than
				# the first positive distance to an upper bound it replaces it as it is more restrictive.
				# later on the code will check to see if it is more restrictive than the other ub distances
				
				# If we have to move a negative amount in the direction of the unit vector the code
				# also just checks to see if that is more restrictive than the first negative distance to an
				# upper bound
				#	rel_u_to_add = max_add_from_u/np.sum(u)
				#	if rel_u_to_add > 0 and rel_u_to_add < pos_u_rel_dis_ub[0]:
				#		pos_u_rel_dis_ub[0] = rel_u_to_add
				#	elif rel_u_to_add < 0 and rel_u_to_add > neg_u_rel_dis_ub[0]:
				#		neg_u_rel_dis_ub[0] = rel_u_to_add
				
				# Find the smallest allowed step size in both the positive and negative directions
				neg_u_rel_dis_lb = -1 * np.compress(negValidDir, distLb, axis=0) / np.compress(negValidDir, u, axis=0)
				pos_u_rel_dis_lb = -1 * np.compress(posValidDir, distLb, axis=0) / np.compress(posValidDir, u, axis=0)
				
				# find the current value of the point
				
				# Find the smallest and largest allowed step sizes
				maxStep = np.min(np.hstack((pos_u_rel_dis_ub, neg_u_rel_dis_lb)))
				minStep = np.max(np.hstack((pos_u_rel_dis_lb, neg_u_rel_dis_ub)))
				
				if ((abs(minStep) < maxMinTol) & (abs(maxStep) < maxMinTol)) or (minStep > maxStep):
					print(f"Warning: \nMin step = {minStep}\n Max step = {maxStep}")
					continue
				
				# grab random value from pre-generated list for the step distance
				stepDist = randVector[step_count - 1] * (maxStep - minStep) + minStep
				curPoint = prev_point + stepDist * u
				
				if total_steps % 100 == 0:
					curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				
				if total_steps % 200000 == 0 and False:
					print(f"Error {max(curPoint - ub)}, {max(lb - curPoint)}")
				# Both values should be negative
				
				time_elapsed = time.time() - t0
				
				if step_count % 500000 == 0:
					time_per_step = time_elapsed / step_count
					print(
						f"{point_count}, {step_count}, {time_elapsed / 60}, {(total_steps_needed - step_count) * time_per_step / 60}")
				
				if np.any((ub - curPoint) < 0) or np.any((curPoint - lb) < 0):
					indexer = np.arange(0, np.size(curPoint))
					overInd = indexer[curPoint > ub]
					underInd = indexer[curPoint < lb]
					# print(overInd)
					# print(underInd)
					# scaling might be useful here
					for j in overInd:
						curPoint[j] = ub[j]
					for j in underInd:
						curPoint[j] = lb[j]
				if total_steps % 2000000 == 0:
					print(f"fidErr {np.max(np.abs(np.matmul(S, curPoint)))}")
				prev_point = curPoint
				step_count += 1
				
				total_steps += 1
				if np.floor(total_steps / total_steps_needed * 100) > prev:
					prev = total_steps / total_steps_needed * 100
					print(total_steps / total_steps_needed)
				# print(np.sum(points))
				# print(np.shape(points))
				# print(time_elapsed/(total_steps / total_steps_needed)/60)
				
				if time_elapsed > maxTime:
					points[point_count] = curPoint
					print(f"Time constraint reached after {point_count} points")
					return points
			int_point = self.positive_flux_point_to_bi_dir(curPoint[:-1], simp_neg_ind, comp_neg_ind, comp_perm, cutter,
		                                                          exchange_cutter=exchange_cutter)
			if self.generate_therm_model(NS,int_point) == 2:
				#print("accepted point")
				points[point_count - 1] = curPoint
				point_count += 1
				first_pass = True
			else:
				#print("rejected point")
				first_pass = False
			
		
		return points
	
	def precomp_positive_flux_point_to_bi_dir(self,pos_rxn_list,prune_exchange = True, prune_specific = None):
		# pos_rxn_list is used in stead of self.rxn_dict.keys() so that a desired order of keys may be passed
		if prune_specific is None:
			prune_specific = []
		elif not prune_exchange:
			print("prune_exchange is False, but prune specific contains entries, this is a contradiction")
		
		if np.alltrue(np.array(["_inverted" not in i for i in pos_rxn_list])):
			print("WARNING: NO INVERTED RXNS FOUND, WRONG MODEL MAY BE IN USE")
		zeros = np.zeros(len(pos_rxn_list), dtype=float)
		ones = np.ones(len(pos_rxn_list), dtype=float)
		
		
		
		cutter = np.nonzero(np.where(
			np.array([("_inverted" in i and i.replace("_inverted", "") in pos_rxn_list) for i in pos_rxn_list]),
			zeros, ones))[0]

		simp_neg_ind = np.where(
			np.array([("_inverted" in i and i.replace("_inverted", "") not in pos_rxn_list) for i in pos_rxn_list]),
			-1 * ones, ones)
		comp_neg_ind = []
		comp_perm = []
		for i in range(len(pos_rxn_list)):
			# slowly read through to make sure this is doing what it is supposed to do
			rxn_name = pos_rxn_list[i]
			# if reaction is forward but its inverse is in list
			# make the comp list 0, but tell the comp perm list what the index of the inverted (-1) reaction is
			# in that way if you do comp_neg_ind[comp_perm] the comp perm will request the -1 from the comp neg ind to come to the
			# location of its positive value. Furthermore if you modify comp_neg_ind to have -1* its flux and run comp_neg_ind[comp_perm]
			# you will have the negative fluxes in the right place to simply add the list and apply the cutter
			if "_inverted" not in rxn_name and rxn_name + "_inverted" in pos_rxn_list:
				comp_neg_ind.append(0.0)
				comp_perm.append(pos_rxn_list.index(rxn_name + "_inverted"))
			elif "_inverted" in rxn_name and rxn_name.replace("_inverted", "") in pos_rxn_list:
				comp_neg_ind.append(-1.0)
				comp_perm.append(-5)
			else:
				comp_neg_ind.append(0.0)
				comp_perm.append(-5)
		for i in range(len(comp_perm)):
			if comp_perm[i] == -5:
				comp_perm[i] = np.argmax(comp_neg_ind)
		comp_perm = np.array(comp_perm)
		comp_neg_ind = np.array(comp_neg_ind)
		pos_rxn_array = np.array(pos_rxn_list)
		cut_rxn_array = pos_rxn_array[cutter]
		if prune_exchange:
			exchange_cutter = np.nonzero(np.where(
				np.array([(i in prune_specific or len(self.rxn_dict[i]["S"].keys()) == 1) for i in cut_rxn_array]),
				zeros[cutter], ones[cutter]))[0]
		else:
			exchange_cutter = None

		
		pos_rxn_array_post_cut = pos_rxn_array[cutter]
		return simp_neg_ind, comp_neg_ind, comp_perm, cutter, exchange_cutter
	
	def positive_flux_point_to_bi_dir(self,point,simp_neg_ind,comp_neg_ind,comp_perm,cutter,exchange_cutter = None):
		bi_point = copy.deepcopy(point)
		comp_neg_ind_cp = copy.deepcopy(comp_neg_ind)
		bi_point *= simp_neg_ind
		comp_neg_ind_cp *= bi_point
		comp_neg_ind_cp = comp_neg_ind_cp[comp_perm]
		bi_point += comp_neg_ind_cp
		nbi_point = bi_point[cutter]
		if exchange_cutter is not None:
			nbi_point = nbi_point[exchange_cutter]
		return nbi_point
	
	def positive_S_to_bi_dir(self, S, simp_neg_ind, cutter, exchange_cutter=None):
		bi_dir_S = copy.deepcopy(S)
		for col_ind in range(np.shape(bi_dir_S)[1]):
			bi_dir_S[:,col_ind] *= simp_neg_ind[col_ind]
		bi_dir_S = bi_dir_S[:,cutter]
		if exchange_cutter is not None:
			bi_dir_S = bi_dir_S[:,exchange_cutter]
		return bi_dir_S
	
	def convert_model_to_positive_flux(self):
		# Only works for getting S
		reaction_names = list(self.rxn_dict.keys())
		for rxn_name in reaction_names:
			if self.rxn_dict[rxn_name]['ub'] <= 0:
				new_rxn = rxn_name+"_inverted"
				nrxn_dict = {}
				for key in self.rxn_dict[rxn_name].keys():
					if key == "S":
						nrxn_dict["S"] = {i:-1*self.rxn_dict[rxn_name]['S'][i] for i in self.rxn_dict[rxn_name]['S'].keys()}
					elif key == "lb":
						nrxn_dict["ub"] = self.rxn_dict[rxn_name]["lb"]*-1
					elif key == "ub":
						nrxn_dict["lb"] = self.rxn_dict[rxn_name]["ub"]*-1
					else:
						nrxn_dict[key] = self.rxn_dict[rxn_name][key]
				del self.rxn_dict[rxn_name]
				self.rxn_dict[new_rxn] = copy.deepcopy(nrxn_dict)
			elif self.rxn_dict[rxn_name]['lb'] < 0:
				new_rxn = rxn_name + "_inverted"
				nrxn_dict = {}
				for key in self.rxn_dict[rxn_name].keys():
					if key == "S":
						nrxn_dict["S"] = {i: -1 * self.rxn_dict[rxn_name]['S'][i] for i in
						                  self.rxn_dict[rxn_name]['S'].keys()}
					elif key == "lb":
						nrxn_dict["ub"] = self.rxn_dict[rxn_name]["lb"] * -1
					elif key == "ub":
						nrxn_dict["lb"] = 0
					else:
						nrxn_dict[key] = self.rxn_dict[rxn_name][key]
				self.rxn_dict[new_rxn] = copy.deepcopy(nrxn_dict)
				self.update_reaction_bounds(rxn_name,0,'keep')
	
	def convert_model_to_bidirectional_flux(self):
		# This currently should only be used for S, its not fully tested for other entries
		rxn_names = list(self.rxn_dict.keys())
		for rxn_name in rxn_names:
			if "_inverted" in rxn_name:
				forward_rxn_name = rxn_name.replace("_inverted","")
				if forward_rxn_name in rxn_names:
					if self.rxn_dict[forward_rxn_name]["ub"] == self.rxn_dict[forward_rxn_name]["lb"] and self.rxn_dict[forward_rxn_name]["lb"] == 0:
						self.update_reaction_bounds(forward_rxn_name, -1 * self.rxn_dict[rxn_name]["ub"], -1 * self.rxn_dict[rxn_name]["lb"])
					elif self.rxn_dict[rxn_name]["ub"] == self.rxn_dict[rxn_name]["lb"] and self.rxn_dict[rxn_name]["lb"] == 0:
						pass
					else:
						self.update_reaction_bounds(forward_rxn_name,-1*self.rxn_dict[rxn_name]["ub"],"keep")
					del self.rxn_dict[rxn_name]
				else:
					self.rxn_dict[forward_rxn_name] = self.rxn_dict.pop(rxn_name)
					lb_ph = copy.deepcopy(self.rxn_dict[forward_rxn_name]["lb"])
					ub_ph = copy.deepcopy(self.rxn_dict[forward_rxn_name]["ub"])
					rxn_S = copy.deepcopy(self.rxn_dict[forward_rxn_name]['S'])
					self.rxn_dict[forward_rxn_name]["S"] = {i: -1 * rxn_S[i] for i in rxn_S.keys()}
					self.update_reaction_bounds(forward_rxn_name,-1*ub_ph,-1*lb_ph)

	
	
	def weight_fluxes(self, eq_bound_tol=0.0, new_bound_range_tol=10):
		max_constrained_flux = -1
		max_constrained_flux_id = -1
		for i in range(len(self.reaction_names)):
			if self.lb[i] >= self.ub[i] - eq_bound_tol:
				if abs(self.lb[i]) > max_constrained_flux:
					max_constrained_flux = abs(self.lb[i])
					max_constrained_flux_id = i
		print(f"New extreme bound is {max_constrained_flux * new_bound_range_tol}")
		print(f"Set based on {self.reaction_names[max_constrained_flux_id]}")
		for i in range(len(self.reaction_names)):
			if self.lb[i] < -1 * new_bound_range_tol * max_constrained_flux:
				self.lb[i] = -1 * new_bound_range_tol * max_constrained_flux
			elif self.ub[i] < 1 * new_bound_range_tol * max_constrained_flux:
				self.ub[i] = 1 * new_bound_range_tol * max_constrained_flux
	
	def find_essential_fluxes_fast(self, list_to_return="index", print_progress=False, **kwargs):
		# A more efficient system to find the essential fluxes
		# The old method is kept for comparison
		
		lb,ub,S,b,rxn_list,met_list = self.dicts_to_mats()
		essential_reactions = []
		true_lb = copy.deepcopy(lb)
		true_ub = copy.deepcopy(ub)
		
		true_lb_2 = copy.deepcopy(lb)
		true_ub_2 = copy.deepcopy(ub)
		if self.env != None:
			model = gp.Model("fva_model_" + self.model_name, env=self.env)
		else:
			model = gp.Model("fva_model_" + self.model_name)
		react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux",
		                           lb=copy.deepcopy(lb),
		                           ub=copy.deepcopy(ub))
		model.addConstr(S @ react_flux == b, name="c")
		# obj = np.zeros((1, self.S.shape[1]))
		
		# model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.Params.LogToConsole = 0
		model.Params.ScaleFlag = 0
		model.Params.NumericFocus = 3
		model.Params.Method = 1
		# model.Params.Presolve = 0
		# model.Params.Aggregate = 0
		model.Params.FeasibilityTol = 1e-9
		
		index_size = np.size(lb)
		index_select = [i for i in range(index_size)]
		if "seed" in kwargs.keys():
			random.Random(kwargs["seed"]).shuffle(index_select)
		else:
			random.shuffle(index_select)
		for i in range(index_size):
			if print_progress and i % 100 == 0:
				print((index_size - i) / index_size)
				print(len(essential_reactions))
			reaction_index_to_check = index_select[i]
			if lb[reaction_index_to_check] > 0 or ub[reaction_index_to_check] < 0:
				# This additional line prevents an issue where a preset bound is broken but the system is still feasible
				# If a metabolite intake can range between -100 and 100 and we force it to be between 9 and 10
				# setting it to 0 and testing for feasibility is not enough
				is_feasible = False
			else:
				save_lb = true_lb_2[reaction_index_to_check]
				save_ub = true_ub_2[reaction_index_to_check]
				true_lb[reaction_index_to_check] = 0
				true_ub[reaction_index_to_check] = 0
				react_flux.lb = true_lb
				react_flux.ub = true_ub
				model.optimize()
				is_feasible = (model.Status == 2)
				if not is_feasible:
					if model.Status != 3:
						print(model.Status)
						input("Error, unexpected model status")
				model.reset()
			react_flux.lb = true_lb_2
			react_flux.ub = true_ub_2
			if is_feasible:
				continue
			else:
				if not (lb[reaction_index_to_check] > 0 or ub[reaction_index_to_check] < 0):
					true_lb[reaction_index_to_check] = save_lb
					true_ub[reaction_index_to_check] = save_ub
				if list_to_return == "index":
					essential_reactions.append(reaction_index_to_check)
				elif list_to_return == "names":
					essential_reactions.append(reaction_names[reaction_index_to_check])
		if "dict_count" in kwargs.keys():
			initial_flux_counts_copy = copy.deepcopy(kwargs["dict_count"])
			for i in essential_reactions:
				initial_flux_counts_copy[i] += 1
			return initial_flux_counts_copy
		else:
			print(len(essential_reactions))
			return essential_reactions
	
	def fit_to_experimental_data(self, experimental_data_file, alpha_array, assign_opt_alpha_fluxes=True, tol=1e-7,
	                             search_save_path=""):
		# experiment vector = ei
		# model vector = vi
		# (vi-aei)^2/(aei)^2 = vi^2/(aei)^2-2vi/aei+1
		# returns:
		# model with best experimental fit for given alpha
		# the relative error (objective term)
		# the index values of the experimental data values
		lb,ub,S,b,ordered_rxn,ordered_met = self.dicts_to_mats()
		experimental_fluxes = np.loadtxt(experimental_data_file, delimiter=",", dtype=str)
		model_index_of_experimental_fluxes = []
		e = np.zeros(len(lb))
		amount_to_add_back = 0
		for i in range(np.shape(experimental_fluxes)[0]):
			if i > 0:
				rxn_name = experimental_fluxes[i, 0]
				flux_values = experimental_fluxes[i, 1].astype('float64')
				index_of_reaction = ordered_rxn.index(rxn_name)
				model_index_of_experimental_fluxes.append(index_of_reaction)
				e[index_of_reaction] = flux_values
				amount_to_add_back += 1
		# For very low values of experimental flux gurobi will run into accuracy errors
		# The current solution is to just take the upper limit of the experimental value of the one
		# poorly behaving flux, but for other cases it might be worth establishing a minimal cutoff point
		# And simply do linear minimaztion on the small values error
		array_of_states = np.zeros((np.size(alpha_array), np.size(lb)))
		array_of_obj_val = np.zeros((np.size(alpha_array)))
		
		min_ind = 0
		min_objective = -1
		for alpha_ind in range(np.size(alpha_array)):
			Q = np.zeros((len(e), len(e)))
			for i in range(len(e)):
				if e[i] != 0:
					Q[i, i] = 1 / (alpha_array[alpha_ind] * e[i]) ** 2
				else:
					Q[i, i] = 0
			c = np.zeros(len(e))
			for i in range(len(e)):
				if e[i] != 0:
					c[i] = 2 / (alpha_array[alpha_ind] * e[i])
				else:
					c[i] = 0
			if self.env != None:
				model = gp.Model("fva_model_" + self.model_name, env=self.env)
			else:
				model = gp.Model("fva_model_" + self.model_name)
			model.Params.NumericFocus = 3
			# model.Params.ScaleFlag = 0
			model.Params.BarConvTol = 0
			model.Params.BarHomogeneous = 1
			# model.Params.LogToConsole = 0
			react_flux = model.addMVar(shape=S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=lb,
			                           ub=ub)
			model.setObjective(react_flux @ Q @ react_flux - c @ react_flux + amount_to_add_back, GRB.MINIMIZE)
			rhs = np.transpose(b)
			model.addConstr(S @ react_flux == rhs, name="c")
			model.optimize()
			if model.Status == 2 or model.Status == 13:
				array_of_states[alpha_ind] = model.X
				array_of_obj_val[alpha_ind] = model.objVal
				if model.objVal < min_objective or min_objective == -1:
					min_objective = model.objVal
					min_ind = alpha_ind
			else:
				array_of_states[alpha_ind] = np.empty_like(model.X) * np.nan
				array_of_obj_val[alpha_ind] = np.nan
		
		# print(np.matmul(np.matmul(model.X, Q), model.X) - np.matmul(c, model.X) + amount_to_add_back)
		if assign_opt_alpha_fluxes:
			for i in model_index_of_experimental_fluxes:
				self.update_reaction_bounds(ordered_rxn[i], array_of_states[min_ind, i] - tol, array_of_states[min_ind, i] + tol)
		if search_save_path != "":
			result_dict = {}
			result_dict["alpha_list"] = alpha_array
			result_dict["saved_rxn_name_order"] = ordered_rxn
			result_dict["saved_met_order"] = ordered_met
			result_dict["best_alpha_index"] = min_ind
			result_dict["model_alignment"] = array_of_states
			result_dict["model_alignment_performances"] = array_of_obj_val
			result_dict["experimental_flux_index"] = model_index_of_experimental_fluxes
			result_dict["experimental_fluxes"] = e
			with open(search_save_path.parent/(search_save_path.stem+".pkl"), "wb") as outfile:
				pickle.dump(result_dict, outfile)
	
	def fit_to_experimental_data_l1_norm(self, experimental_data_file, alpha_array, assign_opt_alpha_fluxes=True, tol=1e-7,
	                             search_save_path=""):
		# experiment vector = ei
		# model vector = vi
		# |vi-aei|
		# (vi-aei)^2/(aei)^2 = vi^2/(aei)^2-2vi/aei+1
		# returns:
		# model with best experimental fit for given alpha
		# the relative error (objective term)
		# the index values of the experimental data values
		
		experimental_fluxes = np.loadtxt(experimental_data_file, delimiter=",", dtype=str)
		model_index_of_experimental_fluxes = []
		e = np.zeros(len(self.lb))
		amount_to_add_back = 0
		for i in range(np.shape(experimental_fluxes)[0]):
			if i > 0:
				rxn_name = experimental_fluxes[i, 0]
				flux_values = experimental_fluxes[i, 1].astype('float64')
				index_of_reaction = self.reaction_names.index(rxn_name)
				model_index_of_experimental_fluxes.append(index_of_reaction)
				e[index_of_reaction] = flux_values
				amount_to_add_back += 1
		# For very low values of experimental flux gurobi will run into accuracy errors
		# The current solution is to just take the upper limit of the experimental value of the one
		# poorly behaving flux, but for other cases it might be worth establishing a minimal cutoff point
		# And simply do linear minimaztion on the small values error
		array_of_states = np.zeros((np.size(alpha_array), np.size(self.lb)))
		array_of_obj_val = np.zeros((np.size(alpha_array)))
		
		min_ind = 0
		min_objective = -1
		for alpha_ind in range(np.size(alpha_array)):
			if self.env != None:
				model = gp.Model("fva_model_" + self.model_name, env=self.env)
			else:
				model = gp.Model("fva_model_" + self.model_name)
			c = np.zeros(len(e))
			for i in range(len(e)):
				if e[i] != 0:
					c[i] = (alpha_array[alpha_ind] * e[i])
				else:
					c[i] = 0
			model.optimize()
			model.Params.NumericFocus = 3
			model.Params.BarConvTol = 0
			model.Params.BarHomogeneous = 1
			print(gp.sys.version)
			react_flux = model.addMVar(shape=self.S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=self.lb,ub=self.ub)

			model.setObjective(b, GRB.MINIMIZE)
			rhs = np.transpose(self.b)
			model.addConstr(self.S @ react_flux == rhs, name="c")
			model.optimize()
			input()
			if model.Status == 2 or model.Status == 13:
				array_of_states[alpha_ind] = model.X
				array_of_obj_val[alpha_ind] = model.objVal
				if model.objVal < min_objective or min_objective == -1:
					min_objective = model.objVal
					min_ind = alpha_ind
			else:
				array_of_states[alpha_ind] = np.empty_like(model.X) * np.nan
				array_of_obj_val[alpha_ind] = np.nan
		
		# print(np.matmul(np.matmul(model.X, Q), model.X) - np.matmul(c, model.X) + amount_to_add_back)
		if assign_opt_alpha_fluxes:
			for i in model_index_of_experimental_fluxes:
				self.update_reaction_bounds(i, array_of_states[min_ind, i] - tol, array_of_states[min_ind, i] + tol)
		if search_save_path != "":
			result_dict = {}
			result_dict["alpha_list"] = alpha_array
			result_dict["best_alpha_index"] = min_ind
			result_dict["model_alignment"] = array_of_states
			result_dict["model_alignment_performances"] = array_of_obj_val
			result_dict["experimental_flux_index"] = model_index_of_experimental_fluxes
			result_dict["experimental_fluxes"] = e
			with open(search_save_path.parent/(search_save_path.stem+".pkl"), "wb") as outfile:
				pickle.dump(result_dict, outfile)
	
	def fit_to_experimental_data_abs_error(self, experimental_data_file, alpha_array, assign_opt_alpha_fluxes=True, tol=1e-7,
	                             search_save_path=""):
		# experiment vector = ei
		# model vector = vi
		# (vi-aei)^2 = vi^2-2vi*ei*a+a^2ei^2
		# returns:
		# model with best experimental fit for given alpha
		# the relative error (objective term)
		# the index values of the experimental data values
		
		experimental_fluxes = np.loadtxt(experimental_data_file, delimiter=",", dtype=str)
		model_index_of_experimental_fluxes = []
		e = np.zeros(len(self.lb))
		for i in range(np.shape(experimental_fluxes)[0]):
			if i > 0:
				rxn_name = experimental_fluxes[i, 0]
				flux_values = experimental_fluxes[i, 1].astype('float64')
				index_of_reaction = self.reaction_names.index(rxn_name)
				model_index_of_experimental_fluxes.append(index_of_reaction)
				e[index_of_reaction] = flux_values
		# For very low values of experimental flux gurobi will run into accuracy errors
		# The current solution is to just take the upper limit of the experimental value of the one
		# poorly behaving flux, but for other cases it might be worth establishing a minimal cutoff point
		# And simply do linear minimaztion on the small values error
		array_of_states = np.zeros((np.size(alpha_array), np.size(self.lb)))
		array_of_obj_val = np.zeros((np.size(alpha_array)))
		
		min_ind = 0

		min_objective = -1
		for alpha_ind in range(np.size(alpha_array)):
			Q = np.zeros((len(e), len(e)))
			for i in range(len(e)):
				if e[i] != 0:
					Q[i, i] = 1
				else:
					Q[i, i] = 0
			c = np.zeros(len(e))
			for i in range(len(e)):
				if e[i] != 0:
					c[i] = 2*e[i]*alpha_array[alpha_ind]
				else:
					c[i] = 0
			if self.env != None:
				model = gp.Model("fva_model_" + self.model_name, env=self.env)
			else:
				model = gp.Model("fva_model_" + self.model_name)
			model.Params.NumericFocus = 3
			# model.Params.ScaleFlag = 0
			model.Params.BarConvTol = 0
			model.Params.BarHomogeneous = 1
			# model.Params.LogToConsole = 0
			react_flux = model.addMVar(shape=self.S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=self.lb,
			                           ub=self.ub)
			const_add = np.sum((c/2)**2)
			model.setObjective(react_flux @ Q @ react_flux - c @ react_flux + const_add, GRB.MINIMIZE)
			rhs = np.transpose(self.b)
			model.addConstr(self.S @ react_flux == rhs, name="c")
			model.optimize()
			if model.Status == 2 or model.Status == 13:
				array_of_states[alpha_ind] = model.X
				array_of_obj_val[alpha_ind] = model.objVal
				if model.objVal < min_objective or min_objective == -1:
					min_objective = model.objVal
					min_ind = alpha_ind
			else:
				array_of_states[alpha_ind] = np.empty_like(model.X) * np.nan
				array_of_obj_val[alpha_ind] = np.nan
		# print(np.matmul(np.matmul(model.X, Q), model.X) - np.matmul(c, model.X) + amount_to_add_back)
		if assign_opt_alpha_fluxes:
			for i in model_index_of_experimental_fluxes:
				self.update_reaction_bounds(i, array_of_states[min_ind, i] - tol, array_of_states[min_ind, i] + tol)
		if search_save_path != "":
			result_dict = {}
			result_dict["alpha_list"] = alpha_array
			result_dict["best_alpha_index"] = min_ind
			result_dict["model_alignment"] = array_of_states
			result_dict["model_alignment_performances"] = array_of_obj_val
			result_dict["experimental_flux_index"] = model_index_of_experimental_fluxes
			result_dict["experimental_fluxes"] = e
			with open(search_save_path.parent/(search_save_path.stem+".pkl"), "wb") as outfile:
				pickle.dump(result_dict, outfile)
	
	def pinch_restricted_exchange_reactions(self, allowed_exchange_reaction_file, restore_essential=True,true_exchange_signifier = "",restore_previously_pinched = False,restore_value = 1000):
		# Leave option to keep excrete free
		met_rxn_dict = self.get_exchanged_reaction_info(true_exchange_signifier = true_exchange_signifier)
		reactions_to_allow = list(np.loadtxt(allowed_exchange_reaction_file, dtype="str", delimiter=",")[1:, 0])
		reactions_to_allow = [i.replace(" ", "") for i in reactions_to_allow]
		exchange_type = list(np.loadtxt(allowed_exchange_reaction_file, dtype="str", delimiter=",")[1:, 1])
		exchange_type = [i.replace(" ", "") for i in exchange_type]
		for rxn in met_rxn_dict.keys():
			if rxn not in reactions_to_allow:
				# right now this still allows excretion
				if restore_essential:
					ub_save = self.rxn_dict[rxn]['ub']
					lb_save = self.rxn_dict[rxn]['ub']
				self.update_reaction_bounds(rxn, 0, 0)
				if restore_essential:
					if not self.test_feasibility():
						print(f"removing reaction {rxn} breaks feasibility")
						self.update_reaction_bounds(rxn, lb_save, ub_save)
			elif rxn in reactions_to_allow:
				rxn_index = reactions_to_allow.index(rxn)
				rxn_exchange_type = exchange_type[rxn_index]
				if rxn_exchange_type == "intake":
					if restore_previously_pinched:
						self.update_reaction_bounds(rxn, -1*restore_value, 0)
					else:
						self.update_reaction_bounds(rxn, "keep", 0)
				elif rxn_exchange_type == "excrete":
					if restore_previously_pinched:
						self.update_reaction_bounds(rxn, 0, restore_value)
					else:
						self.update_reaction_bounds(rxn, 0, "keep")
				elif rxn_exchange_type == "both" and restore_previously_pinched:
					self.update_reaction_bounds(rxn, -1 * restore_value, restore_value)
				elif rxn_exchange_type != "both":
					print("Exchange Reaction direction not understood")

	
	def formula_to_dict(self, formula):
		chemical_dict = {}
		initial_list = [i for i in range(len(formula)) if formula[i].isupper()]
		len_list = [formula[initial_list[i]:initial_list[i + 1]] for i in range(len(initial_list) - 1)]
		len_list.append(formula[initial_list[-1]:])
		for i in len_list:
			first_number_location = 0
			for alnum in range(len(i)):
				if i[alnum].isalpha():
					first_number_location += 1
			element = i[:first_number_location]
			if i[first_number_location:] != "":
				number = int(i[first_number_location:])
			else:
				number = 1
			if element in chemical_dict.keys():
				chemical_dict[element] += number
			else:
				chemical_dict[element] = number
		return chemical_dict
	
	def print_nonpinched_exchange(self):
		exch_dict = self.get_exchanged_reaction_info()
		
		for e_rxn in exch_dict.keys():
			if not (self.rxn_dict[e_rxn]["lb"] == 0 and self.rxn_dict[e_rxn]["ub"] == 0):
				if self.rxn_dict[e_rxn]["lb"] == 0:
					print(e_rxn, "excrete")
				elif self.rxn_dict[e_rxn]["ub"] == 0:
					print(e_rxn, "intake")
				else:
					print(e_rxn, "both")
	
	def test_essential_list(self, reaction_index_list, use_gp_env=True):
		# test against set too
		lb, ub, S, b, ordered_rxn, ordered_met = self.dicts_to_mats()
		if isinstance(reaction_index_list[0], str):
			reaction_index_list = [ordered_rxn.index(name) for name in reaction_index_list]
		if len(set(reaction_index_list)) != len(reaction_index_list):
			return False, "Duplicate Essential Flux"
		save_ub_full = copy.deepcopy(ub)
		save_lb_full = copy.deepcopy(lb)
		essential_flux_list = np.zeros(len(ordered_rxn))
		for rxn_ind in reaction_index_list:
			essential_flux_list[rxn_ind] = 1
		for essential_matrix_index in range(len(essential_flux_list)):
			# if the reaction is essential then lb and ub will be the same, if it is not then lb and ub will be 0
			self.update_reaction_bounds(ordered_rxn[essential_matrix_index],
			                            essential_flux_list[essential_matrix_index] * lb[essential_matrix_index],
			                            essential_flux_list[essential_matrix_index] * ub[essential_matrix_index])
			lb[essential_matrix_index] = essential_flux_list[essential_matrix_index] * lb[essential_matrix_index]
			ub[essential_matrix_index] = essential_flux_list[essential_matrix_index] * ub[essential_matrix_index]
		if not self.test_feasibility(lb, ub, S, b):
			return False, "not feasible enough"
		for essential_matrix_index in range(len(essential_flux_list)):
			if essential_matrix_index % 1000 == 0:
				print(essential_matrix_index / len(essential_flux_list))
			if essential_flux_list[essential_matrix_index] == 1:
				if not (lb[essential_matrix_index] > 0 or ub[essential_matrix_index] < 0):
					save_ub = copy.deepcopy(ub[essential_matrix_index])
					save_lb = copy.deepcopy(lb[essential_matrix_index])
					self.update_reaction_bounds(ordered_rxn[essential_matrix_index], 0, 0)
					lb[essential_matrix_index] = 0
					ub[essential_matrix_index] = 0
					if self.test_feasibility(lb, ub, S, b):
						return False, "Too feasible"
					else:
						self.update_reaction_bounds(ordered_rxn[essential_matrix_index], save_lb, save_ub)
						lb[essential_matrix_index] = save_ub
						ub[essential_matrix_index] = save_lb
		for i in range(len(ub)):
			self.update_reaction_bounds(ordered_rxn[i], save_lb_full[i], save_ub_full[i])
		return True, "Just right"
	
	def get_used_genes(self, nomenclature="HGNC"):
		used_genes = []
		for i in self.grRules:
			if nomenclature == "HGNC":
				if i != "[]":
					genes_found = re.findall("(HGNC:\d+)", i)
					for j in genes_found:
						if j not in used_genes:
							used_genes.append(j)
		return used_genes
	
	def generate_essential_flux_dataframe(self, n, output_dir, file_name_base, seed_list=[], print_progress=False, seed_name = False, **kwargs):
		lb, ub, S, b, ordered_rxn, ordered_met = self.dicts_to_mats()
		seed_empty = False
		if seed_list == []:
			seed_list = [random.random() for i in range(n)]
			seed_empty = True
		if seed_name:
			file_name = file_name_base+str(seed_list[0]).replace(".","")
		value_list = np.zeros((n, len(self.rxn_dict)))
		count = 0
		seed_fail = []
		failed = False
		for seed_val in seed_list:
			print(seed_val)
			if seed_empty:
				essential_flux = self.find_essential_fluxes_fast(print_progress=print_progress)
			else:
				essential_flux = self.find_essential_fluxes_fast(print_progress=print_progress, seed=seed_val)
			is_essential, fail_desc = self.test_essential_list(essential_flux)
			if is_essential:
				for i in essential_flux:
					value_list[count, i] += 1
			else:
				failed = True
				seed_fail.append([seed_val, fail_desc])
			count += 1
		if failed:
			print(f"Failed, description of failures:\n{seed_fail}")
		else:
			np.save(output_dir/file_name,pd.DataFrame(value_list).to_numpy())
			np.save(output_dir / (file_name_base+"header"), ordered_rxn)
	
	def get_met_comp_dict(self,met_name):
		chem_count_dict = {"C": 0, "H": 0, "N": 0, "O": 0, "P": 0, "S": 0, "Ca":0,"Cl":0, "K":0, "Na":0,"Se":0,"Co":0,"I":0,"Fe":0,"R": 0,"Ra": 0,"Rb": 0,"Rc": 0,"X": 0,"Y": 0}
		comp = self.met_dict[met_name]['comp']
		flag = False
		if "R" in comp:
			flag = True
		comp = comp.replace("FULLR3", "Rc")
		comp = comp.replace("FULLR2", "Rb")
		comp = comp.replace("FULLR", "Ra")
		search = re.findall("[A-Z][a-z]?\d*", comp)
		for i in search:
			chem = re.findall("[A-Z][a-z]?", i)[0]
			count = 1
			if len(re.findall("\d+", i)) == 1:
				count = re.findall("\d+", i)[0]
			try:
				chem_count_dict[chem] += int(count)
			except KeyError:
				print("Key not found")
				print(chem)
				return {}
		return chem_count_dict
	
	def get_met_mol_weight(self,met_name):
		mol_weight_dict = {"H":1.008,"C":12.011,"N":14.007,"O":15.999,"P":30.97376,"S":32.06,"Ca":40.078,"Cl":35.45,"K":39.0983,"Na":22.98977,"Se":78.971,"Co":58.93319,"I":126.9045,"Fe":55.845,"R":12.011*15+1.008*31,"Ra":(12.011*15+1.008*31)*1,"Rb":(12.011*15+1.008*31)*1,"Rc":(12.011*15+1.008*31)*1,"X": 0,"Y": 0}
		chem_count_dict = self.get_met_comp_dict(met_name)
		mw = 0
		for chem in chem_count_dict.keys():
			mw += int(chem_count_dict[chem])*mol_weight_dict[chem]
		return mw
	
	def get_rxn_mass_change(self,rxn_name):
		st = self.rxn_dict[rxn_name]["S"]
		mass_change = 0
		for met in st.keys():
			mass_change += st[met]*self.get_met_mol_weight(met)
		return mass_change
		
	def get_rxn_element_change(self,rxn_name):
		chem_count_dict = {"C": 0, "H": 0, "N": 0, "O": 0, "P": 0, "S": 0,"Ca":0,"Cl":0,"K":0, "Na":0,"Se":0,"Co":0,"I":0,"Fe":0, "R": 0,"Ra": 0,"Rb": 0,"Rc": 0,"X": 0,"Y": 0}
		st = self.rxn_dict[rxn_name]["S"]
		for met in st.keys():
			chem_count = self.get_met_comp_dict(met)
			for chem in chem_count.keys():
				chem_count_dict[chem] += st[met]*chem_count[chem]
		return chem_count_dict
		

		

	def build_essential_flux_matrix_from_dataframes(self, list_of_dataframes):
		
		for dataframe_ind in range(len(list_of_dataframes)):
			if dataframe_ind == 0:
				array = list_of_dataframes[dataframe_ind].to_numpy()
			else:
				array = np.vstack((array, list_of_dataframes[dataframe_ind].to_numpy()))
		return array
	
	def generate_essential_flux_model_from_essential_flux_matrix(self, essential_flux_matrix, use_gp_env=True):
		essential_flux_scores = np.average(essential_flux_matrix, axis=0)
		least_essential_flux_indices = np.argsort(essential_flux_scores)
		
		most_essential_flux_indices = np.flip(least_essential_flux_indices)
		
		saved_ub = copy.deepcopy(self.ub)
		saved_lb = copy.deepcopy(self.lb)
		
		for i in most_essential_flux_indices:
			if (self.ub[i] < 0 or self.lb[i] > 0):
				continue
			else:
				self.update_reaction_bounds(i, 0, 0)
		essential_flux_main = []
		for i in most_essential_flux_indices:
			if (self.ub[i] < 0 or self.lb[i] > 0):
				essential_flux_main.append(i)
				continue
			else:
				self.update_reaction_bounds(i, saved_lb[i], saved_ub[i])
				if self.test_feasibility():
					essential_flux_main.append(i)
					break
				else:
					essential_flux_main.append(i)
		
		most_essential_flux_names = [self.reaction_names[i] for i in most_essential_flux_indices]
		essential_flux_main_names = [self.reaction_names[i] for i in essential_flux_main]
		not_essential_fluxes = list(set(most_essential_flux_names) - set(essential_flux_main_names))
		
		self.del_reaction(not_essential_fluxes)
		self.purge_metabolites()
		flux_limits = self.fva()
		for i in range(len(flux_limits)):
			self.lb[i] = flux_limits[i][0]
			self.ub[i] = flux_limits[i][1]

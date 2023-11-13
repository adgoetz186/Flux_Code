import pickle
import sys
from pympler import asizeof
import gurobipy as gp
import copy
import re
import pandas as pd
from gurobipy import GRB
import scipy as sc
import scipy.linalg
import time
import scipy.sparse as sp
import numpy as np
import scipy.io as spio
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
	
	def __init__(self, model_name=None, S=None, lb=None, ub=None, b=None, **kwargs):
		self.model_name = model_name
		self.S = S
		self.lb = lb
		self.ub = ub
		self.b = b
		
		self.env = None
		
		self.pinched_reactions = None
		# free this up by tying it to a dict
		self.grRules = None
		self.genes = None
		self.metabolite_comp = None
		if self.S != None:
			self.metabolite_names = [f"metabolite_{i}" for i in range(np.shape(S)[0])]
			self.reaction_names = [f"Reaction_{i}" for i in range(np.shape(S)[1])]
		else:
			self.metabolite_names = None
			self.reaction_names = None
		
		if "pinched_reactions" in kwargs.keys():
			self.pinched_reactions = kwargs["pinched_reactions"]
		
		if "reaction_names" in kwargs.keys():
			self.reaction_names = kwargs["reaction_names"]
		
		if "grRules" in kwargs.keys():
			self.grRules = kwargs["grRules"]
		
		if "genes" in kwargs.keys():
			print("genes are not currently directly used by model, genes are added and handled through grRules")
			self.genes = kwargs["genes"]
		
		if "metabolite_names" in kwargs.keys():
			self.metabolite_names = kwargs["metabolite_names"]
		
		if "metabolite_comp" in kwargs.keys():
			self.metabolite_comp = kwargs["metabolite_comp"]
	
	def add_gp_key_env_to_model(self, TokenServer_name):
		# for hpg, TokenServer_name should be 'grb-ts.ufhpc'
		env = gp.Env(empty=True)
		env.setParam('TokenServer', TokenServer_name)
		env.start()
		self.env = env
	
	def save_model_as_pkl(self, filepath):
		dict_of_model = {}
		dict_of_model["model_name"] = self.model_name
		dict_of_model["S"] = self.S
		dict_of_model["lb"] = self.lb
		dict_of_model["ub"] = self.ub
		dict_of_model["b"] = self.b
		dict_of_model["pinched_reactions"] = self.pinched_reactions
		dict_of_model["reaction_names"] = self.reaction_names
		dict_of_model["grRules"] = self.grRules
		dict_of_model["genes"] = self.genes
		dict_of_model["metabolite_names"] = self.metabolite_names
		dict_of_model["metabolite_comp"] = self.metabolite_comp
		with open(filepath.parent/(filepath.stem+".pkl"), "wb") as outfile:
			pickle.dump(dict_of_model, outfile)
	
	def load_pkl_model(self, filepath):
		# Loads a pickled model which was saved using the "save_model_as_pkl" function from this class
		with open(filepath, "rb") as readfile:
			model_dict = pickle.load(readfile)
		self.model_name = model_dict["model_name"]
		self.S = model_dict["S"]
		self.lb = model_dict["lb"]
		self.ub = model_dict["ub"]
		self.b = model_dict["b"]
		self.pinched_reactions = model_dict["pinched_reactions"]
		self.reaction_names = model_dict["reaction_names"]
		self.grRules = model_dict["grRules"]
		self.genes = model_dict["genes"]
		self.metabolite_names = model_dict["metabolite_names"]
		self.metabolite_comp = model_dict["metabolite_comp"]
	
	def load_mat_model(self, filename, read_model_name="", stored_model_name="", model_comp=None):
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
		mat = loadmat(filename)
		if read_model_name == "":
			size_list = [asizeof.asizeof(mat[key]) for key in mat.keys()]
			read_model_name = list(mat.keys())[size_list.index(max([asizeof.asizeof(mat[key]) for key in mat.keys()]))]
		if stored_model_name == "":
			stored_model_name = read_model_name
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
			raise ValueError("model_comp is missing keys")
		
		try:
			self.model_name = stored_model_name
			tried_key = model_comp["S"]
			self.S = sp.csr_matrix(mat[read_model_name][tried_key]).astype(np.float)
			tried_key = model_comp["lb"]
			self.lb = mat[read_model_name][tried_key].astype(np.float)
			tried_key = model_comp["ub"]
			self.ub = mat[read_model_name][tried_key].astype(np.float)
			tried_key = model_comp["b"]
			self.b = mat[read_model_name][tried_key].astype(np.float)
			tried_key = model_comp["pinched_reactions"]
			if tried_key == None:
				self.pinched_reactions = {}
			else:
				self.pinched_reactions = mat[read_model_name][tried_key].astype(np.float)
			
			tried_key = model_comp["reaction_names"]
			if tried_key == None:
				if self.S != None:
					self.reaction_names = [f"Reaction_{i}" for i in range(np.shape(self.S)[1])]
				else:
					self.reaction_names = None
			else:
				self.reaction_names = list(mat[read_model_name][tried_key])
			
			tried_key = model_comp["grRules"]
			if tried_key == None:
				self.grRules = None
			else:
				self.grRules = list(mat[read_model_name][tried_key])
			
			tried_key = model_comp["genes"]
			if tried_key == None:
				self.genes = None
			else:
				self.genes = list(mat[read_model_name][tried_key])
			
			tried_key = model_comp["metabolite_names"]
			if tried_key == None:
				if self.S != None:
					self.metabolite_names = [f"metabolite_{i}" for i in range(np.shape(self.S)[0])]
				else:
					self.metabolite_names = None
			else:
				self.metabolite_names = list(mat[read_model_name][tried_key])
			
			tried_key = model_comp["metabolite_comp"]
			if tried_key == None:
				self.metabolite_comp = None
			else:
				self.metabolite_comp = list(mat[read_model_name][tried_key])
		
		except KeyError:
			print(f"The model uses a different key for {tried_key}")
			print(f"possible keys: {mat[read_model_name].keys()}")
	
	def update_reaction_bounds(self, name, new_lb, new_ub):
		# Reaction name can be the name of the reaction or the index location of the reaction
		# If new_lb or new_ub is "keep" keeps old value
		if isinstance(name, str):
			name = self.reaction_names.index(name)
		if new_lb == "keep":
			new_lb = self.lb[name]
		if new_ub == "keep":
			new_ub = self.ub[name]
		if new_ub < new_lb:
			print(self.reaction_info(name))
			print(new_ub, new_lb)
			print(f"WARNING! Upper bound smaller than lower bound at {name}")
			input()
			new_ub = new_lb
		self.ub[name] = new_ub
		self.lb[name] = new_lb
	
	def limit_flux_bounds(self, absolute_max_value):
		for rxn in self.reaction_names:
			self.update_reaction_bounds(rxn, -1 * absolute_max_value, absolute_max_value)
	
	def add_reaction(self, S_col, lb, ub, **kwargs):
		# Reaction name can be the name of the reaction or the index location of the reaction
		print(np.shape(self.S))
		
		self.S = sp.csr_matrix(sp.hstack((self.S, S_col)))
		self.ub = np.hstack((self.ub, ub))
		self.lb = np.hstack((self.lb, lb))
		if "reaction_name" in kwargs.keys():
			self.reaction_names.append(kwargs["reaction_name"])
		else:
			self.reaction_names = [f"Reaction_{i}" for i in range(1, np.shape(self.S)[1] + 1)]
		
		if "Gene_Reaction_Rule" in kwargs.keys():
			if isinstance(self.grRules, list):
				self.grRules.append(kwargs["Gene_Reaction_Rule"])
		if "Gene_Reaction_Rule" not in kwargs.keys():
			if isinstance(self.grRules, list):
				print(
					"Gene Reaction Rule should be included, currently assuming no gene association with this reaction")
				self.grRules.append("[]")
	
	def purge_metabolites(self):
		# Finds and removes all metabolites which do not have any associated reactions
		# This can lead to issues if used after reactions are pinched
		metabolite_usage = np.sum(np.abs(self.S), axis=1)
		metabolites_to_remove = []
		for i in range(np.size(metabolite_usage)):
			if metabolite_usage[i, 0] == 0:
				metabolites_to_remove.append(i)
		metabolites_to_remove.sort(reverse=True)
		s_placeholder = self.S.todense()
		for i in metabolites_to_remove:
			s_placeholder = np.delete(s_placeholder, i, 0)
			if self.metabolite_comp != "None":
				del self.metabolite_comp[i]
			if self.metabolite_names != "None":
				del self.metabolite_names[i]
			self.b = np.delete(self.b, i, 0)
		self.S = sp.csr_matrix(s_placeholder)
	
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
	
	def del_reaction(self, name):
		# name can be the name of the reaction or its index
		# deletes reaction from either normal reactions or pinched reactions
		# only deletes pinched reactions if specified by name of reaction
		# capability to handle list added to improve performance as going to and from dense matrix is expensive
		if not isinstance(name, list):
			if isinstance(name, str):
				try:
					name = self.reaction_names.index(name)
				except ValueError:
					try:
						del self.pinched_reactions[name]
					except ValueError:
						print(f"The reaction you tried to delete called {name}, never existed in the first place.")
			self.S = sp.csr_matrix(np.delete(self.S.todense(), name, 1))
			self.lb = np.delete(self.lb, name, 0)
			self.ub = np.delete(self.ub, name, 0)
			if isinstance(self.grRules, list):
				del self.grRules[name]
			del self.reaction_names[name]
		else:
			s_placeholder = self.S.todense()
			for name_entry_index in range(len(name)):
				if isinstance(name[name_entry_index], str):
					try:
						name[name_entry_index] = self.reaction_names.index(name[name_entry_index])
					except ValueError:
						try:
							del self.pinched_reactions[name[name_entry_index]]
						except ValueError:
							print(
								f"The reaction you tried to delete called {name[name_entry_index]}, never existed in the first place.")
			name.sort(reverse=True)
			for name_entry in name:
				s_placeholder = np.delete(s_placeholder, name_entry, 1)
				self.lb = np.delete(self.lb, name_entry, 0)
				self.ub = np.delete(self.ub, name_entry, 0)
				if isinstance(self.grRules, list):
					del self.grRules[name_entry]
				del self.reaction_names[name_entry]
			self.S = sp.csr_matrix(s_placeholder)
	
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
	
	def fva(self, **kwargs):
		# Returns list of [min,max] elements
		# If reactions are edited this fva no longer applies
		# This is valid for removed reactions, but could lead to issues if 2 reactions are simply switched around
		if self.env != None:
			model = gp.Model("fva_model_" + self.model_name, env=self.env)
		else:
			model = gp.Model("fva_model_" + self.model_name)
		react_flux = model.addMVar(shape=self.S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=self.lb,
		                           ub=self.ub)
		model.addConstr(self.S @ react_flux == self.b, name="c")
		if "extra_constraints" in kwargs.keys():
			cons = 0
			for constraint in kwargs["extra_constraints"]:
				eval(f"model.addConstr({constraint},name='{cons}')")
				cons += 1
		model.Params.LogToConsole = 0
		obj = np.zeros((1, self.S.shape[1]))
		variability_list = [[0, 0] for i in range(len(self.lb))]
		if "index_list" in kwargs.keys():
			i_list = kwargs["index_list"]
		else:
			i_list = range(len(self.lb))
		for i in i_list:
			print(len(self.lb) - i)
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
			min_max.append(-1 * model.objVal)
			min_max.sort()
			variability_list[i] = min_max
		return variability_list
	
	def test_feasibility(self, **kwargs):
		if self.env != None:
			model = gp.Model("fva_model_" + self.model_name, env=self.env)
		else:
			model = gp.Model("fva_model_" + self.model_name)
		if "unique_bounds" in kwargs.keys():
			tflb = kwargs["unique_bounds"][0]
			tfub = kwargs["unique_bounds"][1]
		else:
			tflb = self.lb
			tfub = self.ub
		
		react_flux = model.addMVar(shape=self.S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=tflb,
		                           ub=tfub)
		model.addConstr(self.S @ react_flux == self.b, name="c")
		if "extra_constraints" in kwargs.keys():
			cons = 0
			for constraint in kwargs["extra_constraints"]:
				eval(f"model.addConstr({constraint},name='{cons}')")
				cons += 1
		# model.Params.LogToConsole = 0
		obj = np.zeros((1, self.S.shape[1]))
		model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.Params.LogToConsole = 0
		model.Params.NumericFocus = 3
		# model.Params.BarHomogeneous = 1
		model.Params.ScaleFlag = 0
		# model.Params.NumericFocus = 3
		# model.Params.Presolve = 0
		# model.optimize()
		model.Params.FeasibilityTol = 1e-9
		
		# model.Params.Aggregate = 0
		model.optimize()
		# if model.Status == 2:
		# model.printQuality()
		# print(model.Status)
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
	
	def get_exchanged_reaction_info(self):
		indx = np.reshape(np.argwhere((np.ravel(np.sum(np.abs(self.S), axis=0)) == 1)), (-1))
		indc = np.zeros_like(indx)
		for i in range(np.size(indx)):
			indc[i] = self.S[:, indx[i]].nonzero()[0][0]
		exchange_reaction_dictionary = {}
		for i in range(len(indx)):
			exchange_reaction_dictionary[self.metabolite_names[indc[i]]] = self.reaction_names[indx[i]]
		return exchange_reaction_dictionary
	
	def metabolite_info(self, name, **kwargs):
		if isinstance(name, str):
			name = self.metabolite_names.index(name)
		stndrdth = "th"
		if str(name)[-1] == "1":
			stndrdth = "st"
		elif str(name)[-1] == "2":
			stndrdth = "nd"
		elif str(name)[-1] == "3":
			stndrdth = "rd"
		print(f"{self.metabolite_names[name]} is the {name}{stndrdth} metabolite")
		print(f"The reactions in which {self.metabolite_names[name]} is produced or consumed are:")
		denseS = self.S.todense()[name]
		print("start")
		for i in list(np.argwhere(denseS != 0)[:, 1]):
			print(self.reaction_names[i])
		print("stop")
	
	# if np.size(self.feasible_space_points) != 0:
	#    point_to_use = 0
	#    if "point_number" in kwargs.keys():
	#        point_to_use = kwargs["point_number"]
	#    flux = np.ravel(np.multiply(denseS[name],self.feasible_space_points[point_to_use]))
	# total_flux = 0
	# for i in col_values:
	#    total_flux+= flux[i]
	#    if "hide_small_flux" in kwargs.keys():
	#        if kwargs["hide_small_flux"] < abs(flux[i]):
	#            print(f"Reaction is {self.reaction_names[i]} (index of {i}), metabolite production rate is {flux[i]} (stochiometric constant {self.S[name,i]})")
	#    else:
	#        print(f"Reaction is {self.reaction_names[i]} (index of {i}), metabolite production rate is {flux[i]} (stochiometric constant {self.S[name, i]})")
	# print(f"The total flux is {total_flux} (should be 0)")
	
	def reaction_info(self, name):
		# Reaction name can be the name of the reaction or the index location of the reaction
		if isinstance(name, str):
			name = self.reaction_names.index(name)
		if str(name) in ["11", "12", "13"]:
			stndrdth = "th"
		elif str(name)[-1] == "1":
			stndrdth = "st"
		elif str(name)[-1] == "2":
			stndrdth = "nd"
		elif str(name)[-1] == "3":
			stndrdth = "rd"
		else:
			stndrdth = "th"
		print(f"{self.reaction_names[name]} is the {name}{stndrdth} reaction")
		print(f"Lower Bound: {self.lb[name]}, Upper Bound {self.ub[name]}")
		print(f"The relevant chemical species participating in the reaction are:")
		row_values = self.S[:, name].nonzero()[0]
		for i in row_values:
			print(
				f"Chemical Species is {self.metabolite_names[i]} (index of {i}), stochiometric constant is {self.S[i, name]}")
		if isinstance(self.grRules, list):
			print(f"The reaction gene rule is:\n{self.grRules[name]}")
	
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
	
	def lin_dep_np_ctp(self, old_basis, origin, new_point):
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
		true_old_rank = np.linalg.matrix_rank(old_basis - origin)
		
		# make new basis set
		new_basis = np.vstack((old_basis, new_point))
		
		# find the dimensionality of the new basis set
		new_rank = np.linalg.matrix_rank(new_basis - origin)
		
		# if the new basis has a higher dimensionality, the added vector is indeed spanning and is kept
		if new_rank <= true_old_rank:
			return old_basis, true_old_rank
		return new_basis, new_rank
	
	def generate_warmup_gb(self, max_search_number,required_rank = np.inf, **kwargs):
		# Description:
		#   Generates warmup vectors for coordinate direction hit and run
		# Input:
		#   max_search_number: integer value which determines how many proposed vectors are subsequently rejected before
		#       the algorithm terminates
		# Output:
		#   2d array with each row giving a point in the feasible space the top row provides an origin point
		#       for the space while the following points gives a spanning set of vectors wrt the origin point
		
		# Defines Model
		env = gp.Env(empty=True)
		#env.setParam('TokenServer', 'grb-ts.ufhpc')
		env.start()
		model = gp.Model("warmup_model_" + self.model_name, env=env)
		react_flux = model.addMVar(shape=self.S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=self.lb,
		                           ub=self.ub)
		model.addConstr(self.S @ react_flux == self.b, name="c")
		
		# Suppress printing to console
		model.Params.LogToConsole = 0
		
		# Defines origin and points as empty
		origin = np.array([])
		points = np.array([])
		
		# Initializes number of subsequent vectors which have been rejected
		search_number = 0
		
		# Initializes rank
		rank = -1
		while search_number < max_search_number and rank < required_rank:
			# Finds random vertex
			while model.Status == 1:
				obj = (2 * np.random.random_sample(self.S.shape[1]) - 1)
				model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
				model.optimize()
			
			# The below ensures the points fall inside the bounds
			# Does seem to potentially break constraints
			# entries of points are just being moved, the points should be scaled
			# if this becomes a problem, define scaling function
			point_to_use = np.minimum(react_flux.X, self.ub)
			point_to_use = np.maximum(point_to_use, self.lb)
			# Saves rank to later check for rank change
			old_rank = rank
			if search_number%100 == 1:
				print(search_number,rank)
			# The first point becomes and origin point
			# The rank is 0 with the origin point since no vectors have been defined
			if np.size(origin) == 0:
				origin = point_to_use
				rank += 1
			# Next a slightly special bit of code creates the first basis vector
			elif np.size(points) == 0:
				test_point, test_rank = self.lin_dep_np_ctp(origin, origin, point_to_use)
				if test_rank == 1:
					points = np.reshape(test_point[1], (1, -1))
				rank = test_rank
			# Finally points are added in through a straightforward approach
			else:
				points, rank = self.lin_dep_np_ctp(points, origin, point_to_use)
			
			# Reset the model to discard the current solution
			model.reset()
			
			# Reset the search_number if a vector is accepted as basis
			if rank != old_rank:
				search_number = 0
			else:
				search_number += 1
		print(rank)
		# returns set of points with top point being basis
		return np.vstack((origin, points))
	
	def generate_vertex_samples(self, n,max_combined_flux = -1, **kwargs):
		# Description:
		#   Generates warmup vectors for coordinate direction hit and run
		# Input:
		#   max_search_number: integer value which determines how many proposed vectors are subsequently rejected before
		#       the algorithm terminates
		# Output:
		#   2d array with each row giving a point in the feasible space the top row provides an origin point
		#       for the space while the following points gives a spanning set of vectors wrt the origin point
		
		# Defines Model
		env = gp.Env(empty=True)
		#env.setParam('TokenServer', 'grb-ts.ufhpc')
		env.start()
		model = gp.Model("warmup_model_" + self.model_name, env=env)
		react_flux = model.addMVar(shape=self.S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=self.lb,
		                           ub=self.ub)
		model.addConstr(self.S @ react_flux == self.b, name="c")
		if max_combined_flux != -1:
			all_flux = self.objective_function_generator(self.reaction_names,
			                                              [1 for i in range(np.size(self.reaction_names))])
			model.addConstr(react_flux @ all_flux <= max_combined_flux, name="max_flux")
		
		# Suppress printing to console
		model.Params.LogToConsole = 0
		
		# Initializes number of subsequent vectors which have been rejected
		search_number = 0
		
		# Initializes rank
		point_matrix = np.zeros((n,len(self.lb)))
		for i in range(n):
			if i % 100 == 0:
				print(i)
			# Finds random vertex
			while model.Status == 1:
				obj = (2 * np.random.random_sample(self.S.shape[1]) - 1)
				model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
				model.optimize()
			point_matrix[i] = react_flux.X
			model.reset()
			
		return point_matrix
	
	def ACHRSampler_Display(self, c, warmupPoints, fileName, number_of_points, stepsPerPoint, store_points,
	                        save_fig=False, dir="", file_name="", **kwargs):
		maxMinTol = 1e-9
		uTol = 1e-9
		dTol = 1e-9
		prev = 0
		
		# find center
		centerPoint = np.mean(warmupPoints, axis=0)
		print("centerPoint")
		# flux_balance_test(centerPoint,S,b,lb,ub,c)
		
		if "maxTime" in kwargs.keys():
			maxTime = kwargs["maxTime"]
		else:
			maxTime = 1000 * 3600
		
		if "initial_point" in kwargs.keys():
			prev_point = kwargs["initial_point"]
		else:
			prev_point = centerPoint
		
		Sa = self.S.toarray()
		N = sc.linalg.null_space(Sa)
		
		nWrmup = np.shape(warmupPoints)[0]
		
		total_steps = 0
		total_steps_needed = number_of_points * stepsPerPoint
		
		# random_vector = np.random.choice(total_steps_needed)
		print(warmupPoints)
		t0 = time.time()
		
		val = np.argpartition(np.max(warmupPoints, axis=0) - np.min(warmupPoints, axis=0), -2)[-2:]
		for j in warmupPoints:
			plt.scatter(j[val[0]], j[val[1]], color="gray", alpha=.3)
		plt.hlines([self.lb[val[1]], self.ub[val[1]]], self.lb[val[0]] - 10, self.ub[val[0]] + 10)
		plt.vlines([self.lb[val[0]], self.ub[val[0]]], self.lb[val[1]] - 10, self.ub[val[1]] + 10)
		plt.scatter(prev_point[val[0]], prev_point[val[1]])
		oc = np.mean(warmupPoints, axis=0)
		plt.scatter(oc[val[0]], oc[val[1]], color="gray")
		plt.scatter(centerPoint[val[0]], centerPoint[val[1]], color="black")
		x_range = self.ub[val[0]] - self.lb[val[0]]
		y_range = self.ub[val[1]] - self.lb[val[1]]
		plt.ylim([self.lb[val[1]] - .05 * y_range, self.ub[val[1]] + .05 * y_range])
		plt.xlim([self.lb[val[0]] - .05 * x_range, self.ub[val[0]] + .05 * x_range])
		if save_fig:
			sfh = 1
			plt.savefig(f'{dir}{file_name}/{sfh}.png')
			plt.clf()
			sfh += 1
		else:
			plt.show()
		
		points = np.zeros((number_of_points, len(warmupPoints[0])))
		
		point_count = 1
		
		while point_count <= number_of_points:
			randVector = np.random.random_sample(stepsPerPoint)
			# handle this better
			
			step_count = 1
			
			while step_count <= stepsPerPoint:
				# print(step_count)
				# Pick a random warmup point
				random_row_index = np.random.choice(np.shape(warmupPoints)[0])
				rand_point = warmupPoints[random_row_index]
				
				# Get a direction from the center point to the warmup point
				u = self.normalize(rand_point - centerPoint)
				
				if "display" in kwargs.keys():
					val = np.argpartition(np.max(warmupPoints, axis=0) - np.min(warmupPoints, axis=0), -2)[-2:]
					for j in warmupPoints:
						plt.scatter(j[val[0]], j[val[1]], color="gray", alpha=.3)
					oc = np.mean(warmupPoints, axis=0)
					plt.hlines([self.lb[val[1]], self.ub[val[1]]], self.lb[val[0]] - 10, self.ub[val[0]] + 10)
					plt.vlines([self.lb[val[0]], self.ub[val[0]]], self.lb[val[1]] - 10, self.ub[val[1]] + 10)
					plt.scatter(oc[val[0]], oc[val[1]], color="gray")
					plt.scatter(prev_point[val[0]], prev_point[val[1]])
					plt.scatter(rand_point[val[0]], rand_point[val[1]])
					plt.scatter(centerPoint[val[0]], centerPoint[val[1]], color="black")
					for i in range(point_count):
						plt.scatter(points[i][val[0]], points[i][val[1]], color="gray", alpha=0.3)
					x_range = self.ub[val[0]] - self.lb[val[0]]
					y_range = self.ub[val[1]] - self.lb[val[1]]
					plt.ylim([self.lb[val[1]] - .05 * y_range, self.ub[val[1]] + .05 * y_range])
					plt.xlim([self.lb[val[0]] - .05 * x_range, self.ub[val[0]] + .05 * x_range])
					if kwargs["display"][0] == "save":
						plt.savefig(f'{kwargs["display"][1]}/{kwargs["display"][2]}{sfh}.png')
						plt.clf()
						sfh += 1
					else:
						plt.show()
					val = np.argpartition(np.max(warmupPoints, axis=0) - np.min(warmupPoints, axis=0), -2)[-2:]
					for j in warmupPoints:
						plt.scatter(j[val[0]], j[val[1]], color="gray", alpha=.3)
					print(val)
					oc = np.mean(warmupPoints, axis=0)
					plt.hlines([self.lb[val[1]], self.ub[val[1]]], self.lb[val[0]] - 10, self.ub[val[0]] + 10)
					plt.vlines([self.lb[val[0]], self.ub[val[0]]], self.lb[val[1]] - 10, self.ub[val[1]] + 10)
					plt.scatter(oc[val[0]], oc[val[1]], color="gray")
					plt.scatter(prev_point[val[0]], prev_point[val[1]])
					plt.scatter(rand_point[val[0]], rand_point[val[1]])
					plt.scatter(centerPoint[val[0]], centerPoint[val[1]], color="black")
					for i in range(point_count):
						plt.scatter(points[i][val[0]], points[i][val[1]], color="gray", alpha=0.3)
					plt.plot([prev_point[val[0]] - 1000 * u[val[0]], prev_point[val[0]],
					          prev_point[val[0]] + 1000 * u[val[0]]],
					         [prev_point[val[1]] - 1000 * u[val[1]], prev_point[val[1]],
					          prev_point[val[1]] + 1000 * u[val[1]]])
					x_range = self.ub[val[0]] - self.lb[val[0]]
					y_range = self.ub[val[1]] - self.lb[val[1]]
					plt.ylim([self.lb[val[1]] - .05 * y_range, self.ub[val[1]] + .05 * y_range])
					plt.xlim([self.lb[val[0]] - .05 * x_range, self.ub[val[0]] + .05 * x_range])
					if kwargs["display"][0] == "save":
						plt.savefig(f'{kwargs["display"][1]}/{kwargs["display"][2]}{sfh}.png')
						plt.clf()
						sfh += 1
					else:
						plt.show()
				
				# Figure out the distances to upper and lower bounds
				distUb = (self.ub - prev_point)
				distLb = (prev_point - self.lb)
				
				# Figure out if we are too close to a boundary
				validDir = ((distUb > dTol) & (distLb > dTol))
				
				posValidDir = ((distUb > dTol) & (distLb > dTol)) & (u > uTol)
				
				negValidDir = ((distUb > dTol) & (distLb > dTol)) & (u < (-1 * uTol))
				
				posmaxStepTemp = np.compress(posValidDir, distUb, axis=0) / np.compress(posValidDir, u, axis=0)
				negmaxStepTemp = np.compress(negValidDir, distUb, axis=0) / np.compress(negValidDir, u, axis=0)
				
				negminStepTemp = -1 * np.compress(negValidDir, distLb, axis=0) / np.compress(negValidDir, u, axis=0)
				posminStepTemp = -1 * np.compress(posValidDir, distLb, axis=0) / np.compress(posValidDir, u, axis=0)
				
				# print(posmaxStepTemp)
				# print(negmaxStepTemp)
				
				# print(negminStepTemp)
				# print(posminStepTemp)
				# input()
				
				maxStep = np.min(np.hstack((posmaxStepTemp, negminStepTemp)))
				minStep = np.max(np.hstack((posminStepTemp, negmaxStepTemp)))
				
				# print(f"maxStep = {maxStep}")
				# print(f"minStep = {minStep}")
				# input()
				
				if ((abs(minStep) < maxMinTol) & (abs(maxStep) < maxMinTol)) or (minStep > maxStep):
					print(f"Warning: \nMin step = {minStep}\n Max step = {maxStep}")
					continue
				
				# grab random value from pregenerated list for the step distance
				
				stepDist = randVector[step_count - 1] * (maxStep - minStep) + minStep
				curPoint = prev_point + stepDist * u
				
				if "objective_constraint" in kwargs.keys():
					ocl = kwargs["objective_constraint"]
					if ocl[0] == ">":
						cnt = 0
						while np.matmul(curPoint, np.transpose(c)) < ocl[1] * .999:
							if "display" in kwargs.keys():
								print(np.max(warmupPoints, axis=0))
								print(np.min(warmupPoints, axis=0))
								
								print(np.max(warmupPoints, axis=0) - np.min(warmupPoints, axis=0))
								val = np.argpartition(np.max(warmupPoints, axis=0) - np.min(warmupPoints, axis=0),
								                      -2)[
								      -2:]
								for j in warmupPoints:
									plt.scatter(j[val[0]], j[val[1]], color="gray", alpha=.3)
								print(val)
								plt.hlines([self.lb[val[1]], self.ub[val[1]]], self.lb[val[0]] - 10,
								           self.ub[val[0]] + 10)
								plt.vlines([self.lb[val[0]], self.ub[val[0]]], self.lb[val[1]] - 10,
								           self.ub[val[1]] + 10)
								oc = np.mean(warmupPoints, axis=0)
								plt.scatter(oc[val[0]], oc[val[1]], color="gray")
								plt.scatter(prev_point[val[0]], prev_point[val[1]])
								plt.scatter(rand_point[val[0]], rand_point[val[1]])
								plt.scatter(curPoint[val[0]], curPoint[val[1]], color="red")
								plt.scatter(centerPoint[val[0]], centerPoint[val[1]], color="black")
								plt.plot([prev_point[val[0]] - 1000 * u[val[0]], prev_point[val[0]],
								          prev_point[val[0]] + 1000 * u[val[0]]],
								         [prev_point[val[1]] - 1000 * u[val[1]], prev_point[val[1]],
								          prev_point[val[1]] + 1000 * u[val[1]]])
								x_range = self.ub[val[0]] - self.lb[val[0]]
								y_range = self.ub[val[1]] - self.lb[val[1]]
								plt.ylim([self.lb[val[1]] - .05 * y_range, self.ub[val[1]] + .05 * y_range])
								plt.xlim([self.lb[val[0]] - .05 * x_range, self.ub[val[0]] + .05 * x_range])
								for i in range(point_count):
									plt.scatter(points[i][val[0]], points[i][val[1]], color="gray", alpha=0.3)
								if kwargs["display"][0] == "save":
									plt.savefig(f'{kwargs["display"][1]}/{kwargs["display"][2]}{sfh}.png')
									plt.clf()
									sfh += 1
								else:
									plt.show()
							if cnt % 100000 == 0 and cnt > 0:
								print(
									f"Several point failures, current value is: {np.matmul(curPoint, np.transpose(c))} objective is {ocl[1]}, {np.matmul(curPoint, np.transpose(c)) / ocl[1]}")
							stepDist = np.random.random(1)[0] * (maxStep - minStep) + minStep
							curPoint = prev_point + stepDist * u
							cnt += 1
					elif ocl[0] == "<":
						while curPoint * np.transpose(c) > ocl[1] * 1.001:
							print("count2")
							stepDist = np.random.random(1)[0] * (maxStep - minStep) + minStep
							curPoint = prev_point + stepDist * u
				
				# Advance to the next point
				
				if "display" in kwargs.keys():
					val = np.argpartition(np.max(warmupPoints, axis=0) - np.min(warmupPoints, axis=0), -2)[-2:]
					for j in warmupPoints:
						plt.scatter(j[val[0]], j[val[1]], color="gray", alpha=.3)
					plt.hlines([self.lb[val[1]], self.ub[val[1]]], self.lb[val[0]] - 10, self.ub[val[0]] + 10)
					plt.vlines([self.lb[val[0]], self.ub[val[0]]], self.lb[val[1]] - 10, self.ub[val[1]] + 10)
					oc = np.mean(warmupPoints, axis=0)
					plt.scatter(oc[val[0]], oc[val[1]], color="gray")
					plt.scatter(prev_point[val[0]], prev_point[val[1]])
					plt.scatter(rand_point[val[0]], rand_point[val[1]])
					plt.scatter(curPoint[val[0]], curPoint[val[1]])
					plt.scatter(centerPoint[val[0]], centerPoint[val[1]], color="black")
					plt.plot([prev_point[val[0]] - 1000 * u[val[0]], prev_point[val[0]],
					          prev_point[val[0]] + 1000 * u[val[0]]],
					         [prev_point[val[1]] - 1000 * u[val[1]], prev_point[val[1]],
					          prev_point[val[1]] + 1000 * u[val[1]]])
					x_range = self.ub[val[0]] - self.lb[val[0]]
					y_range = self.ub[val[1]] - self.lb[val[1]]
					plt.ylim([self.lb[val[1]] - .05 * y_range, self.ub[val[1]] + .05 * y_range])
					plt.xlim([self.lb[val[0]] - .05 * x_range, self.ub[val[0]] + .05 * x_range])
					for i in range(point_count):
						plt.scatter(points[i][val[0]], points[i][val[1]], color="gray", alpha=0.3)
					if kwargs["display"][0] == "save":
						plt.savefig(f'{kwargs["display"][1]}/{kwargs["display"][2]}{sfh}.png')
						plt.clf()
						sfh += 1
					else:
						plt.show()
					for j in warmupPoints:
						plt.scatter(j[val[0]], j[val[1]], color="gray", alpha=.3)
					print(val)
					plt.hlines([self.lb[val[1]], self.ub[val[1]]], self.lb[val[0]] - 10, self.ub[val[0]] + 10)
					plt.vlines([self.lb[val[0]], self.ub[val[0]]], self.lb[val[1]] - 10, self.ub[val[1]] + 10)
					oc = np.mean(warmupPoints, axis=0)
					plt.scatter(oc[val[0]], oc[val[1]], color="gray")
					plt.scatter(curPoint[val[0]], curPoint[val[1]])
					plt.scatter(centerPoint[val[0]], centerPoint[val[1]], color="black")
					x_range = self.ub[val[0]] - self.lb[val[0]]
					y_range = self.ub[val[1]] - self.lb[val[1]]
					plt.ylim([self.lb[val[1]] - .05 * y_range, self.ub[val[1]] + .05 * y_range])
					plt.xlim([self.lb[val[0]] - .05 * x_range, self.ub[val[0]] + .05 * x_range])
					for i in range(point_count):
						plt.scatter(points[i][val[0]], points[i][val[1]], color="gray", alpha=0.3)
					if kwargs["display"][0] == "save":
						plt.savefig(f'{kwargs["display"][1]}/{kwargs["display"][2]}{sfh}.png')
						plt.clf()
						
						sfh += 1
					else:
						plt.show()
				
				# I dont understand this well
				
				if total_steps % 100 == 0:
					# print("pt")
					# print(np.max(np.abs(np.matmul(Sa,curPoint))))
					if np.max(np.abs(np.matmul(Sa, curPoint))) > 1e-10:
						curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				
				# print out errors
				
				if total_steps % 200000 == 0 and False:
					print(f"Error {max(curPoint - self.ub)}, {max(self.lb - curPoint)}")
				# Both values should be negative
				
				time_elapsed = time.time() - t0
				
				# print step information
				
				if step_count % 500000 == 0:
					time_per_step = time_elapsed / step_count
					print(
						f"{point_count}, {step_count}, {time_elapsed / 60}, {(total_steps_needed - step_count) * time_per_step / 60}")
				
				indexer = np.arange(0, np.size(curPoint))
				overInd = indexer[curPoint > self.ub]
				
				underInd = indexer[curPoint < self.lb]
				# print(overInd)
				# print(underInd)
				# for i in overInd:
				#    print("over" , self.reaction_names[i], self.lb[i], curPoint[i], self.ub[i],step_count)
				# for i in underInd:
				#    print("under ", self.reaction_names[i], self.lb[i], curPoint[i], self.ub[i],step_count)
				if np.any((self.ub - curPoint) < 0) or np.any((curPoint - self.lb) < 0):
					# print("before")
					# self.flux_balance_test(curPoint, 0.1)
					# for i in range(len(u)):
					#    if abs(u[i]) < uTol:
					#        print(self.reactionnames[i])
					# This is where issue is happening
					# print("oi9")
					for j in overInd:
						curPoint[j] = self.ub[j]
					for j in underInd:
						curPoint[j] = self.lb[j]
				# print("stop")
				# print("After")
				# self.flux_balance_test(curPoint, 0.1)
				
				if total_steps % 2000000 == 0:
					print(f"fidErr {np.max(np.abs(np.matmul(Sa, curPoint)))}")
				# self.flux_balance_test(curPoint, 0.1)
				prev_point = curPoint
				step_count += 1
				
				total_steps += 1
				if np.floor(total_steps / total_steps_needed * 100) > prev:
					prev = total_steps / total_steps_needed * 100
					print(total_steps / total_steps_needed)
					print(np.sum(points))
					print(np.shape(points))
				
				if time_elapsed > maxTime:
					points[point_count] = curPoint
					np.savetxt(f"{fileName}.txt", points, delimiter=',')
					print(f"Time constraint reached after {point_count} points")
					return
			points[point_count - 1] = curPoint
			
			# recalculate the centerpart
			if "artificial_centering" in kwargs.keys():
				if kwargs["artificial_centering"]:
					centerPoint = (np.sum(warmupPoints, axis=0) + np.sum(points, axis=0)) / (
							np.shape(warmupPoints)[0] + point_count)
			
			point_count += 1
		if store_points:
			self.feasible_space_points = points
		np.savetxt(f"{fileName}.txt", points, delimiter=',')
		return points
	
	def HRSampler(self, origin_and_warmup_points, number_of_points, stepsPerPoint, **kwargs):
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
		
		maxMinTol = 1e-9
		uTol = 1e-9
		dTol = 1e-9
		u_sum_tol = 1e-9
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
		
		Sa = self.S.toarray()
		N = sc.linalg.null_space(Sa)
		
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
				distUb = (self.ub - prev_point)
				distLb = (prev_point - self.lb)
				
				# Figure out if we are too close to a boundary
				validDir = ((distUb > dTol) & (distLb > dTol))
				
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
				
				# if "objective_constraint" in kwargs.keys():
				#    ocl = kwargs["objective_constraint"]
				#    if ocl[0] == ">":
				#        cnt = 0
				#        while np.matmul(curPoint, np.transpose(c)) < ocl[1] * .999:
				#            if cnt % 100000 == 0 and cnt > 0:
				#                print(f"Several point failures, current value is: {np.matmul(curPoint, np.transpose(c))} objective is {ocl[1]}, {np.matmul(curPoint, np.transpose(c)) / ocl[1]}")
				#            stepDist = np.random.random(1)[0] * (maxStep - minStep) + minStep
				#            curPoint = prev_point + stepDist * u
				#            cnt += 1
				#    elif ocl[0] == "<":
				#        while curPoint * np.transpose(c) > ocl[1] * 1.001:
				#            print("count2")
				#            stepDist = np.random.random(1)[0] * (maxStep - minStep) + minStep
				#            curPoint = prev_point + stepDist * u
				
				if total_steps % 100 == 0:
					if np.max(np.abs(np.matmul(Sa, curPoint))) > 1e-10:
						curPoint = np.matmul(N, np.matmul(np.transpose(N), curPoint))
				
				if total_steps % 200000 == 0 and False:
					print(f"Error {max(curPoint - self.ub)}, {max(self.lb - curPoint)}")
				# Both values should be negative
				
				time_elapsed = time.time() - t0
				
				if step_count % 500000 == 0:
					time_per_step = time_elapsed / step_count
					print(
						f"{point_count}, {step_count}, {time_elapsed / 60}, {(total_steps_needed - step_count) * time_per_step / 60}")
				
				if np.any((self.ub - curPoint) < 0) or np.any((curPoint - self.lb) < 0):
					indexer = np.arange(0, np.size(curPoint))
					overInd = indexer[curPoint > self.ub]
					underInd = indexer[curPoint < self.lb]
					#print(overInd)
					#print(underInd)
					# scaling might be useful here
					for j in overInd:
						curPoint[j] = self.ub[j]
					for j in underInd:
						curPoint[j] = self.lb[j]
				
				if total_steps % 2000000 == 0:
					print(f"fidErr {np.max(np.abs(np.matmul(Sa, curPoint)))}")
				prev_point = curPoint
				step_count += 1
				
				total_steps += 1
				if np.floor(total_steps / total_steps_needed * 100) > prev:
					prev = total_steps / total_steps_needed * 100
					print(total_steps / total_steps_needed)
					print(np.sum(points))
					print(np.shape(points))
				
				if time_elapsed > maxTime:
					points[point_count] = curPoint
					print(f"Time constraint reached after {point_count} points")
					return points
			points[point_count - 1] = curPoint
			
			point_count += 1
		
		return points
	
	
	
	def convert_model_to_positive_flux(self):
		Sa = self.S.todense()
		print(Sa)
		ubc = copy.deepcopy(self.ub)
		# This does NOT handle c, handle c with the section defined to generate c
		neg_ub = np.reshape(np.argwhere(self.ub <= 0), (-1))
		for i in neg_ub:
			self.ub[i] = -1 * self.lb[i]
			self.lb[i] = -1 * ubc[i]
			Sa[:, i] *= -1
			self.reaction_names[i] = self.reaction_names[i] + "[inverted]"
		neg_lb = np.reshape(np.argwhere(self.lb < 0), (-1))
		new_ub = np.zeros(np.size(neg_lb))
		new_lb = np.zeros(np.size(neg_lb))
		new_S_col = np.zeros((np.shape(Sa)[0], np.size(neg_lb)))
		count = 0
		for i in neg_lb:
			new_ub[count] = -1 * self.lb[i]
			new_S_col[:, count] = -1 * np.reshape(Sa[:, i], (-1))
			self.lb[i] = 0
			self.reaction_names.append(self.reaction_names[i] + "[inverted]")
			count += 1
		self.lb = np.concatenate((self.lb, new_lb))
		self.ub = np.concatenate((self.ub, new_ub))
		self.S = sp.csr_matrix(np.hstack((Sa, new_S_col)))
	
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
		essential_reactions = []
		true_lb = copy.deepcopy(self.lb)
		true_ub = copy.deepcopy(self.ub)
		
		true_lb_2 = copy.deepcopy(self.lb)
		true_ub_2 = copy.deepcopy(self.ub)
		if self.env != None:
			model = gp.Model("fva_model_" + self.model_name, env=self.env)
		else:
			model = gp.Model("fva_model_" + self.model_name)
		react_flux = model.addMVar(shape=self.S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux",
		                           lb=copy.deepcopy(self.lb),
		                           ub=copy.deepcopy(self.ub))
		model.addConstr(self.S @ react_flux == self.b, name="c")
		# obj = np.zeros((1, self.S.shape[1]))
		
		# model.setObjective(obj @ react_flux, GRB.MAXIMIZE)
		model.Params.LogToConsole = 0
		model.Params.ScaleFlag = 0
		model.Params.NumericFocus = 3
		model.Params.Method = 1
		# model.Params.Presolve = 0
		# model.Params.Aggregate = 0
		model.Params.FeasibilityTol = 1e-9
		
		index_size = np.size(self.lb)
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
			if self.lb[reaction_index_to_check] > 0 or self.ub[reaction_index_to_check] < 0:
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
				if not (self.lb[reaction_index_to_check] > 0 or self.ub[reaction_index_to_check] < 0):
					true_lb[reaction_index_to_check] = save_lb
					true_ub[reaction_index_to_check] = save_ub
				if list_to_return == "index":
					essential_reactions.append(reaction_index_to_check)
				elif list_to_return == "names":
					essential_reactions.append(self.reaction_names[reaction_index_to_check])
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
			react_flux = model.addMVar(shape=self.S.shape[1], vtype=GRB.CONTINUOUS, name="react_flux", lb=self.lb,
			                           ub=self.ub)
			model.setObjective(react_flux @ Q @ react_flux - c @ react_flux + amount_to_add_back, GRB.MINIMIZE)
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
	
	def pinch_restricted_exchange_reactions(self, allowed_exchange_reaction_file, restore_essential=True,true_exchange_signifier = ""):
		# Leave option to keep excrete free
		met_rxn_dict = self.get_exchanged_reaction_info()
		reactions_to_allow = list(np.loadtxt(allowed_exchange_reaction_file, dtype="str", delimiter=",")[1:, 0])
		reactions_to_allow = [i.replace(" ", "") for i in reactions_to_allow]
		exchange_type = list(np.loadtxt(allowed_exchange_reaction_file, dtype="str", delimiter=",")[1:, 1])
		exchange_type = [i.replace(" ", "") for i in exchange_type]
		for rxn in met_rxn_dict.values():
			if rxn not in reactions_to_allow and true_exchange_signifier in rxn:
				print(rxn)
				# right now this still allows excretion
				if restore_essential:
					rxn_index = self.reaction_names.index(rxn)
					ub_save = self.ub[rxn_index]
					lb_save = self.lb[rxn_index]
				self.update_reaction_bounds(rxn, 0, 0)
				if restore_essential:
					if not self.test_feasibility():
						print(f"removing reaction {rxn} breaks feasibility")
						self.update_reaction_bounds(rxn, lb_save, ub_save)
			elif rxn in reactions_to_allow:
				rxn_index = reactions_to_allow.index(rxn)
				rxn_exchange_type = exchange_type[rxn_index]
				if rxn_exchange_type == "intake":
					self.update_reaction_bounds(rxn, "keep", 0)
				elif rxn_exchange_type == "excrete":
					self.update_reaction_bounds(rxn, 0, "keep")
				elif rxn_exchange_type != "both":
					print("Exchange Reaction direction not understood")
	
	def metabolite_search(self, formula, require_whole_match=True, elements_to_ignore=[]):
		# Currently a WIP
		
		formula_dict = {}
		while len(formula) > 0:
			ph = 14
		possible_metabolites = {}
		for i in range(len(self.metabolite_names)):
			d = 2
		
		return possible_metabolites
	
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
	
	def test_essential_list(self, reaction_index_list, use_gp_env=True):
		# test against set too
		if isinstance(reaction_index_list[0], str):
			reaction_index_list = [self.reaction_names.index(name) for name in reaction_index_list]
		if len(set(reaction_index_list)) != len(reaction_index_list):
			return False, "Duplicate Essential Flux"
		save_ub_full = copy.deepcopy(self.ub)
		save_lb_full = copy.deepcopy(self.lb)
		essential_flux_list = np.zeros(len(self.reaction_names))
		for rxn_ind in reaction_index_list:
			essential_flux_list[rxn_ind] = 1
		for essential_matrix_index in range(len(essential_flux_list)):
			# if the reaction is essential then lb and ub will be the same, if it is not then lb and ub will be 0
			self.update_reaction_bounds(essential_matrix_index,
			                            essential_flux_list[essential_matrix_index] * self.lb[essential_matrix_index],
			                            essential_flux_list[essential_matrix_index] * self.ub[essential_matrix_index])
		if not self.test_feasibility():
			return False, "not feasible enough"
		for essential_matrix_index in range(len(essential_flux_list)):
			if essential_matrix_index % 1000 == 0:
				print(essential_matrix_index / len(essential_flux_list))
			if essential_flux_list[essential_matrix_index] == 1:
				if not (self.lb[essential_matrix_index] > 0 or self.ub[essential_matrix_index] < 0):
					save_ub = copy.deepcopy(self.ub[essential_matrix_index])
					save_lb = copy.deepcopy(self.lb[essential_matrix_index])
					self.update_reaction_bounds(essential_matrix_index, 0, 0)
					if self.test_feasibility():
						return False, "Too feasible"
					else:
						self.update_reaction_bounds(essential_matrix_index, save_lb, save_ub)
		for i in range(len(self.ub)):
			self.update_reaction_bounds(i, save_lb_full[i], save_ub_full[i])
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
		seed_empty = False
		if seed_list == []:
			seed_list = [random.random() for i in range(n)]
			seed_empty = True
		if seed_name:
			file_name_base = file_name_base+str(seed_list[0]).replace(".","")
		value_list = np.zeros((n, len(self.reaction_names)))
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
			np.save(output_dir/file_name_base,pd.DataFrame(value_list).to_numpy())
	
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

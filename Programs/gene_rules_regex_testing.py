import re
import numpy as np


def size_order_key(e):
	return len(e)


def get_reaction_activity_scores(string_to_handle,gene_value_dictionary):
	# additional arguements will be needed for different notation systems
	# gene value dictionary should be in the form "gene name" and "abundance"
	string_to_handle = string_to_handle.replace("(", "{").replace(")", "}")
	vart = True
	while vart:
		m = re.search('{([^{}]+)}', string_to_handle)
		if m is not None:
			if "or" in m.group(0):
				new_str = m.group(0).replace("or", "+").replace("{", "(").replace("}", ")")
			if "and" in m.group(0):
				new_str = m.group(0).replace("and", ",").replace("{", "min(").replace("}", ")")
			string_to_handle = string_to_handle.replace(m.group(0), new_str)
		else:
			if "or" in string_to_handle:
				string_to_handle = "(" + string_to_handle.replace("or", "+") + ")"
			if "and" in string_to_handle:
				string_to_handle = "min(" + string_to_handle.replace("and", ",") + ")"
			vart = False
	gene_list = list(gene_value_dictionary.keys())

	gene_list.sort(reverse=True, key=size_order_key)
	
	for i in gene_list:
		string_to_handle = string_to_handle.replace(i, str(gene_value_dictionary[i]))
	
	return eval(string_to_handle)


'(HGNC:10606 and HGNC:121 and HGNC:2754 and (HGNC:3247 or HGNC:5213)) or (HGNC:121 and HGNC:2754 and (HGNC:3247 or HGNC:5213) and HGNC:82)'

for i in range(10000):
	rand_pick_gene = np.round(np.random.rand(16) * 100).astype(int)
	gene_values = {}
	gene_values["HGNC:10606"] = rand_pick_gene[0]
	gene_values["HGNC:121"] = rand_pick_gene[1]
	gene_values["HGNC:2754"] = rand_pick_gene[2]
	gene_values["HGNC:3247"] = rand_pick_gene[3]
	gene_values["HGNC:5213"] = rand_pick_gene[4]
	gene_values["HGNC:82"] = rand_pick_gene[5]
	
	gene_values["HGNC:4801"] = rand_pick_gene[6]
	gene_values["HGNC:4803"] = rand_pick_gene[7]
	gene_values["HGNC:88"] = rand_pick_gene[8]
	gene_values["HGNC:89"] = rand_pick_gene[9]
	gene_values["HGNC:90"] = rand_pick_gene[10]
	gene_values["HGNC:91"] = rand_pick_gene[11]
	
	gene_values["HGNC:42"] = rand_pick_gene[12]
	gene_values["HGNC:53"] = rand_pick_gene[13]
	gene_values["HGNC:74"] = rand_pick_gene[14]
	
	gene_values["HGNC:40"] = rand_pick_gene[15]
	
	ft = min(gene_values["HGNC:10606"], gene_values["HGNC:121"], gene_values["HGNC:2754"],
	         gene_values["HGNC:3247"] + gene_values["HGNC:5213"])
	st = min(gene_values["HGNC:121"], gene_values["HGNC:2754"], gene_values["HGNC:3247"] + gene_values["HGNC:5213"],
	         gene_values["HGNC:82"])
	expected_equation_1 = ft + st
	
	ft = min(gene_values["HGNC:4801"], gene_values["HGNC:4803"], gene_values["HGNC:88"])
	ft += min(gene_values["HGNC:4801"], gene_values["HGNC:4803"], gene_values["HGNC:89"])
	ft += min(gene_values["HGNC:4801"], gene_values["HGNC:4803"], gene_values["HGNC:90"])
	ft += min(gene_values["HGNC:4801"], gene_values["HGNC:4803"], gene_values["HGNC:91"])
	expected_equation_2 = ft
	
	expected_equation_3 = gene_values["HGNC:42"] + gene_values["HGNC:53"] + gene_values["HGNC:74"]
	
	expected_equation_4 = gene_values["HGNC:40"]
	
	sol_list = [expected_equation_1, expected_equation_2, expected_equation_3,expected_equation_4]
	string_to_handle_1 = '(HGNC:10606 and HGNC:121 and HGNC:2754 and (HGNC:3247 or HGNC:5213)) or (HGNC:121 and HGNC:2754 and (HGNC:3247 or HGNC:5213) and HGNC:82)'
	string_to_handle_2 = "(HGNC:4801 and HGNC:4803 and HGNC:88) or (HGNC:4801 and HGNC:4803 and HGNC:89) or (HGNC:4801 and HGNC:4803 and HGNC:90) or (HGNC:4801 and HGNC:4803 and HGNC:91)"
	string_to_handle_3 = "HGNC:42 or HGNC:53 or HGNC:74"
	
	
	
	string_list = [string_to_handle_1, string_to_handle_2,string_to_handle_3]
	for string_to_handle_p in range(len(string_list)):
		string_to_handle = string_list[string_to_handle_p].replace("(", "{").replace(")", "}")
		string_to_handle = get_reaction_activity_scores(string_to_handle,gene_values)
		#print(string_to_handle)
		print(string_to_handle)

		if string_to_handle != sol_list[string_to_handle_p]:
			print(string_to_handle, sol_list[string_to_handle_p])
			input()

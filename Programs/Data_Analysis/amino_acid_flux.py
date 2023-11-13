import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

os.chdir("../..")
print(os.getcwd())
file_array = np.loadtxt("Data/Input/Experimental_Data/Amino_Acid_Raw/Dixit_Amino Acids_Sum.csv",delimiter=",",dtype=str)
row_header = list(file_array[1:,0])
mouse_names = []
days = []
replicate_numbers = []
for i in row_header:
	if i.split(" ")[0] not in mouse_names:
		mouse_names.append(i.split(" ")[0])
	if i.split(" ")[1] not in replicate_numbers:
		replicate_numbers.append(i.split(" ")[1])
	if i.split(" ")[3] not in days:
		days.append(i.split(" ")[3])
		
	
print(mouse_names)
column_header = list(file_array[0,1:])
data_array_string = file_array[1:,1:]
data_array = np.ones((np.shape(data_array_string)))
for i in range(np.shape(data_array)[0]):
	for j in range(np.shape(data_array)[1]):
		if data_array_string[i,j].replace(".","",1).isdigit():
			data_array[i,j] = float(data_array_string[i,j])
		else:
			data_array[i, j] = np.nan

unique_row_headers = []

for i in row_header:
	if i in unique_row_headers:
		continue
	else:
		unique_row_headers.append(i)

list_of_header_dicts = []
for i in unique_row_headers:
	header_dict = {}
	header_dict["m_type"] = i.split(" ")[0]
	header_dict["r_number"] = i.split(" ")[1]
	header_dict["day"] = i.split(" ")[3]
	list_of_header_dicts.append(header_dict)

sample_ID_arrays = [data_array[4*i:4*(i+1),:] for i in range(len(unique_row_headers))]

list_of_timecourse_results = []
list_of_timecourse_names = []
for m_type in mouse_names:
	for replicate_number in replicate_numbers:
		timecourse_array = np.zeros((len(days),len(column_header)))
		for day_index in range(len(days)):
			for header_index in range(len(list_of_header_dicts)):
				if list_of_header_dicts[header_index]["m_type"] == m_type and list_of_header_dicts[header_index]["r_number"] == replicate_number and list_of_header_dicts[header_index]["day"] == days[day_index]:
					timecourse_array[day_index] = np.nanmean(sample_ID_arrays[header_index],axis=0)
		if np.sum(timecourse_array) != 0:
			list_of_timecourse_results.append(timecourse_array)
			list_of_timecourse_names.append({"m_type":m_type,"r_number":replicate_number})
print(list_of_timecourse_results)
print(list_of_timecourse_names)
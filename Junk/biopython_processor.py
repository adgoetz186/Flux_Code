import sys
import os
import json


print(os.listdir())
with open("biopython_retrived_annotation.json","r") as bpra:
	data_dict = json.load(bpra)
print(data_dict["DocumentSummarySet"]["DocumentSummary"][0].keys())
print(data_dict["DocumentSummarySet"]["DocumentSummary"][0]["OtherAliases"])
print(data_dict["DocumentSummarySet"]["DocumentSummary"][1]["OtherAliases"])
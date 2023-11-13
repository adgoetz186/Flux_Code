import numpy as np
import gzip

f=gzip.open("human_chimpanzee_hcop_fifteen_column.txt.gz",'rb')
file_content=f.readlines()
print(str(file_content[0]).split("\\t"))

gene = np.loadtxt("human_mouse_hcop_fifteen_column.txt.gz",dtype=str,delimiter="\t")
print(gene[0])
new_method_list = ['HGNC:249', 'HGNC:250', 'HGNC:663', 'HGNC:1795', 'HGNC:417', 'HGNC:4336', 'HGNC:7850', 'HGNC:8807', 'HGNC:8724', 'HGNC:9463', 'HGNC:25313', 'HGNC:26891', 'HGNC:341', 'HGNC:2279', 'HGNC:7427', 'HGNC:11038', 'HGNC:637', 'HGNC:16270', 'HGNC:11061', 'HGNC:7455', 'HGNC:7456', 'HGNC:7458', 'HGNC:7459', 'HGNC:7460', 'HGNC:7461', 'HGNC:7462', 'HGNC:2287', 'HGNC:7419', 'HGNC:7421', 'HGNC:7422', 'HGNC:10971', 'HGNC:16029']
old_method_list = ['HGNC:249', 'HGNC:250', 'HGNC:663', 'HGNC:1795', 'HGNC:417', 'HGNC:4336', 'HGNC:7850', 'HGNC:8807', 'HGNC:8724', 'HGNC:25313', 'HGNC:26891', 'HGNC:2279', 'HGNC:7427', 'HGNC:11038', 'HGNC:637', 'HGNC:16270', 'HGNC:11061', 'HGNC:7455', 'HGNC:7456', 'HGNC:7458', 'HGNC:7459', 'HGNC:7460', 'HGNC:7461', 'HGNC:7462', 'HGNC:2287', 'HGNC:7419', 'HGNC:7421', 'HGNC:7422', 'HGNC:10971', 'HGNC:16029']

print(len(new_method_list),len(old_method_list))
for i in old_method_list:
	if i not in new_method_list:
		print(i)
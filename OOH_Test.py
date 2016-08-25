import os

SMILES_files = []
for root, dirs, files in os.walk("../RMG-models", topdown=False):
	for name in files:
		if str(os.path.join(root,name)).endswith('SMILES.txt'):
			SMILES_files.append(os.path.join(root,name))


smiles_dict = {}
for name in SMILES_files:
	try:
		e,d,c,b,a = name.split('/')
		smiles_entry = c + "/" + b
	except ValueError:
		d,c,b,a = name.split('/')
		smiles_entry = b
	smiles_dict[name] = smiles_entry

#print smiles_dict

for entry in SMILES_files:
	smiles_file = str(entry)
	print smiles_file
	print smiles_dict[entry]

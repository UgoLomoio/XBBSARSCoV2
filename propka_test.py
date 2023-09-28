import os
from propka.run import single
dir = "Data/PDB/"
pdb_list = [file for file in os.listdir(dir) if file[-4:] == ".pdb"]
for pdb in pdb_list:
	print(pdb)
	output = single(dir+pdb)

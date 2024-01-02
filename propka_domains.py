import os
from propka.run import single
import os 
import numpy as np

spikes = ["WT","Delta","BA.2.75", "XBB.1","XBB.1.16","XBB.1.5","XBB.1.9.1", "EG.5.1"]
domains = {
"NTD": np.arange(14, 305, 1),
"RBD":  np.arange(331, 527, 1)
}
conformations = ["closed", "open", "complex"]


cwd = os.getcwd()
sep = os.sep
pdb_dir = cwd+sep+"PDBs"+sep

for domain, residues in domains.items():
    
    pdb_domain_dir = pdb_dir + sep + domain
    resi_to_leave = "{}-{}".format(residues[0], residues[-1])
    if not os.path.isdir(pdb_domain_dir):
        os.mkdir(pdb_domain_dir)
    
    for conformation in conformations:

        for variant  in spikes:
            print(variant, conformation, domain)
            pdb_domain_file = pdb_domain_dir + sep + variant + "_" + domain + "_" + conformation + ".pdb"
            output = single(pdb_domain_file)

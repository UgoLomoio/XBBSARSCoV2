import pymol 
from pymol import cmd
import os 
import numpy as np

#spikes = ["XBB.1","XBB.1.16","XBB.1.5","XBB.1.9.1", "EG.5.1"]
spikes = ["WT","Delta","BA.2.75", "XBB.1","XBB.1.16","XBB.1.5","XBB.1.9.1", "EG.5.1"]


domains = {
"NTD": np.arange(14, 305, 1),
"RBD":  np.arange(331, 527, 1)
}
conformations = ["closed", "open", "complex"]
cwd = os.getcwd()
sep = os.sep
pdb_dir = cwd+sep+"PDBs"

for domain, residues in domains.items():
    
    pdb_domain_dir = pdb_dir + sep + domain
    resi_to_leave = "{}-{}".format(residues[0], residues[-1])
    if not os.path.isdir(pdb_domain_dir):
        os.mkdir(pdb_domain_dir)
    
    for conformation in conformations:
        for pdb in spikes:
        
            pdb_file = pdb_dir + sep + conformation + sep + pdb +  "_" + conformation + ".pdb"
            pdb_domain_file = pdb_domain_dir + sep + pdb + "_" + domain + "_" + conformation + ".pdb"
            print(pdb_file, conformation, domain)
        
            cmd.do("delete all")
            cmd.load(pdb_file)
            cmd.do("remove chain B+C")
            cmd.do("remove not resi {} and chain A".format(resi_to_leave))
            cmd.do("hide licorice")
            cmd.do("show cartoon")
            cmd.do("save {}, {}".format(pdb_domain_file, pdb))
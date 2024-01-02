import os
from propka.run import single
import os 
import numpy as np

spike_muts_variants = {
                        "XBB.1": np.array([19, 24, 83, 142, 146, 183, 213, 252, 339, 346, 368, 371, 373, 375, 376, 405, 408, 417, 440, 445, 446,
                                              460, 477, 478, 484, 490, 498, 501, 505, 614, 655, 679, 681, 764, 796, 954, 969]),
                        "XBB.1.9.1": np.array([19, 24, 83, 142, 146, 183, 213, 252, 339, 346, 368, 371, 373, 375, 376, 405, 408, 417, 440, 445, 446,
                                              460, 477, 478, 484, 486, 490, 498, 501, 505, 614, 655, 679, 681, 764, 796, 954, 969]),
                       "XBB.2.3": np.array([19, 24, 83, 142, 146, 183, 213, 253, 339, 346, 368, 371, 373, 375, 376, 405, 408, 417, 440, 445, 446,
                                              460, 477, 478, 484, 486, 490, 498, 501, 505, 521, 614, 655, 679, 681, 764, 796, 954, 969]),
                       "XBB.1.5": np.array([19, 24, 83, 142, 146, 183, 213, 252, 339, 346, 368, 371, 373, 375, 376, 405, 408, 417, 440, 445, 446,
                                              460, 477, 478, 484, 486, 490, 498, 501, 505, 614, 655, 679, 681, 764, 796, 954, 969]),
                       "XBB.1.16": np.array([19, 24, 83, 142, 146, 180, 183, 213, 252, 339, 346, 368, 371, 373, 375, 376, 405, 408, 417, 440, 445, 446,
                                              460, 477, 478, 484, 486, 490, 498, 501, 505, 614, 655, 679, 681, 764, 796, 954, 969]),
                       "EG.5.1": np.array([19, 24, 52, 83, 142, 146, 180, 183, 213, 252, 339, 346, 368, 371, 373, 375, 376, 405, 408, 417, 440, 445, 446,
                                            456, 460, 477, 478, 484, 486, 490, 498, 501, 505, 614, 655, 679, 681, 764, 796, 954, 969])
                      }
                      
                      
domains = {
"NTD": np.arange(14, 305, 1),
"RBD":  np.arange(331, 527, 1)
}

variants = list(spike_muts_variants.keys())
spikes = []
conformations = ["closed", "open", "complex"]
for conformation in conformations:
    for variant in variants:
        spikes.append("{}_{}".format(variant, conformation))
cwd = os.getcwd()
sep = os.sep
pdb_dir = cwd

for domain, residues in domains.items():
    
    pdb_domain_dir = pdb_dir + sep + domain
    resi_to_leave = "{}-{}".format(residues[0], residues[-1])
    if not os.path.isdir(pdb_domain_dir):
        os.mkdir(pdb_domain_dir)
    
    for pdb in spikes:
    
        pdb_domain_file = pdb_domain_dir + sep + pdb + "_"+domain+".pdb"
        output = single(pdb_domain_file)

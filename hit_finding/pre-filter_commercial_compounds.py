# filter a list of commercially available compounds

import os, gzip, numpy
from collections import defaultdict
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Crippen

path = os.getcwd() + '/'
inpath = path + '../data/'

# read, filter and write the commercial compounds
count = 0
outfile = gzip.open(path+'commercial_cmps_cleaned.dat.gz', 'w')
outfile.write("#Identifier\tSMILES\n")
for line in gzip.open(inpath+'parent.smi.gz', 'r'):
    if line[0] == "#": continue
    line = line.rstrip().split()
    # contains: [smiles, identifier]
    m = Chem.MolFromSmiles(line[0])
    if m is None: continue
    # number of heavy atoms
    num_ha = m.GetNumHeavyAtoms()
    if num_ha < 15 or num_ha > 50: continue
    # molecular weight
    mw = AllChem.CalcExactMolWt(m)
    if mw < 200 or mw > 700: continue
    # number of rotatable bonds
    num_rb = AllChem.CalcNumRotatableBonds(m)
    if num_rb > 8: continue
    # number of H-bond donors and acceptors
    num_hba = AllChem.CalcNumHBA(m)
    num_hbd = AllChem.CalcNumHBD(m)
    if num_hba > 10 or num_hbd > 5: continue
    # keep the molecule
    outfile.write("%s\t%s\n" % (line[0], line[1]))
    count += 1
outfile.close()
print "number of molecules that passed the filters:", count

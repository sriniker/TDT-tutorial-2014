# global variables and
# common functions used for training of ML models

import numpy, cPickle
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors as rdmd

# global variables
num_act = 1528
num_dcy = 293606

num_rep = 50
num_percent = 0.1

# fingerprint dictionary
fp_dict = {}
fp_dict['morgan2'] = lambda m: rdmd.GetMorganFingerprintAsBitVect(m, 2, nBits=1024)
fp_dict['ap'] = lambda m: rdmd.GetHashedAtomPairFingerprintAsBitVect(m, nBits=2048)
fp_dict['rdk5'] = lambda m: Chem.RDKFingerprint(m, maxPath=5, fpSize=2048, nBitsPerHash=2)

def getNumpyFP(smiles, fpname, fptype):
    m = Chem.MolFromSmiles(smiles)
    if m is not None:
        # calculate fingerprint
        fp = fp_dict[fpname](m)
        # convert to numpy array
        if fptype == 'bool':
            arr = numpy.zeros((1,), numpy.bool)
        elif fptype == 'float':
            arr = numpy.zeros((1,), numpy.float32)
        else:
            raise ValueError('fptype not known', fptype)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
    else:
        return None

def getFP(smiles, fpname):
    m = Chem.MolFromSmiles(smiles)
    if m is not None:
        # calculate fingerprint
        fp = fp_dict[fpname](m)
        return fp
    else:
        return None

def writeActiveRanks(scores, outfile, num_act):
    # scores contains: [proba, simil, info]
    # get indices
    indices = [(p[0], p[1], p[2], i) for i,p in enumerate(scores)]
    # sort by descending probability
    indices.sort(reverse=True)
    # add ranks (highest = best) and sort again by ascending index
    num = len(indices)
    indices = [(p[3], num-i) for i,p in enumerate(indices)]
    indices.sort()
    # keep only the actives
    indices = [r for idx,r in indices[:num_act]]
    # write to file
    cPickle.dump(indices, outfile, 2)

def assignRanksWithInfo(plist): # contains [proba, simil, info]
    num = len(plist)
    tmplist = [(p[0], p[1], p[2], j) for j,p in enumerate(plist)] # add index
    tmplist.sort(reverse=True) # sort by descending probability
    # add the rank: highest = best
    # contains: [index, rank, proba, simil, info]
    tmplist = [(t[3], num-j, t[0], t[1], t[2]) for j,t in enumerate(tmplist)]
    tmplist.sort() # sort by ascending index to return to original order
    # remove the index
    return [(rank, p, s, info) for j,rank,p,s,info in tmplist]

def assignRanks(plist, slist): # input: list of probabilities and list of similarities
    num = len(plist)
    # combine the lists and add the index
    tmplist = [(plist[j], slist[j][0], j) for j in range(num)]
    tmplist.sort(reverse=True) # sort by descending probability
    # add the rank: highest = best
    # contains: [index, rank, proba, simil]
    tmplist = [(t[2], num-j, t[0], t[1]) for j,t in enumerate(tmplist)]
    tmplist.sort() # sort by ascending index to return to original order
    # remove the index
    return [(rank, p, s) for j,rank,p,s in tmplist]



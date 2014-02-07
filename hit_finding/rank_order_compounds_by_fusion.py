# rank orders a list of commercially available compounds
# using a fusion model

import os, gzip, numpy, cPickle, sys
from collections import defaultdict
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors as rdMD

# common functions
sys.path.insert(0, os.getcwd()+'/../')
import common_functions as cf

path = os.getcwd() + '/'
inpath = path + '../data/'

# load the ML models
lr_rdk5 = cPickle.load(gzip.open(path+'../final_models/lr_rdk5_model.pkl.gz', 'r'))
rf_rdk5 = cPickle.load(gzip.open(path+'../final_models/rf_rdk5_model.pkl.gz', 'r'))
rf_morgan2 = cPickle.load(gzip.open(path+'../final_models/rf_morgan2_model.pkl.gz', 'r'))
print "rf models loaded"

# loop over commercial products
proba_lr_rdk5 = []
proba_rf_rdk5 = []
proba_rf_morgan2 = []
mols = []
for line in gzip.open(path+'commercial_cmps_cleaned.dat.gz', 'r'):
    if line[0] == "#": continue
    line = line.rstrip().split()
    # contains: [smiles, identifier]
    # RDK5
    fp = cf.getNumpyFP(line[0], 'rdk5', 'float')
    proba_lr_rdk5.append(lr_rdk5.predict_proba(fp)[0][1])
    proba_rf_rdk5.append(rf_rdk5.predict_proba(fp)[0][1])
    fp = cf.getNumpyFP(line[0], 'morgan2', 'float')
    proba_rf_morgan2.append(rf_morgan2.predict_proba(fp)[0][1])
    mols.append((line[1], line[0]))
print "probabilities calculated"

# load similarities
scores_rdk5 = cPickle.load(gzip.open(path+'scores_rdk5.pkl.gz' , 'r'))
scores_morgan2 = cPickle.load(gzip.open(path+'scores_morgan2.pkl.gz' , 'r'))
"similarities loaded"

# assign ranks
scores_lr_rdk5 = cf.assignRanks(proba_lr_rdk5, scores_rdk5)
scores_rf_rdk5 = cf.assignRanks(proba_rf_rdk5, scores_rdk5)
scores_rf_morgan2 = cf.assignRanks(proba_rf_morgan2, scores_morgan2)
print "ranks assigned"

# fusion
fusion_scores = []
for m1, m2, m3, m in zip(scores_lr_rdk5, scores_rf_rdk5, scores_rf_morgan2, mols):
    rank = max([m1[0], m2[0], m3[0]]) # maximum rank
    pp = max([m1[1], m2[1], m3[1]]) # maximum probability
    # store: [max rank, max proba, similarity, identifier, smiles]
    fusion_scores.append([rank, pp, m1[2], m[0], m[1]])
# sort by descending rank
fusion_scores.sort(reverse=True)
print "fusion done"

# write out
outfile = gzip.open(path+'ranked_list_top10K_commercial_cmps.dat.gz', 'w')
outfile.write("#Identifier\tSMILES\tMax_Rank\tMax_Proba\tSimilarity\n")
for r, pp, s, idx, smiles in fusion_scores[:10000]:
    outfile.write("%s\t%s\t%i\t%.5f\t%.5f\n" % (idx, smiles, r, pp, s))
outfile.close()
print "list written"

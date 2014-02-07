# trains and evaluates a RF model with a given fingerprint

import os, gzip, numpy, cPickle, sys
from collections import defaultdict
from rdkit import Chem, DataStructs
from sklearn.ensemble import RandomForestClassifier, forest
from optparse import OptionParser

# common functions
sys.path.insert(0, os.getcwd()+'/../')
import common_functions as cf

# monkey path for random forest
import random_forest_functions as rfmeth
forest._parallel_build_trees = rfmeth._balanced_parallel_build_trees

path = os.getcwd() + '/'
inpath = path + '../data/'

# prepare command-line option parser
usage = "usage: %prog [options] arg"
parser = OptionParser(usage)
parser.add_option("-f", "--fingerprint", dest="fp", help="fingerprint name (options are: morgan2, ap, rdk5)")

# read in command line options
(options, args) = parser.parse_args()
# required arguments
if options.fp:
    fpname = options.fp
else:
    raise RuntimeError('fingerprint name missing')
print "ML model is trained with", fpname

# read the actives
fps_act = []
for line in open(inpath+'training_actives_cleaned.dat', 'r'):
    line = line.strip().split()
    # contains: [sample_id, hit, pec50, smiles]
    fp = cf.getNumpyFP(line[3], fpname, 'float')
    if fp is not None:
        fps_act.append(fp)
num_actives = len(fps_act)
print "actives read and fingerprints calculated:", num_actives

# read the inactives
fps_inact = []
for line in open(inpath+'training_inactives_cleaned.dat', 'r'):
    line = line.strip().split()
    # contains: [sample_id, hit, pec50, smiles]
    fp = cf.getNumpyFP(line[3], fpname, 'float')
    if fp is not None:
        fps_inact.append(fp)
num_inactives = len(fps_inact)
print "inactives read and fingerprints calculated:", num_inactives

# training 
train_fps = fps_act + fps_inact
ys_fit = [1]*num_actives + [0]*num_inactives
# train the model
ml = RandomForestClassifier(n_estimators=100, max_depth=100, min_samples_split=2, min_samples_leaf=1)
ml.fit(train_fps, ys_fit)

# write the model to file
outfile = gzip.open(path+'rf_'+fpname+'_model.pkl.gz', 'wb')
cPickle.dump(ml, outfile, 2)
outfile.close()

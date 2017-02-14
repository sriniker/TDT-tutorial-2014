# does similarity search (MAX fusion) with a given fingerprint
# and does the prediction for the test molecules
from __future__ import print_function

try:
    import cPickle as pickle
except ImportError:
    import pickle
import gzip
import os
import sys
from optparse import OptionParser

from rdkit import DataStructs

# common functions
sys.path.insert(0, os.getcwd()+'/../')
import common_functions as cf

path = os.getcwd() + '/'
inpath = path + '../data/'

# prepare command-line option parser
usage = "usage: %prog [options] arg"
parser = OptionParser(usage)
parser.add_option("-f", "--fingerprint", dest="fp",
                  help="fingerprint name (options are: %s)" % ', '.join(sorted(cf.fp_dict)))

# read in command line options
(options, args) = parser.parse_args()
# required arguments
if options.fp:
    fpname = options.fp
else:
    raise RuntimeError('fingerprint name missing')
print("similarity search with", fpname)

# read the training actives
fps_act = []
for line in gzip.open(inpath+'training_actives_cleaned.dat.gz', 'rt'):
    if line[0] == "#":
        continue
    line = line.rstrip().split()
    # contains: [sample_id, hit, pec50, smiles]
    fp = cf.getFP(line[3], fpname)
    if fp is not None:
        fps_act.append(fp)
print("training actives read and fingerprints calculated:", len(fps_act))

# read the commercial compounds
scores = []
for line in gzip.open(path+'commercial_cmps_cleaned.dat.gz', 'rt'):
    if line[0] == "#":
        continue
    line = line.rstrip().split()
    # contains: [smiles, identifier]
    fp = cf.getFP(line[0], fpname)
    simil = DataStructs.BulkTanimotoSimilarity(fp, fps_act)
    simil.sort(reverse=True)
    scores.append([simil[0], line[1]])
print("commercial compounds read and fingerprints calculated")

# write the list
rankfile = gzip.open(path+'scores_'+fpname+'.pkl.gz', 'wb')
pickle.dump(scores, rankfile, 2)
rankfile.close()
print("scoring done")

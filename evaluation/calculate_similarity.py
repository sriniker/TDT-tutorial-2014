# does similarity search (MAX fusion) with a given fingerprint
# and does the prediction for the test molecules

import os, gzip, numpy, cPickle, sys
from collections import defaultdict
from rdkit import Chem, DataStructs
from rdkit.ML.Scoring import Scoring
from optparse import OptionParser

# common functions
sys.path.insert(0, os.getcwd()+'/')
sys.path.insert(0, os.getcwd()+'/../')
import common_functions as cf

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
print "similarity search with", fpname

# read the actives
fps_act = []
for line in open(inpath+'training_actives_cleaned.dat', 'r'):
    line = line.rstrip().split()
    # contains: [sample_id, hit, pec50, smiles]
    fp = cf.getFP(line[3], fpname)
    if fp is not None:
        fps_act.append(fp)
print "actives read and fingerprints calculated:", len(fps_act)

# read the inactives
fps_dcy = []
for line in open(inpath+'training_inactives_cleaned.dat', 'r'):
    line = line.rstrip().split()
    # contains: [sample_id, hit, pec50, smiles]
    fp = cf.getFP(line[3], fpname)
    if fp is not None:
        fps_dcy.append(fp)
print "inactives read and fingerprints calculated:", len(fps_dcy)

# test lists
test_indices_act = []
for line in gzip.open(path+'../test_lists/test_lists_10percent_actives.dat.gz', 'r'):
    line = line.rstrip().split()
    test_indices_act.append([int(l) for l in line])
test_indices_dcy = []
for line in gzip.open(path+'../test_lists/test_lists_10percent_inactives.dat.gz', 'r'):
    line = line.rstrip().split()
    test_indices_dcy.append([int(l) for l in line])
print "test lists read in:", len(test_indices_act), len(test_indices_dcy)

# loop over the repetitions
outfile = open(path+'analysis/output_analysis_simil_'+fpname+'.dat', 'w')
rankfile = gzip.open(path+'scores/lists_'+fpname+'.pkl.gz', 'wb')
for i in range(cf.num_rep):
    print "repetition", i

    # divide the fps into training and testing
    train_fps = []
    test_fps = []
    for j in range(cf.num_act): # actives
        if j in test_indices_act[i]:
            test_fps.append([fps_act[j], 1])
        else:
            train_fps.append(fps_act[j])
    for j in range(cf.num_dcy): # inactives
        if j in test_indices_dcy[i]:
            test_fps.append([fps_dcy[j], 0])

    # ranking
    scores = []
    for fp, j in test_fps:
        simil = DataStructs.BulkTanimotoSimilarity(fp, train_fps)
        simil.sort(reverse=True)
        scores.append([simil[0], j])
    # write the list before the final ranking
    cPickle.dump(scores, rankfile, 2)

    scores.sort(reverse=True)

    # evaluation
    auc = Scoring.CalcAUC(scores, -1)
    ef = Scoring.CalcEnrichment(scores, -1, [0.05])

    # write out
    outfile.write("%i\t%.10f\t%.10f\n" % (i, auc, ef[0]))

rankfile.close()
outfile.close()

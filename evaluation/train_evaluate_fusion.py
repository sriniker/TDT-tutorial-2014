# trains and evaluates a fusion model

import os, gzip, numpy, cPickle, sys
from collections import defaultdict
from rdkit import Chem, DataStructs
from rdkit.ML.Scoring import Scoring
from sklearn.ensemble import RandomForestClassifier, forest
from sklearn.linear_model import LogisticRegression

# common functions
sys.path.insert(0, os.getcwd()+'/')
sys.path.insert(0, os.getcwd()+'/../')
import common_functions as cf

# monkey path for random forest
import random_forest_functions as rfmeth
forest._parallel_build_trees = rfmeth._balanced_parallel_build_trees

path = os.getcwd() + '/'
inpath = path + '../data/'

fpname1 = 'rdk5'
fpname2 = 'morgan2'

# read the actives
fps_act_rdk5 = []
fps_act_morgan2 = []
for line in gzip.open(inpath+'training_actives_cleaned.dat.gz', 'r'):
    line = line.rstrip().split()
    # contains: [sample_id, hit, pec50, smiles]
    fp = cf.getNumpyFP(line[3], fpname1, 'float')
    if fp is not None:
        fps_act_rdk5.append(fp)
    fp = cf.getNumpyFP(line[3], fpname2, 'float')
    if fp is not None:
        fps_act_morgan2.append(fp)
print "actives read and fingerprints calculated"
num_actives = len(fps_act_rdk5)

# read the inactives
fps_inact_rdk5 = []
fps_inact_morgan2 = []
for line in gzip.open(inpath+'training_inactives_cleaned.dat.gz', 'r'):
    line = line.rstrip().split()
    # contains: [sample_id, hit, pec50, smiles]
    fp = cf.getNumpyFP(line[3], fpname1, 'float')
    if fp is not None:
        fps_inact_rdk5.append(fp)
    fp = cf.getNumpyFP(line[3], fpname2, 'float')
    if fp is not None:
        fps_inact_morgan2.append(fp)
print "inactives read and fingerprints calculated"
num_inactives = len(fps_inact_rdk5)

# test lists
test_indices_act = []
for line in gzip.open(path+'../test_lists/test_lists_10percent_actives.dat.gz', 'r'):
    line = line.rstrip().split()
    test_indices_act.append(set([int(l) for l in line]))
test_indices_inact = []
for line in gzip.open(path+'../test_lists/test_lists_10percent_inactives.dat.gz', 'r'):
    line = line.rstrip().split()
    test_indices_inact.append(set([int(l) for l in line]))
print "test lists read in:", len(test_indices_act), len(test_indices_inact)

# loop over the repetitions
infile1 = gzip.open(path+'scores/lists_'+fpname1+'.pkl.gz' , 'r')
infile2 = gzip.open(path+'scores/lists_'+fpname2+'.pkl.gz' , 'r')
outfile = open(path+'analysis/output_analysis_fusion.dat', 'w')
for i in range(cf.num_rep):
    print "repetition", i
    # indices of training molecules
    train_indices_act = set(range(num_actives)).difference(test_indices_act[i])
    train_indices_inact = set(range(num_inactives)).difference(test_indices_inact[i])
    ys_fit = [1]*len(train_indices_act) + [0]*len(train_indices_inact)

    ### Models with RDK5 
    train_fps = [fps_act_rdk5[j] for j in train_indices_act]
    train_fps += [fps_inact_rdk5[j] for j in train_indices_inact]
    # train the RF
    rf_rdk5 = RandomForestClassifier(n_estimators=100, max_depth=100, min_samples_split=2, min_samples_leaf=1)
    rf_rdk5.fit(train_fps, ys_fit)
    # train the LR
    lr_rdk5 = LogisticRegression()
    lr_rdk5.fit(train_fps, ys_fit)
    # chemical similarity
    simil = cPickle.load(infile1)
    # ranking
    test_fps = [fps_act_rdk5[j] for j in test_indices_act[i]] 
    test_fps += [fps_inact_rdk5[j] for j in test_indices_inact[i]]
    scores_rf_rdk5 = [[pp[1], s[0], s[1]] for pp,s in zip(rf_rdk5.predict_proba(test_fps), simil)]
    scores_lr_rdk5 = [[pp[1], s[0], s[1]] for pp,s in zip(lr_rdk5.predict_proba(test_fps), simil)]

    ### Model with Morgan2
    train_fps = [fps_act_morgan2[j] for j in train_indices_act]
    train_fps += [fps_inact_morgan2[j] for j in train_indices_inact]
    # train the RF
    rf_morgan2 = RandomForestClassifier(n_estimators=100, max_depth=100, min_samples_split=2, min_samples_leaf=1)
    rf_morgan2.fit(train_fps, ys_fit)
    # chemical similarity
    simil = cPickle.load(infile2)
    # ranking
    test_fps = [fps_act_morgan2[j] for j in test_indices_act[i]] 
    test_fps += [fps_inact_morgan2[j] for j in test_indices_inact[i]]
    scores_rf_morgan2 = [[pp[1], s[0], s[1]] for pp,s in zip(rf_morgan2.predict_proba(test_fps), simil)]

    # assign ranks
    scores_rf_rdk5 = cf.assignRanksWithInfo(scores_rf_rdk5)
    scores_lr_rdk5 = cf.assignRanksWithInfo(scores_lr_rdk5)
    scores_rf_morgan2 = cf.assignRanksWithInfo(scores_rf_morgan2)

    # fusion
    fusion_scores = []
    for m1,m2,m3 in zip(scores_rf_rdk5, scores_lr_rdk5, scores_rf_morgan2):
        rank = max([m1[0], m2[0], m3[0]]) # max. rank
        proba = max([m1[1], m2[1], m3[1]]) # max. rank
        # store: [max rank, max proba, simil, info]
        fusion_scores.append([rank, proba, m1[2], m1[3]])
    fusion_scores.sort(reverse=True)

    # evaluation
    auc = Scoring.CalcAUC(fusion_scores, -1)
    ef = Scoring.CalcEnrichment(fusion_scores, -1, [0.05])

    # write out
    outfile.write("%i\t%.10f\t%.10f\n" % (i, auc, ef[0]))

infile1.close()
infile2.close()
outfile.close()

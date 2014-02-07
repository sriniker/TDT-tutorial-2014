# read the ranks of the actives from different models
# and calculates the Pearson correlation coefficient

import os, gzip, numpy, cPickle, sys
from collections import defaultdict
from rdkit import Chem, DataStructs
from scipy.stats import pearsonr

# common functions
sys.path.insert(0, os.getcwd()+'/../')
import common_functions as cf

path = os.getcwd() + '/'
inpath = path + '../data/'

#methods = ['lr', 'rf', 'nb']
methods = ['rf', 'lr']
fpnames = ['ap', 'rdk5', 'morgan2']

# read the ranks
ranks = defaultdict(list)
names = []
for m in methods:
    for fp in fpnames:
        name = m+'_'+fp
        print name
        names.append(name)
        infile = gzip.open(path+'scores/ranks_actives_'+name+'.pkl.gz', 'r')
        for i in range(cf.num_rep):
            ranks[name].append(cPickle.load(infile))
        infile.close()
print "number of models =", len(ranks)

# loop over the models
num_names = len(names)
names2 = []
ave_pearson = {}
for i in range(num_names):
    n1 = names[i]
    for j in range(i+1, num_names):
        n2 = names[j]
        pearson = []
        for r in range(cf.num_rep):
            tmp = pearsonr(ranks[n1][r], ranks[n2][r]) # [pearson, p-value]
            pearson.append(tmp[0]) # keep only pearson coefficient
        ave_pearson[n1+'-'+n2] = [numpy.average(pearson), numpy.std(pearson)]
        names2.append(n1+'-'+n2)

# write out
outfile = open(path+'analysis/pearson_correlation.dat', 'w')
outfile.write("#Model1-Model2\tAve_Pearson\tSTD\n")
outfile.writelines("%s\t%.3f\t%.3f\n" % (k, ave_pearson[k][0], ave_pearson[k][1]) for k in names2)
outfile.close()

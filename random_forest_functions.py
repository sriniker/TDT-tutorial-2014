from rdkit.ML.Data import DataUtils
import numpy
from sklearn.ensemble import RandomForestClassifier, forest
from sklearn.tree import tree

### FOR SKLEARN VERSION 0.13 ###

# HELPER FUNCTIONS FOR RANDOM FOREST
def _balanced_parallel_build_trees(n_trees, forest, X, y, sample_weight, sample_mask, X_argsorted, seed, verbose):
    """Private function used to build a batch of trees within a job"""
    from sklearn.utils import check_random_state
    from sklearn.utils.fixes import bincount
    import random
    MAX_INT = numpy.iinfo(numpy.int32).max
    random_state = check_random_state(seed)

    trees = []
    for i in xrange(n_trees):
        if verbose > 1:
            print("building tree %d of %d" % (i+1, n_trees))
        seed = random_state.randint(MAX_INT)

        tree = forest._make_estimator(append = False)
        tree.set_params(compute_importances=forest.compute_importances)
        tree.set_params(random_state = check_random_state(seed))

        if forest.bootstrap:
            n_samples = X.shape[0]
            if sample_weight is None:
                curr_sample_weight = numpy.ones((n_samples,), dtype=numpy.float64)
            else:
                curr_sample_weight = sample_weight.copy()

            ty = list(enumerate(y))
            indices = DataUtils.FilterData(ty, val=1, frac=0.5, col=1, indicesToUse=0, indicesOnly=1)[0]
            indices2 = random_state.randint(0, len(indices), len(indices))
            indices = [indices[j] for j in indices2]
            sample_counts = bincount(indices, minlength=n_samples)

            curr_sample_weight *= sample_counts
            curr_sample_mask = sample_mask.copy()
            curr_sample_mask[sample_counts==0] = False

            tree.fit(X, y, sample_weight=curr_sample_weight, sample_mask=curr_sample_mask, X_argsorted=X_argsorted, check_input=False)
            tree.indices = curr_sample_mask
        else:
            tree.fit(X, y, sample_weight=sample_weight, sample_mask=sample_mask, X_argsorted=X_argsorted, check_input=False)
        trees.append(tree)
    return trees

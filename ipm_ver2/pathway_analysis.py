import numpy as np

# Functions:
#       cluster_pathways
#       score_align
#       make_clusters
#       get_solution_subset
#       load_pathway_scores
#       filter_good_scoring_solutions

def cluster_pathways(solutions, threshold=0.2):
    from alignment.sequence import Sequence
    from alignment.vocabulary import Vocabulary
    from alignment.sequencealigner import SimpleScoring, GlobalSequenceAligner
    import scipy.spatial.distance as distance
    import scipy.cluster.hierarchy as sch


    def score_align(x, y):
        a = Sequence(x)
        b = Sequence(y)
        v = Vocabulary()
        aEncoded = v.encodeSequence(a)
        bEncoded = v.encodeSequence(b)
        scoring = SimpleScoring(2, -1)
        aligner = GlobalSequenceAligner(scoring, -2)
        score, encodeds = aligner.align(aEncoded, bEncoded, backtrace=True)
        pI = 0.0
        for e in encodeds:
            alignment = v.decodeSequenceAlignment(e)
            pI = max(pI, alignment.percentIdentity())
        return 1 - pI/100.0


    def make_clusters(pathstrs, threshold=0.2):
        distances = np.zeros((len(pathstrs), len(pathstrs)))
        for i in range(len(pathstrs)):
            for j in range(len(pathstrs)):
                distances[i][j] = score_align(pathstrs[i], pathstrs[j])
        
        z = sch.distance.squareform(distances)
        
        clusters = sch.linkage(z, method='single')
        flatclusters = sch.fcluster(clusters, threshold, criterion='distance')
        return flatclusters


    splitsols = [s.split(' -> ') for s in solutions['strrepr']]
    #splitsols = [s.split(' -> ') for s in solutions]
    flatcl = make_clusters(splitsols, threshold=threshold)
    flatcl_dict = {}
    for i in range(len(flatcl)):
        flatcl_dict[i] = flatcl[i]
    porder = sorted(range(len(splitsols)), key=lambda j: flatcl_dict[j])
    #sub_groups = {k:[] for k in range(1, 1+max(flatcl))}
    sub_groups = {}
    for k in range(1, 1+max(flatcl)):
        sub_groups[k] = []
    reprclusters = {}
    for i in porder:
        path = solutions['strrepr'][i]
        clustergroup = flatcl_dict[i]
        sub_groups[clustergroup].append(path)
    for j in range(1, 1+max(flatcl)):
        clustersubset = solutions[flatcl==j]
        maxidx = np.argmax(clustersubset['obj'])
        reprclusters[j] = clustersubset[maxidx]['strrepr']
    print '%d paths clustered into %d groups' % (len(splitsols), max(flatcl))
    return sub_groups, reprclusters


def get_solution_subset(ref_proteins, solutions, keepsub=False):
    subset_solutions = []
    for s in solutions:
        components = s.split(' -> ')
        pathway_proteins = components[1::2]
        ok = True
        for pp in pathway_proteins:
            if pp not in ref_proteins:
                ok = False
        if ok:
            subset_solutions.append(s)
            
    # Remove sub-solutions
    if not keepsub:
        count = 0
        sublist = []
        for s in range(len(subset_solutions)):
            for r in range(s, len(subset_solutions)):
                if subset_solutions[s] in subset_solutions[r] and subset_solutions[s] != subset_solutions[r]:
                    count += 1
                    sublist.append(subset_solutions[s])
                    break
        for subsol in sublist:
            subset_solutions.remove(subsol)
    return subset_solutions

def load_pathway_scores(h5out):
    from tables import open_file
    
    path = np.dtype([('obj', np.float64), ('strrepr', np.str_, 200)])
    with open_file(h5out, 'r') as h5handle:
        table = h5handle.root.paths.pathTable
        prs = [(row['obj'], row['strrepr']) for row in table.iterrows()]
    np_path_array = np.array(prs, dtype=path)
    np_path_array = np.sort(np_path_array, order='obj')[::-1]
    return np_path_array

def filter_good_scoring_solutions(patharray, num_stdevs=1.):
    st = np.std(patharray['obj'])
    maxscore = np.max(patharray['obj'])
    cu = maxscore - num_stdevs*st
    good_solutions = patharray[patharray['obj'] >= cu]
    print '%d paths filtered down to %d within %.1f standard deviation(s)' % (len(patharray), len(good_solutions), num_stdevs)
    good_solutions = np.sort(good_solutions, order='obj')[::-1]
    return good_solutions


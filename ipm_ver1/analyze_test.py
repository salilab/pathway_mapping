import numpy as np
import pickle
import random as rand
import scipy.spatial.distance as distance
import scipy.cluster.hierarchy as sch
from tables import open_file
import time
#import get_uniq_paths
from itertools import islice
from itertools import cycle
import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite
import colorsys
import matplotlib.pyplot as plt
#plt.rcParams['text.usetex'] = True
#plt.rcParams['font.size'] = 18
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib.collections import PatchCollection

def roundrobin(*iterables):
    "roundrobin('ABC', 'D', 'EF') --> A D E B F C"
    # Recipe credited to George Sakkis
    pending = len(iterables)
    nexts = cycle(iter(it).next for it in iterables)
    while pending:
        try:
            for next in nexts:
                yield next()
        except StopIteration:
            pending -= 1 
            nexts = cycle(islice(nexts, pending))


pathdt = np.dtype([('obj', np.float32), ('strrepr', np.str_, 200), ('run', np.int32)])

class pathens(object):

    def __init__(self, datafile=None, fileprefix=None):
        self.enzymes = []
        self.true_enzymes = []
        self.ligands = []

        self.uniqpathstrs = []
        self.best_path = []
        self.best_score = None
        self.pathdict = {}
        self.np_path_array = None # All paths stored as a numpy array
                                  # with objective score and string representation stored
                                  # and run number
        self.flatclusters = None
        self.reprclusters = None
        self.num_runs = None
        self.cutoff = None

        self.conv_test_means = None
        self.conv_test_stds = None

        self.datafile = datafile
        self.fileprefix = fileprefix

        self.interaction_df = None
        self.substrate_interaction_df = None
        self.product_interaction_df = None
        self.ligands_df = None
        self.enzyme_df = None

        self.true_path = ''


    def set_datafile(self, datafile):
        self.datafile = datafile

    def set_fileprefix(self, fileprefix):
        self.fileprefix = fileprefix

    def set_all_enzymes(self, enzymelist):
        self.enzymes = enzymelist

    def set_true_enzymes(self, enzymelist):
        self.true_enzymes = enzymelist

    def set_true_ligands(self, ligandlist):
        self.ligands = ligandlist

    def set_true_path(self, enzymelist=[], ligandlist=[], startposition=0):
        if ligandlist:
            self.ligands = ligandlist
        if enzymelist:
            self.true_enzymes = enzymelist
        true_path = [None]*(len(self.ligands)+len(self.true_enzymes))
        true_path[startposition::2] = self.ligands
        true_path[1-startposition::2] = self.true_enzymes
        self.true_path = ' -> '.join(true_path)

    def filter_by_top_paths(self, topnum=10):
        self.np_path_array.sort(order='obj')
        self.np_path_array = self.np_path_array[::-1][0:topnum]
        self.uniqpathstrs = [ps.split(' -> ') for ps in self.np_path_array['strrepr']]

    def load_data_from_analysis_pickle(self, picklefile):
        with open(picklefile, 'r') as handle:
            pens = pickle.load(handle)
            self.uniqpathstrs = pens.uniqpathstrs
            self.best_path = pens.best_path
            self.best_score = pens.best_score
            self.np_path_array = pens.np_path_array

            self.flatclusters = pens.flatclusters
            self.reprclusters = pens.reprclusters
            self.conv_test_means = pens.conv_test_means
            self.conv_test_stds = pens.conv_test_stds


    def get_uniq_counts(self, runs, cutoff):
        #trs = []
        sofar = {}
        counts = []
        for run in runs:
            h5file = '%s%d.h5' % (self.fileprefix, run)

            print run
            with open_file(h5file, 'r') as h5handle:
                table = h5handle.root.paths.pathTable
                #condition = 'obj2 < %.5f' % cutoff
                condition = 'obj >= %.5f' % cutoff
                for row in table.where(condition):
                    pathstring = row['strrepr']
                    if not sofar.get(pathstring):
                        sofar[pathstring] = 1
                counts.append(len(sofar))
            #print len(trs)
        return counts

    def get_uniq_paths(self, runs, cutoff):
        trs = []
        uniqpathstrs = []
        sofar = {}
        for run in runs:
            runsofar = {}
            rununiqpaths = []
            rappend = rununiqpaths.append
            print 'Run %d' % run
            h5file = '%s%d.h5' % (self.fileprefix, run)
            try:
                with open_file(h5file, 'r') as h5handle:
                    table = h5handle.root.paths.pathTable
                    condition = 'obj >= %.5f' % cutoff
                    for row in table.where(condition):
                        pathstring = row['strrepr']
                        if not runsofar.get(pathstring):
                            strlist = pathstring.split(' -> ')
                            strlist = [s.strip() for s in strlist]
                            rappend(strlist)
                            runsofar[pathstring] = 1

                            if not sofar.get(pathstring):
                                sofar[pathstring] = 1
                                uniqpathstrs.append(strlist)
                                trs.append((row['obj'], row['strrepr'], run))
                    self.pathdict[run] = rununiqpaths
            except:
                print '%d FILE NOT FOUND/RUN INCOMPLETE' % run
        self.uniqpathstrs = uniqpathstrs
        npar = np.array(trs, dtype=pathdt)
        self.np_path_array = npar
        maxidx = np.argmax(self.np_path_array['obj'])
        self.best_path = self.np_path_array['strrepr'][maxidx]
        self.best_score = self.np_path_array['obj'][maxidx]
        self.num_runs = np.max(runs)
        self.cutoff = cutoff


        #########################################################

        #t, u, p = get_uniq_paths.get_uniq_paths(runs, cutoff, self.fileprefix)
        #self.uniqpathstrs = u
        #npar = np.array(t, dtype=pathdt)
        #self.pathdict = p
        #
        #self.np_path_array = npar
        #maxidx = np.argmax(self.np_path_array['obj'])
        #self.best_path = self.np_path_array['strrepr'][maxidx]
        #self.best_score = self.np_path_array['obj'][maxidx]
        #self.num_runs = np.max(runs)
        #self.cutoff = cutoff

        #########################################################

    def count_clusters(self, runs, cutoff):
        uniqpaths = []
        sofar = {}
        countdicts = []
        for run in runs:
            uniqpaths = uniqpaths + filter(lambda x : x not in uniqpaths, self.pathdict[run])
            countdict = count_paths(uniqpaths, self.flatclusters, self.uniqpathstrs)
            countdicts.append(countdict)
        return countdicts

    def make_clusters(self, threshold=0.2, pickleit=False):
    #    f, r = get_uniq_paths.make_clusters(self.uniqpathstrs, self.np_path_array, threshold=threshold)
    #    self.reprclusters = r
    #    self.flatclusters = f
        print 'Calculating distance matrix'
        distances = np.zeros((len(self.uniqpathstrs), len(self.uniqpathstrs)))
        for i in range(len(self.uniqpathstrs)):
            for j in range(len(self.uniqpathstrs)):
                distances[i][j] = compare_lists(self.uniqpathstrs[i], self.uniqpathstrs[j])
        print 'Average distance: %.3f' % np.mean(distances)
        print 'Minimum distance: %.3f' % np.min(distances)
        print 'Maximum distance: %.3f' % np.max(distances)
        print 'Standard deviation distance: %.3f' % np.std(distances)
        print 'Clustering paths'
        z = sch.distance.squareform(distances)
        clusters = sch.linkage(z,  method='single')
        self.flatclusters = sch.fcluster(clusters, threshold, criterion='distance')
        self.reprclusters = {}
        for i in range(1, np.max(self.flatclusters)+1):
            clustersubset = self.np_path_array[self.flatclusters==i]
            maxidx = np.argmax(clustersubset['obj'])
            self.reprclusters[i] = clustersubset[maxidx]['strrepr']

    def simple_cluster(self, threshold=0.3):
        reprclusters = {}
        reprscores = {}
        clusters = []
        totalclnum = 0
        
        for up in np.sort(self.np_path_array, order='obj')[::-1]:
            uniqpath = up['strrepr']
            added = False
            for clnum in reprclusters.keys():
                reprpath = reprclusters[clnum]
                if compare_lists(uniqpath, reprpath) < 0.1:
                    clusters.append(clnum)
                    added = True
                    break
            if not added:
                reprclusters[totalclnum + 1] = uniqpath
                reprscores[totalclnum + 1] = up['obj']
                totalclnum += 1
                clusters.append(totalclnum)
        print 'Number of clusters: %d' % len(reprclusters.keys())
        if len(reprclusters.keys()) < 25:
            print 'Representative pathway models:'
            for k in reprclusters.keys():
                print reprclusters[k], reprscores[k]
        self.reprclusters = reprclusters
        return reprclusters, clusters

    
    def make_position_matrices(self, startposition=0):
        pathlength = len(self.np_path_array['strrepr'][0].split('->'))
        position = np.empty((len(self.np_path_array['strrepr']), pathlength), dtype=('str_', 32))

        uniqpaths = []
        for i, pathwaystring in enumerate(self.np_path_array['strrepr']):
            strlist = pathwaystring.split('->')
            strlist = [s.strip() for s in strlist]
            position[i] = strlist
            uniqpaths.append(strlist)

        position = np.transpose(position)
        self.get_enzyme_order(position, startposition)
        self.set_ligand_counts(uniqpaths, startposition)
        self.set_interactions(uniqpaths, startposition)

    def get_enzyme_order(self, position, startposition=0):

        def countit(enzyme, enzymelist):
            enz_count = 0
            for elem in enzymelist:
                if elem == enzyme:
                    enz_count += 1 
            return enz_count

        count = {}
        for e in self.enzymes:
            count[e] = [countit(e, position[2*i+1-startposition]) for i in range((position.shape[0]+startposition)/2)]
        for e in self.enzymes:
            count[e] = [100.*x/len(self.np_path_array) for x in count[e]]
        self.enzymes_df = pd.DataFrame(count, index=range(len(self.true_enzymes))).T

    def set_interactions(self, pathlist, startposition=0):
        # assumes path starts with ligand
        #ligs = self.uniqpaths[0][::2]
        all_ligs = []
        for u in pathlist:
            all_ligs.extend(u[startposition::2])
        all_ligs.extend(self.ligands[:])
        all_ligs = list(set(all_ligs))
        ldict = {}
        for l in all_ligs:
            ldict[l] = pd.Series(0, index=self.enzymes)
        df = pd.DataFrame(ldict)
        df_sub = pd.DataFrame(ldict)
        df_prod = pd.DataFrame(ldict)

        n = len(pathlist)
        for p in pathlist:
            for i in range(startposition, (len(p))/2):
                idx = 2*i-startposition
                if idx >= 0:
                    df[p[idx]][p[idx+1]] += 1
                    df_sub[p[idx]][p[idx+1]] += 1
                # for product-enz interactions
                if (idx+2)< len(p):
                    df[p[idx+2]][p[idx+1]] += 1
                    df_prod[p[idx+2]][p[idx+1]] += 1
        # normalize
        self.interaction_df = df.apply(lambda x: x/float(n))
        self.substrate_interaction_df = df_sub.apply(lambda x: x/float(n))
        self.product_interaction_df = df_prod.apply(lambda x: x/float(n))

    def set_ligand_counts(self, pathlist, startposition=0):
        lposition = range(len(self.ligands))
        all_ligs = []
        for u in pathlist:
            all_ligs.extend(u[startposition::2])
        all_ligs.extend(self.ligands[:])
        all_ligs = list(set(all_ligs))

        ldict = {}
        for l in all_ligs:
            ldict[l] = pd.Series(0, index=lposition)

        df = pd.DataFrame(ldict)
        for p in pathlist:
            for i in lposition:
                df[p[2*i + startposition]][i] += 1
        n = len(self.np_path_array)
        self.ligands_df = df.apply(lambda x: x/float(n))

    def enzyme_pos_html(self):
        steps = range(len(self.enzymes_df.columns))
        stepdict = {}
        for s in steps:
            stepdict[s] = s + 1
        tempdf = self.enzymes_df.T.reindex_axis(self.enzymes, axis=1)
        tempdf.rename(index=stepdict, inplace=True)
        string = '<strong>Columns: Enzyme, Rows: Position</strong>'
        return string + tempdf.to_html(float_format='{0:.1f}'.format)

    def ligand_rankings_html(self, startposition=0):
        dockranks = self.get_substrate_rank_docking(startposition)
        dockprodranks = self.get_product_rank_docking(startposition)
        intranks = self.get_ligand_rank_interaction(startposition)
        subranks = self.get_substrate_rank_interaction(startposition)
        prodranks = self.get_product_rank_interaction(startposition)
        htmlstring = '''<table><tr>
                        <th>Enzyme</th>
                        <th>Substrate rank from individual docking run</th>
                        <th>Product rank from individual docking run</th>
                        <th>Substrate rank by integrated approach</th>
                        <th>Product rank by integrated approach</th>
                        <th>Interaction rank by integrated approach</th>
                        </tr>'''
        for i in range(len(self.enzymes)-startposition):
            htmlstring += '''<tr><td>%s</td>
                                 <td>%d</td>
                                 <td>%d</td>
                                 <td>%d</td>
                                 <td>%d</td>
                                 <td>%d</td>
                             </tr>''' % (self.true_enzymes[i+startposition],
                                         dockranks[i],
                                         dockprodranks[i],
                                         subranks[i],
                                         prodranks[i],
                                         intranks[i])
        htmlstring += '</table>'
        return htmlstring

    def get_ligand_ranks(self):
        if self.ligands_df is None:
            print 'Make position matrices first'
        else:
            for pos in range(len(self.ligands)):
                freq = self.ligands_df.get_value(pos, self.ligands[pos])
                rank = (self.ligands_df.ix[pos] > freq).sum() + 1
                print (pos + 1), self.ligands[pos], rank

    def get_substrate_rank_docking(self, startposition=0):
        ranks = []
        dock = pd.read_hdf(self.datafile, 'dock')
        smiles = pd.read_hdf(self.datafile, 'smiles')
        allligs = smiles.index
        es = []
        for i in range(len(self.true_enzymes)-startposition):
            e = self.true_enzymes[i+startposition]

            dockscore = dock.get_value(e, self.ligands[i])
            if dockscore == 0:
                rank = -1
            else:
                allscores = [dock.get_value(e, lig) for lig in allligs]
                allscores = np.array(allscores)
                rank = (dockscore <= allscores).sum()
            ranks.append(rank)
            es.append(e)
        return pd.Series(ranks, index=es)

    def get_product_rank_docking(self, startposition=0):
        ranks = []
        dock = pd.read_hdf(self.datafile, 'dock')
        smiles = pd.read_hdf(self.datafile, 'smiles')
        allligs = smiles.index
        es = []
        for i in range(len(self.true_enzymes)-startposition):
            e = self.true_enzymes[i+startposition]
            if i < len(self.ligands)-1:
                dockscore = dock.get_value(e, self.ligands[i+1])
                if dockscore == 0:
                    rank = -1
                else:
                    allscores = [dock.get_value(e, lig) for lig in allligs]
                    allscores = np.array(allscores)
                    rank = (dockscore <= allscores).sum()
                ranks.append(rank)
                es.append(e)
        return pd.Series(ranks, index=es)


    def get_ligand_rank_interaction(self, startposition=0):
        if self.interaction_df is None:
            print 'Make position matrices first'
        else:
            ranks = []
            es = []
            for i in range(len(self.true_enzymes)-startposition):
                e = self.true_enzymes[i+startposition]
                freq = self.interaction_df.get_value(e, self.ligands[i])
                if freq == 0:
                    rank = -1
                else:
                    rank = (self.interaction_df.ix[e] > freq).sum() + 1
                ranks.append(rank)
                es.append(e)
            return pd.Series(ranks, index=es)

    def get_substrate_rank_interaction(self, startposition=0):
        if self.substrate_interaction_df is None:
            print 'Make position matrices first'
        else:
            ranks = []
            es = []
            for i in range(len(self.true_enzymes)-startposition):
                e = self.true_enzymes[i+startposition]
                freq = self.substrate_interaction_df.get_value(e, self.ligands[i])
                if freq == 0:
                    rank = -1
                else:
                    rank = (self.substrate_interaction_df.ix[e] > freq).sum() + 1
                ranks.append(rank)
                es.append(e)
            return pd.Series(ranks, index=es)

    def get_product_rank_interaction(self, startposition=0):
        if self.product_interaction_df is None:
            print 'Make position matrices first'
        else:
            ranks = []
            es = []
            for i in range(len(self.true_enzymes)-startposition):
                e = self.true_enzymes[i+startposition]
                if i < len(self.ligands)-1:
                    freq = self.product_interaction_df.get_value(e, self.ligands[i+1])
                    if freq == 0:
                        rank = -1
                    else:
                        rank = (self.product_interaction_df.ix[e] > freq).sum() + 1
                    ranks.append(rank)
                    es.append(e)
            return pd.Series(ranks, index=es)

    def get_ranks_by_position(self):
        lranks = []
        eranks = []
        counts = 0.0
        erankdf = self.enzymes_df.rank(axis=0, ascending=False, method='min')
        lrankdf = self.ligands_df.rank(axis=1, ascending=False, method='min')
        for i in range(len(self.true_enzymes)):
            r = erankdf.get_value(self.true_enzymes[i], i)
            eranks.append(r)
            if self.enzymes_df.get_value(self.true_enzymes[i], i) > 0:
                counts += 1
        for i in range(len(self.ligands)):
            r = lrankdf.get_value(i, self.ligands[i])
            lranks.append(r)
            if self.ligands_df.get_value(i, self.ligands[i]) > 0:
                counts += 1
        components = [k for k in roundrobin(self.ligands, self.enzymes)]
        ranks = [k for k in roundrobin(lranks, eranks)]
        return components, ranks, counts

    def calculate_average_interaction(self):
        edges = []
        for i in range(len(self.true_enzymes)):
            val = self.substrate_interaction_df.get_value(self.true_enzymes[i], self.ligands[i])
            edges.append(val)
            if i+1 < len(self.ligands):
                val = self.product_interaction_df.get_value(self.true_enzymes[i], self.ligands[i+1])
                edges.append(val)
        return edges

    def calculate_entropies_by_position(self):
        def compute_entropy(column):
            H = 0
            for i in column:
                if i > 0:
                    H -= i*np.log(i)
            return H

        Hls = [compute_entropy(self.ligands_df.ix[i]) for i in self.ligands_df.index]
        Hes = [compute_entropy(self.enzymes_df.T.ix[i]/100) for i in self.enzymes_df.T.index]
        Hs = [k for k in roundrobin(Hls, Hes)]
        return Hs

    def calculate_frequency_by_position(self):
        positdict = {}
        for i in range(len(self.true_enzymes)):
            te = self.true_enzymes[i]
            positdict[i+1] = self.enzymes_df.get_value(te, i)
        posseries  = pd.Series(positdict)
        posseries.index.name = 'position'
        return posseries


    def evaluate_accuracy_precision(self):
        print 'Mean Entropy: %.3f' % np.mean(self.calculate_entropies_by_position())
        print 'Mean Interaction Weight: %.3f' % np.mean(self.calculate_average_interaction())
        c, ranks, counts = self.get_ranks_by_position()
        print 'Mean Ranks over Position: %.3f' % np.mean(ranks)
        print '  Components appearing at position: %.1f' % counts

    def cluster_by_chem_similarity(self, df, h5datafile):
        if h5datafile is None:
            h5datafile = self.datafile
        if h5datafile is not None:
            chemsim = pd.read_hdf(h5datafile, 'chsim')
            pairwise_dists = np.zeros((len(df.columns), len(df.columns)))
            for i, col1 in enumerate(df.columns):
                for j, col2 in enumerate(df.columns):
                    pairwise_dists[i][j] = chemsim.get_value(col1, col2)
        else:
            pairwise_dists = distance.squareform(distance.pdist(df.T))
            print "No data file provided, clustering by interaction profile not chemical similarity"
        clusters = sch.linkage(pairwise_dists, method='complete')
        den = sch.dendrogram(clusters, no_plot=True, orientation='top')
        df = df.ix[:, den['leaves']]
        return df

    def plot_heatmap(self, df, annotate, aspect_ratio, fig, axmatrix):
        cax = axmatrix.matshow(np.log10(df), origin='upper', cmap=cm.Blues)

        if annotate:
            cp = df.copy()
            for lig, col in cp.iteritems():
                cp[lig] = 0
            for i, l in enumerate(self.ligands):
                if i< len(self.enzymes):
                    cp.set_value(self.enzymes[i], l, 1.0)
                if i > 0:
                    cp.set_value(self.enzymes[i-1], l, 1.0)
            counter = 0
            patches = []
            for i, l in cp.iteritems():
                for j, y in enumerate(l):
                    if y > 0:
                        rect = mpatches.Rectangle([counter-0.5, j-0.5], 1, 1)
                        patches.append(rect)
                counter += 1
            collection = PatchCollection(patches, facecolors='none', edgecolors='orange', lw=2, linestyle='-', hatch='')
            supp = axmatrix.add_collection(collection)
        axmatrix.set_aspect(aspect_ratio)
        axmatrix.set_yticklabels(df.index)
        axmatrix.set_yticks(range(len(df.index)))
        axmatrix.set_xticks([])
        axmatrix.set_xlabel('Ligands')
        axmatrix.set_ylabel('Enzymes')

    def plot_heatmap_ligands(self, percent_cutoff, h5datafile=None, annotate=True, aspect_ratio=4, figurename=None, figsize=(16,8)):
        if self.interaction_df is None:
            print 'Make position matrices first'
        else:
            df = self.interaction_df.copy()

            # remove ligands below cutoff
            num_ligs_above_cutoff = 0
            for lig, col in df.iteritems():
                if np.max(col) < percent_cutoff/100. and lig not in self.ligands:
                    del df[lig]
                else:
                    num_ligs_above_cutoff += 1

            df = self.cluster_by_chem_similarity(df, h5datafile)

            fig = plt.figure(figsize=figsize)
            axmatrix = plt.subplot(1, 1, 1)
            self.plot_heatmap(df, annotate, aspect_ratio, fig, axmatrix)
            if figurename is not None:
                fig.savefig(figurename, dpi=300)

    def plot_heatmap_docking(self, percent_cutoff, h5datafile=None, annotate=True, aspect_ratio=4, figurename=None, figsize=(16,8)):
        dock = pd.read_hdf(self.datafile, 'dock')
        chemsim = pd.read_hdf(self.datafile, 'chsim')
        df = dock.copy()

        for lig, col in df.iteritems():
            if np.max(col) < (1 - percent_cutoff/100) and lig not in self.ligands:
                del df[lig]

        df = self.cluster_by_chem_similarity(df, h5datafile)
        tempdf = df.T.reindex_axis(self.enzymes, axis=1)
        df = tempdf.T

        fig = plt.figure(figsize=figsize)
        axmatrix = plt.subplot(1, 1, 1)
        self.plot_heatmap(df, annotate, aspect_ratio, fig, axmatrix)
        if figurename is not None:
            fig.savefig(figurename, dpi=300)

    def plot_heatmaps(self, percent_cutoff, h5datafile=None, annotate=True, aspect_ratio=4, figurename=None, figsize=(16,16)):
        intdf = self.interaction_df.copy()
        dockdf = pd.read_hdf(self.datafile, 'dock').copy()

        for lig, col in intdf.iteritems():
            if np.max(col) < percent_cutoff/100. and lig not in self.ligands:
                del intdf[lig]
                del dockdf[lig]
        for lig, col in dockdf.iteritems():
            if lig not in intdf.columns:
                del dockdf[lig]

        intdf = self.cluster_by_chem_similarity(intdf, h5datafile)
        dockdf = self.cluster_by_chem_similarity(dockdf, h5datafile)
        tempdf = dockdf.T.reindex_axis(self.enzymes, axis=1)
        dockdf = tempdf.T

        fig = plt.figure(figsize=figsize)
        axmatrix = plt.subplot(2, 1, 1)
        self.plot_heatmap(intdf, annotate, aspect_ratio, fig, axmatrix)
        axmatrix2 = plt.subplot(2, 1, 2)
        self.plot_heatmap(dockdf, annotate, aspect_ratio, fig, axmatrix2)
        if figurename is not None:
            fig.savefig(figurename, dpi=300)

    def get_html_table(self, percent_cutoff, enz_percent_cutoff=1.0, 
                       urlstring_to_image='http://www.genome.jp/Fig/compound/%s.gif'):
        if self.ligands_df is None or self.enzymes_df is None:
            print 'Make position matrices'
        else:
            smilesdf = pd.read_hdf(self.datafile, 'smiles')
            lhtml = []
            ehtml = []
            positionarray = range(len(self.ligands))

            htmlstring = '<table class="results">'
            htmlstring += '<tr><th class="correct" colspan="%d">Correct Pathway</th></tr><tr>' % (len(positionarray)*2+1)
            usedligands = [self.ligands[position] for position in positionarray]
            trackligands = []
            erow = {}
            lrow = {}
            maxcol = 0
            for position in positionarray:
                ligandid = self.ligands[position]
                htmlstring += '''<td class="correct"><div class="container"><div id="%s" class="ligand"></div>
                                 %s</div>
                                 </td>''' % (ligandid, ligid_to_img(ligandid, urlstring_to_image))
                                 #</td>''' % (ligandid, ligid_to_svg(ligandid, smilesdf, size=100))
                lrow[position] = []
                frame = self.ligands_df.T
                frame = frame.sort(columns=position, ascending=False)
                ligands_above_co = frame[position][frame[position] > percent_cutoff]
                for j, freq in ligands_above_co.iteritems():
                    lrow[position].append(j)

                maxcol = max(len(lrow[position]), maxcol)

                if position < len(self.true_enzymes):
                    htmlstring += '<td class="correct" id="E%s"><div class="enzyme">%s</div></td>' % (self.enzymes[position], self.enzymes[position])
                    erow[position] = []
                    poslist = []
                    for i in self.enzymes:
                        poslist.append((i, self.enzymes_df.get_value(i, position)))
                    #for i in self.count.keys():
                    #    poslist.append((i, self.count[i][position]))
                    posar = np.array(poslist, dtype=([('enzyme', np.str_, 10), ('occurrence', np.float32)]))
                    posar = np.sort(posar, order='occurrence')
                    for x in range(len(posar)):
                        enz = posar['enzyme'][len(posar)-x-1]
                        freq = posar['occurrence'][len(posar)-x-1]
                        if freq > enz_percent_cutoff:
                            erow[position].append([enz, freq])
                    maxcol = max(len(erow[position]), maxcol)
            htmlstring += '</tr>'
            htmlstring += '<th colspan="%d">Predicted Pathways</th>' % (len(positionarray)*2+1)
            chemsim = pd.read_hdf(self.datafile, 'chsim')
            for i in range(maxcol):
                htmlstring += '<tr>'
                for position in positionarray:

                    ligs = lrow[position]
                    if i < len(ligs):
                        htmlstring += '<td>'
                        if self.ligands[position] == ligs[i]:
                            htmlstring += '<div class="container match"><div class="ligand" id="%s"></div>' % ligs[i]
                        else:
                            htmlstring += '<div class="container"><div id="%s" class="ligand"></div>' % ligs[i]
                        #tc = chemsim.get_value(self.ligands[position], ligs[i])
                        #htmlstring += '''%s</div>%s</td>''' % (ligid_to_img(ligs[i], urlstring_to_image), ligs[i])
                        ligimgid = ligs[i].split('_')[0]
                        htmlstring += '''%s</div></td>''' % (ligid_to_img(ligimgid, urlstring_to_image))
                        #htmlstring += '''%s</div>%s, %.2f</td>''' % (ligid_to_img(ligs[i], urlstring_to_image), ligs[i], tc)
                        #htmlstring += '''%s</div></td>''' % (ligid_to_svg(ligs[i], smilesdf, size=150))
                    else:
                        htmlstring += '<td></td>'

                    if position < len(self.true_enzymes):
                        enzs = erow[position]
                        if i < len(enzs):
                            enz = enzs[i][0]
                            freq = enzs[i][1]
                            htmlstring += '<td>'
                            if enz == self.enzymes[position]:
                                htmlstring += '<div class="enzyme match" id="E%s">' % enz
                            else:
                                htmlstring += '<div class="enzyme">'
                            htmlstring += '%s</div></td>' % enz
                        else:
                            htmlstring += '<td></td>'
                htmlstring += '</tr>'
            htmlstring += '''</table>'''
            htmlstring += '''
                        <style type="text/css">
                        .match {
                            font-size: 1.5em;
                            //font-weight: bold
                        }
                        table.results, .results td, .results th {
                            border: 1px solid gray
                        }
                        .results th {
                            padding: 5px;
                            font-size: 1.5em
                        }
                        table.results {
                            border-collapse: collapse;
                            font-family: "Helvetica"
                        }
                        table.results >> td {
                            height: 100px;
                            width: auto;
                            padding: 0;
                            margin: 0;
                            text-align: center;
                            vertical-align: middle;
                        }
                        div.enzyme {
                            text-align: center;
                            font-size: 1.5em;
                            padding-left: .8em;
                            padding-right: .8em;
                            //border: 2px hidden black;
                        }
                        div.ligand {
                            opacity: 0.4;
                            top:0;
                            width: 100%;
                            height: 100%;
                            position: absolute;
                            border-radius: 15px;
                            z-index: 0
                        }
                        div.container {
                            //text-align: center;
                            height: 100%;
                            width: auto;
                            position:relative
                        }
                        .results img {
                            height: 100px;
                            width: auto;
                            border-radius: 15px
                            
                        }
                        </style>'''

            htmlstring += '<style type="text/css">'
            colors = get_N_HexCol(len(usedligands))
            for i, lig in enumerate(usedligands):
                htmlstring += '''
                #%s {
                    background-color: #%s
                }
                ''' % (lig, colors[i])
            htmlstring += '</style>'
            return htmlstring

    def paths_to_gml(self, percent_cutoff=0.01, graphfilename='mygraph.gml', startposition=0, datafile=None):

        graph = nx.DiGraph()
        proteins = []
        ligands = []
        edges = []

        best = pathstring_to_graph(self.best_path, startposition)
        true = lists_to_graph(self.true_enzymes, self.ligands)
        dock = pd.read_hdf(self.datafile, 'dock')

        # Substrate-Enzyme pairs
        boolinteractions = self.substrate_interaction_df >= percent_cutoff
        for column, series in self.substrate_interaction_df.iteritems():
            for index, row in series.iteritems():
                if boolinteractions[column][index]:
                    protein = index
                    ligand = column
                    weight = self.substrate_interaction_df[column][index]
                    print ligand, protein
                    if ligand not in graph.nodes():
                        add_ligand_node(graph, ligand, datafile=self.datafile)
                        #string = 'http://www.genome.jp/Fig/compound/%s.gif' % ligand
                        #zinc = 'http://zinc.docking.org/img/sub/%s.gif ' % ligand.strip('ZINC')
                        #graph.add_node(ligand, {'group':1, 'kegg':string, 'zinc':zinc})
                    if protein not in graph.nodes():
                        graph.add_node(protein, {'group':0})
                    present = ''
                    if ligand in true.nodes():
                        if protein in true.successors(ligand):
                            present += 't'
                    if ligand in best.nodes():
                        if protein in best.successors(ligand):
                            present += 'b'
                    dockscore = dock.get_value(protein, ligand)
                    rank = (dock.ix[protein] > dockscore).sum() + 1
                    graph.add_edge(ligand, protein, {'weight':weight, 'present':present, 'dockrank':rank})
        #Product-Enzyme pairs
        boolinteractions = self.product_interaction_df >= percent_cutoff
        for column, series in self.product_interaction_df.iteritems():
            for index, row in series.iteritems():
                if boolinteractions[column][index]:
                    protein = index
                    ligand = column
                    weight = self.product_interaction_df[column][index]
                    print protein, ligand

                    present = ''
                    if ligand not in graph.nodes():
                        #string = 'http://www.genome.jp/Fig/compound/%s.gif' % ligand
                        #zinc = 'http://zinc.docking.org/img/sub/%s.gif ' % ligand.strip('ZINC')
                        #graph.add_node(ligand, {'group':1, 'kegg':string, 'zinc':zinc})
                        add_ligand_node(graph, ligand, datafile=self.datafile)
                    if protein not in graph.nodes():
                        graph.add_node(protein, {'group':0})
                    present = ''
                    if protein in true.nodes():
                        if ligand in true.successors(protein):
                            present += 't'
                    if protein in best.nodes():
                        if ligand in best.successors(protein):
                            present += 'b'
                    dockscore = dock.get_value(protein, ligand)
                    rank = (dock.ix[protein] > dockscore).sum() + 1
                    graph.add_edge(protein, ligand, {'weight':weight, 'present':present, 'dockrank':rank})
        groups = nx.get_node_attributes(true, 'group')
        for edge in true.edges_iter():
            if edge not in graph.edges():
                n1 = edge[0]
                n2 = edge[1]
                if n1 in self.enzymes:
                    protein = n1
                    ligand = n2
                else:
                    protein = n2
                    ligand = n1
                for n in [n1, n2]:
                    if n not in graph.nodes():
                        if n not in self.enzymes:
                            string = 'http://www.genome.jp/Fig/compound/%s.gif' % n
                            zinc = 'http://zinc.docking.org/img/sub/%s.gif ' % n.strip('ZINC')
                            graph.add_node(n, {'group':1, 'name':n, 'kegg':string, 'zinc':zinc})
                        else:
                            graph.add_node(n, {'group':0, 'name':n})
                dockscore = dock.get_value(protein, ligand)
                rank = (dock.ix[protein] > dockscore).sum() + 1
                graph.add_edge(n1, n2, {'present':'f', 'dockrank':rank})
        groups = nx.get_node_attributes(best, 'group')
        for edge in best.edges_iter():
            if edge not in graph.edges():
                n1 = edge[0]
                n2 = edge[1]
                if n1 in self.enzymes:
                    protein = n1
                    ligand = n2
                else:
                    protein = n2
                    ligand = n1
                for n in [n1, n2]:
                    if n not in graph.nodes():
                        if n not in self.enzymes:
                            string = 'http://www.genome.jp/Fig/compound/%s.gif' % n
                            zinc = 'http://zinc.docking.org/img/sub/%s.gif ' % n.strip('ZINC')
                            graph.add_node(n, {'group':1, 'name':n, 'kegg':string, 'zinc':zinc})
                        else:
                            graph.add_node(n, {'group':0, 'name':n})
                dockscore = dock.get_value(protein, ligand)
                rank = (dock.ix[protein] > dockscore).sum() + 1
                graph.add_edge(n1, n2, {'present':'b', 'dockrank':rank})

        nx.write_gml(graph, graphfilename)
        return graph

    def plot_convergence_by_runs(self, run_increment=10, figurename=None):
        maxruns = len(self.conv_test_means)

        fig = plt.figure(figsize=(5, 5))
        ax = plt.subplot(1, 1, 1)
        ax.errorbar(1+run_increment*np.arange(len(self.conv_test_means[0:maxruns:run_increment])),
                    np.array(self.conv_test_means[0:maxruns:run_increment]), fmt='o-',
                    yerr=np.array(self.conv_test_stds[0:maxruns:run_increment]))
        ax.axhline(max(self.conv_test_means), ls='--', color='#333333', lw=2)
        ax.set_ylabel('Clusters observed', fontsize=16)
        ax.set_xlabel('Number of runs', fontsize=16)
        ax.set_ylim(bottom=0, top=1.1*np.max(self.conv_test_means))
        ax.xaxis.labelpad = 5
        ax.yaxis.labelpad = 5

        if figurename is not None:
            fig.savefig(figurename, dpi=300)
        else:
            plt.show()

#########################################################################################
#
# Clustering and other related functions
#

def count_paths(paths, flatclusters, alluniqpaths):
    num_clusters = max(flatclusters)
    clusterdict = {}
    for k in range(1, num_clusters+1):
        clusterdict[k] = 0
    pathclusters = [flatclusters[alluniqpaths.index(path)] for path in paths]
    for pathcluster in pathclusters:
        clusterdict[pathcluster] = 1
        if sum(x > 0 for x in clusterdict.values())==num_clusters:
            break
    return clusterdict

def compare_lists(x, y):
    ref = x
    if len(x) > len(y):
        ref = y
    z = [x[i] == y[i] for i in range(len(ref))]
    return 1-float(sum(z))/max(len(x), len(y))

def cluster_pathways(cutoff, outpickle, fileprefix, totalruns=0, thresholdcl=0.2, runvector=[]):
    starttime = time.time()
    pE = pathens()
    pE.set_fileprefix(fileprefix)
    if totalruns > 0:
        allruns = range(1, totalruns+1)
    elif len(runvector) > 0:
        allruns = runvector[:]
    else:
        print 'Error: pass either list of run numbers or total runs'
        sys.exit()
    pE.get_uniq_paths(allruns, cutoff)
    print 'Number of Paths', len(pE.uniqpathstrs)
    print 'Best %.2f' % (pE.best_score)
    
    #pE.make_clusters(threshold=thresholdcl)
    with open(outpickle, 'w') as handle:
        pickle.dump(pE, handle)
    #print 'Number of Clusters', len(pE.reprclusters.keys())
    #print 'Time:', time.time() - starttime

def update_clusters():
    starttime = time.time()
    picklefile = 'gulonate_out.pickle'
    with open(picklefile, 'r') as handle:
        pE = pickle.load(handle)
    pE.simple_cluster(threshold=0.2)
    outpickle = 'gulonate_out_new_clusters.pickle'
    with open(outpickle, 'w') as outhandle:
        pickle.dump(pE, outhandle)
    print 'Time:', time.time() - starttime

def convergence_test(cutoff, outpickle, fileprefix, totalruns, maxgroupsize, numiters=50, thresholdcl=0.2, outfile='clusters.txt'):
    pE = pathens()
    pE.set_fileprefix(fileprefix)
    pE.get_uniq_paths(range(1, totalruns+1), cutoff)
    pE.make_clusters(threshold=thresholdcl)
    
    print 'Number of Paths', len(pE.uniqpathstrs)
    print 'Number of Clusters', max(pE.flatclusters)
    counts = []
    allruns = range(1, totalruns+1)
    for i in range(numiters):
        print i
        runset = [allruns[i]]
        tempruns = allruns[:]
        tempruns.remove(allruns[i])
        runset = runset + rand.sample(tempruns, maxgroupsize-1)
        countdicts = pE.count_clusters(runset, cutoff)
        clustercounts = [sum(x > 0 for x in k.values()) for k in countdicts]
        counts.append(clustercounts)

    npcounts = np.array(counts)
    print npcounts.mean(axis = 0)
    print npcounts.std(axis = 0)
    print ''
    means = npcounts.mean(axis = 0)
    stds = npcounts.std(axis = 0)

    with open(outfile, 'w') as handle:
        handle.write('[')
        for n in range(len(means) - 1):
            handle.write('%.3f, ' % means[n])
        handle.write('%.3f]' % means[-1])
        handle.write('\n')
        handle.write('[')
        for n in range(len(stds) - 1):
            handle.write('%.3f, ' % stds[n])
        handle.write('%.3f]' % stds[-1])

    pE.conv_test_means = means
    pE.conv_test_stds = stds
        
    with open(outpickle, 'w') as handle:
        pickle.dump(pE, handle)

def print_pathways_in_cluster(picklefile, clusternum, threshold):
    with open(picklefile, 'r') as handle:
        pE = pickle.load(handle)
        pE.make_clusters(threshold=threshold)
        print len(pE.uniqpathstrs)
        print max(pE.flatclusters)
        for i in range(len(pE.flatclusters)):
            if pE.flatclusters[i] == clusternum:
                print pE.uniqpathstrs[i]
                
def get_std_from_random(picklefile):
    sl = pickle.load(open(picklefile, 'r'))
    print 'Mean: %.3f' % np.mean(sl)
    print 'Stdev: %.3f' % np.std(sl)
    return np.std(sl)


####################################################################################
#
# Following functions for ipython notebook interactive viewing
#

# Arguments are the ligand id and the url for the image
def ligid_to_img(ligid, url):
    if ligid.startswith('ZINC'):
        ligidclean = ligid.strip('ZINC')
        ligidclean = ligidclean.lstrip('0')
        fullurl = url % ligidclean
    elif '_' in ligid:
        ligidclean = ligid.split('_')
        ligidclean = ligidclean[0].strip()
        fullurl = url % ligidclean
    else:
        fullurl = url % ligid
    return "<img src='%s' alt='%s'/>" % (fullurl, ligid)

# get hexcodes for N different colors
def get_N_HexCol(N=5):
    HSV_tuples = [(x*.95/N, 1, .9) for x in range(N)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    RGB_tuples = map(lambda x: tuple(map(lambda y: int(y * 255),x)),RGB_tuples)
    HEX_tuples = map(lambda x: tuple(map(lambda y: chr(y).encode('hex'),x)), RGB_tuples)
    HEX_tuples = map(lambda x: "".join(x), HEX_tuples)
    return HEX_tuples

def pathstring_to_table(pathstring, datafile=None, title='', url='http://zinc.docking.org/img/sub/%s.gif'):
    if datafile:
        reactpanel = pd.read_hdf(datafile, 'rsim')
        dockdf = pd.read_hdf(datafile, 'dock')
        #tfdf = pd.read_hdf(datafile, 'tfluor')

    path = pathstring.split(' -> ')
    string = '<tr>'

    if len(title) > 0:
        string+='<td style="font-size:1.5em"><strong>%s</strong></td>' % title
    #val = tfdf.get_value(path[0], path[1])
    for i in range(0, len(path)):
        p = path[i]
        if len(p) > 3:
            string += '<td>%s' % (ligid_to_img(p, url))
            if datafile:
                if i + 1 < len(path):
                    val = dockdf.get_value(path[i+1], p)
                else:
                    val = 0
                string += '<br/><br/>%.3f' % val
            string += '</td>'
        else:
            string += '<td style="padding: 15px"><strong>%s</strong>' % p
            if datafile:
                val = reactpanel.get_value(p, path[i-1], path[i+1])
                string += '<br><br>%.3f' % val
            string += '</td>'
    string += '</tr>'
    return string



####################################################################################
#
# Following functions for making networkx graphs
#

def add_ligand_node(g, l, datafile=None, position=0):
    keggstring = 'http://www.genome.jp/Fig/compound/%s.gif' % l
    ligidclean = l.strip('ZINC')
    ligidclean = ligidclean.lstrip('0')
    zinc = 'http://zinc.docking.org/img/sub/%s.gif' % ligidclean
    attributes_dict = {}
    attributes_dict['group'] = 1
    attributes_dict['kegg'] = keggstring
    attributes_dict['zinc'] = zinc
    attributes_dict['position'] = position
    if datafile is not None:
        smilesdf = pd.read_hdf(datafile, 'smiles')
        smiles = smilesdf.get_value(l)
        attributes_dict['smiles'] = smiles
    g.add_node(l, attributes_dict)

def add_protein_node(g, p):
    g.add_node(p, {'group':0})

#startposition indicates where first ligand is
def pathstring_to_graph(pathstring, startposition=0):

    strlist = pathstring.split('->')
    strlist = [s.strip() for s in strlist]
    graph = nx.DiGraph()

    ligandnum = 0

    if startposition > 0:
        prot = strlist[0]
        lig = strlist[0]
        add_ligand_node(graph, lig, position=ligandnum)
        add_protein_node(graph, prot)
        ligandnum += 1

    for i in range(startposition, len(strlist)/2):
        lig = strlist[2*i - startposition]
        prot = strlist[2*i+1 - startposition]
        if lig not in graph.nodes():
            add_ligand_node(graph, lig, position=ligandnum)
            ligandnum += 1
        if prot not in graph.nodes():
            add_protein_node(graph, prot)
        graph.add_edge(lig, prot)

        if (2*i+2-startposition) < len(strlist):
            lig2 = strlist[2*i+2-startposition]
            if lig2 not in graph.nodes():
                add_ligand_node(graph, lig2, position=ligandnum)
                ligandnum += 1
            graph.add_edge(prot, lig2)
    return graph

def lists_to_graph(protlist, liglist):
    graph = nx.DiGraph()
    for p in protlist:
        add_protein_node(graph, p)
    for lig in liglist:
        add_ligand_node(graph, lig)

    for i in range(len(liglist)-1):
        graph.add_edge(liglist[i], protlist[i])
        graph.add_edge(protlist[i], liglist[i+1])
    if len(protlist) == len(liglist):
        graph.add_edge(liglist[-1], protlist[-1])
    return graph

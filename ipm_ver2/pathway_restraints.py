#!/bin/python
#from rxngraphs import rxnGraph

from tables import open_file
import numpy as np
import random as rand
import pandas as pd
from itertools import combinations

class restraint(object):
    def __init__(self, name="restraint"):
        self.name = name
        self.score = None

    def __str__(self):
        string = '%s score: %.3f' % (self.name, self.score)
        return string

    def compute_score(self, rg):
        self.score = 0.0
        return self.score

    def to_json(self, rg):
        return []

class reaction_restraint(restraint):
    def __init__(self, name='chemical transformation'):
        self.name = name
        self.score = None
    
    def compute_score(self, rg):
        rxnscore = 0
        count = 0
        for p in rg.protein_allnodes:
            substrate = rg.graph.predecessors(p)
            product = rg.graph.successors(p)
            for i in substrate:
                for j in product:
                    proteinlabel = rg.graph.node[p]['label']
                    sublabel = rg.graph.node[i]['label']
                    prodlabel = rg.graph.node[j]['label']

                    if proteinlabel in rg.data.rsimpn.items:
                        newterm = rg.data.rsimpn.get_value(proteinlabel, sublabel, prodlabel)
                        rxnscore += newterm
                        count += 1
        if count == 0:
            self.score = 0
        else:
            rxnscore = rxnscore/count
            self.score = rxnscore
        return self.score

    def to_json(self, rg):
        jsrestraints = []
        for p in rg.protein_allnodes:
            rxnrestraint = {}
            substrate = rg.graph.predecessors(p)
            product = rg.graph.successors(p)
            for i in substrate:
                for j in product:
                    proteinlabel = rg.graph.node[p]['label']
                    sublabel = rg.graph.node[i]['label']
                    prodlabel = rg.graph.node[j]['label']
                    score = rg.data.rsimpn.get_value(proteinlabel, sublabel, prodlabel)
                    rxnrestraint['restraint type'] = 'chemical transformation'
                    rxnrestraint['edges'] = [{'source':sublabel, 'target':proteinlabel},
                                             {'source':proteinlabel, 'target':prodlabel},
                                             {'source':sublabel, 'target':prodlabel}]
                    rxnrestraint['score'] = score
                    rxnrestraint['directed'] = True
            jsrestraints.append(rxnrestraint)
        return jsrestraints

class thiolase_rxn_restraint(reaction_restraint):
    def compute_score(self, rg):
        rxnscore = 0
        count = 0
        for p in rg.protein_allnodes:
            substrate = rg.graph.predecessors(p)
            product = rg.graph.successors(p)
            for i in substrate:
                for j in product:
                    proteinlabel = rg.graph.node[p]['label']
                    sublabel = rg.graph.node[i]['label']
                    prodlabel = rg.graph.node[j]['label']

                    if proteinlabel in rg.data.rsimpn.items:
                        newterm = rg.data.rsimpn.get_value(proteinlabel, prodlabel, sublabel)
                        rxnscore += newterm
                        count += 1
        if count == 0:
            self.score = 0
        else:
            rxnscore = rxnscore/count
            self.score = rxnscore
        return self.score


class chemsim_restraint(restraint):
    def __init__(self, name='chemical similarity'):
        self.name = name
        self.score = None

    def compute_score(self, rg):
        tcscore = 0
        count = 0
        for p in rg.protein_allnodes:
            substrate = rg.graph.predecessors(p)
            product = rg.graph.successors(p)
            for i in substrate:
                for j in product:
                    sublabel = rg.graph.node[i]['label']
                    prodlabel = rg.graph.node[i]['label']
                    tcscore += rg.data.chsimdf.get_value(sublabel, prodlabel)
                    count += 1.
        if count == 0:
            self.score = 0
        else:
            self.score = tcscore/count
        return self.score

class sea_restraint(restraint):
    def __init__(self, name='SEA'):
        self.name = name
        self.score = None

    def compute_score(self, rg):
        seascore = 0
        count = 0
        for p in rg.ligand_allnodes:
            prot1 = rg.graph.predecessors(p)
            prot2 = rg.graph.successors(p)
            for i in prot1:
                for j in prot2:
                    label1 = rg.graph.node[i]['label']
                    label2 = rg.graph.node[j]['label']
                    seascore += rg.data.seadf.get_value(label1, label2)
                    count += 1.
        if count == 0:
            self.score=0
        else:
            seascore = (seascore - rg.data.aggmean)/rg.data.aggstd
            self.score = seascore
        return self.score

    def to_json(self, rg):
        jsrestraints = []

        for p in rg.ligand_allnodes:
            prot1 = rg.graph.predecessors(p)
            prot2 = rg.graph.successors(p)
            for i in prot1:
                for j in prot2:
                    label1 = rg.graph.node[i]['label']
                    label2 = rg.graph.node[j]['label']
                    seascore = rg.data.seadf.get_value(label1, label2)
                    seascore = (seascore - rg.data.aggmean)/rg.data.aggstd
                    searestraint  = {}
                    searestraint['restraint type'] = 'SEA'
                    searestraint['edges'] = [{'source':label1, 'target':label2}]
                    searestraint['score'] = seascore
                    searestraint['directed'] = False
                    jsrestraints.append(searestraint)
        return jsrestraints

class dock_restraint(restraint):
    def __init__(self, name='docking'):
        self.name = name
        self.score = None
        self.include_product = True

    def compute_score(self, rg, test=False):
        dockscore = 0
        count = 0
        for p in rg.protein_allnodes:
            substrate = rg.graph.predecessors(p)
            for i in substrate:
                label1 = rg.graph.node[p]['label']
                label2 = rg.graph.node[i]['label']
                newterm = rg.data.dockdf.get_value(label1, label2)
                dockscore += newterm
                count += 1.
            if self.include_product:
                product = rg.graph.successors(p)
                for i in product:
                    label1 = rg.graph.node[p]['label']
                    label2 = rg.graph.node[i]['label']
                    if test:
                        print 'Label 1: %s, Label 2: %s, %.2f' % (label1, label2, newterm)
                    newterm = rg.data.dockdf.get_value(label1, label2)
                    dockscore += newterm
                    count += 1.

        if rg.transporter is not None and rg.transporter in rg.data.dockdf.index:
            first_substrate = rg.get_start()
            transporter_substrate = rg.graph.node[first_substrate]['label']
            dockscore += rg.data.dockdf.get_value(rg.transporter, transporter_substrate)
            #self.transporter_substrate = self.graph.node[first_substrate]['label']
            #dockscore += self.data.dockdf.get_value(self.transporter, self.transporter_substrate)
            count += 1.

        if count == 0:
            self.score = 0
        else:
            #dockscore = dockscore/count

            # When combining z-scores, should sq root of count
            dockscore = dockscore/(count**0.5)
            self.score = dockscore
        return self.score

    def to_json(self, rg):
        jsrestraints = []
        
        for p in rg.protein_allnodes:
            substrate = rg.graph.predecessors(p)
            for i in substrate:
                label1 = rg.graph.node[p]['label']
                label2 = rg.graph.node[i]['label']
                dockscore = rg.data.dockdf.get_value(label1, label2)
                sub_dockrestraint = {}
                sub_dockrestraint['restraint type'] = 'docking'
                sub_dockrestraint['edges'] = [{'source':label2, 'target':label1}]
                sub_dockrestraint['directed'] = True
                sub_dockrestraint['score'] = dockscore
                jsrestraints.append(sub_dockrestraint)

            product = rg.graph.successors(p)
            for i in product:
                label1 = rg.graph.node[p]['label']
                label2 = rg.graph.node[i]['label']
                dockscore = rg.data.dockdf.get_value(label1, label2)
                prod_dockrestraint = {}
                prod_dockrestraint['restraint type'] = 'docking'
                prod_dockrestraint['edges'] = [{'source':label2, 'target':label1}]
                prod_dockrestraint['directed'] = True
                prod_dockrestraint['score'] = dockscore
                jsrestraints.append(prod_dockrestraint)
        return jsrestraints

class evidence_restraint(restraint):
    def __init__(self, name='homology evidence'):
        self.name = name
        self.score = None

    def compute_score(self, rg):
        prots = rg.data.eviddf.columns
        evidscore = 0.
        count = 0
        for prot in prots:
            pnodelist = (n for n in rg.graph if rg.graph.node[n]['label']==prot)
            pnode = list(pnodelist)[0]
            substrate = rg.graph.predecessors(pnode)
            for s in substrate:
                substratelabel = rg.graph.node[s]['label']
                evidscore += rg.data.eviddf.get_value(substratelabel, prot)
                count += 1
        if count == 0:
            self.score = 0.0
        else:
            self.score = evidscore/count
        return self.score

    def to_json(self, rg):
        jsrestraints = []
        for prot in rg.data.eviddf.columns:
            pnodelist = (n for n in rg.graph if rg.graph.node[n]['label']==prot)
            pnode = list(pnodelist)[0]
            substrate = rg.graph.predecessors(pnode)
            for s in substrate:
                substratelabel = rg.graph.node[s]['label']
                evrestraint = {}
                evrestraint['score'] = rg.data.eviddf.get_value(substratelabel, prot)
                evrestraint['restraint type'] = self.name
                evrestraint['directed'] = True
                evrestraint['edges'] = [{'source': substratelabel,
                                         'target': prot}]
                jsrestraints.append(evrestraint)
        return jsrestraints

class cmetab_restraint(restraint):
    def __init__(self, name='central metab endpoint'):
        self.name = name
        self.score = None

    def compute_score(self, rg):
        length = len(rg.graph.nodes())
        #length = len(rg.protein_members) + len(rg.ligand_members)
        last_metabolite = rg.graph.node[length-1]['label']
        cmetabscore = rg.data.cmetabdf.get_value(last_metabolite)
        self.score = cmetabscore
        return self.score

class tfluor_restraint(restraint):
    def __init__(self, name='thermofluor'):
        self.name = name
        self.score = None

    def compute_score(self, rg):
        tfscore = 0
        count = 0
        first_substrate = rg.get_start()
        transporter_substrate = rg.graph.node[first_substrate]['label']

        if rg.transporter in rg.data.tfluordf.index:
            val = rg.data.tfluordf.get_value(rg.transporter, transporter_substrate)
            tfscore += val
            count += 1
        else:
            tfscore += rg.data.tfluordf.mean().mean()

        self.score = tfscore
        return self.score

class genecluster_restraint(restraint):
    def __init__(self, name='gene cluster'):
        self.name = name
        self.score = None

    def compute_score(self, rg):
        self.score = rg.data.genecluster.get_cluster_score(rg.protein_members)
        return self.score

class upreg_restraint(restraint):
    def __init__(self, name='upregulated genes'):
        self.name = name
        self.score = None

    def compute_score(self, rg):
        self.score = rg.data.upregulated_genes.get_upreg_score(rg.ligand_members, rg.protein_members)
        return self.score

class ko_restraint(restraint):
    def __init__(self, name='knockout phenotype'):
        self.name = name
        self.score = None

    def compute_score(self, rg):
        self.score = rg.data.knockouts.get_ko_score(rg.ligand_members, rg.protein_members)
        return self.score

class dataTables(object):
    def __init__(self, proteins=[], ligands=[]):
        self.dockdf = None # Docking scores Data Frame (L-P)
        self.seadf = None  # SEA scores Data Frame (P-P)
        self.rsimpn = None # Reaction chemical similarity Panel (L-P-L)
        self.chsimdf = None # Chemical similarity Data Frame (L-L)
        self.tfluordf = None # Thermofluor hit, chem sim Data Frame (L-P)
        self.gexpdf = None
        self.mcorrdf = None
        self.cmetab = None
        self.genecluster = None # Gene cluster info
        self.eviddf = None # Evidence of known substrate
        self.knockouts = None

        self.proteins = proteins # list of all protein ids
        self.ligands = ligands   # list of ligand ids

        self.aggmean = None
        self.aggstd = None

    def readInTables(self, h5filename):
        #print 'Reading in data from %s' % h5filename
        self.dockdf = pd.read_hdf(h5filename, 'dock')
        self.seadf = pd.read_hdf(h5filename, 'sea')
        self.setBackgroundFromFile(h5filename)
        self.rsimpn = pd.read_hdf(h5filename, 'rsim')
        self.chsimdf = pd.read_hdf(h5filename, 'chsim')

    def readInSelectedTables(self, h5filename, tablelist):
        #print 'Reading in data from %s' % h5filename
        if 'dock' in tablelist:
            self.dockdf = pd.read_hdf(h5filename, 'dock')
        if 'sea' in tablelist:
            self.seadf = pd.read_hdf(h5filename, 'sea')
            self.setBackgroundFromFile(h5filename)
        if 'rsim' in tablelist:
            self.rsimpn = pd.read_hdf(h5filename, 'rsim')
        if 'chsim' in tablelist:
            self.chsimdf = pd.read_hdf(h5filename, 'chsim')
        if 'tfluor' in tablelist:
            self.tfluordf = pd.read_hdf(h5filename, 'tfluor')
        if 'gexp' in tablelist:
            self.gexpdf = pd.read_hdf(h5filename, 'gexp')
        if 'mcorr' in tablelist:
            self.mcorrdf = pd.read_hdf(h5filename, 'mcorr')
        if 'cmetab' in tablelist:
            self.cmetabdf = pd.read_hdf(h5filename, 'cmetab')
        if 'evid' in tablelist:
            self.eviddf = pd.read_hdf(h5filename, 'evidence')

    def readInGeneExpressionData(self, h5filename):
        self.gexpdf = pd.read_hdf(h5filename, 'gexp')

    def readInThermofluor(self, h5filename):
        self.tfluordf = pd.read_hdf(h5filename, 'tfluor')

    def readInMetabolomicsCorr(self, h5filename):
        self.mcorrdf = pd.read_hdf(h5filename, 'mcorr')

    def readInCentralMetabComparison(self, h5filename):
        self.cmetabdf = pd.read_hdf(h5filename, 'cmetab')

    def readInEvidence(self, h5filename):
        self.eviddf = pd.read_hdf(h5filename, 'evidence')

    # Input: filename - name of file that contains ligand ids (assuming id is in second column)
    #        delimiter - non-whitespace delimiter that separates smiles and id
    #        keepmultiples - False (default) reads in molecules, treating those ids with suffix _ 
    #                        as duplicates. True reads in all molecules as unique ids
    def readInLigands(self, filename, delimiter=None, keepmultiples=False):
        with open(filename, 'r') as handle:
            lines = handle.readlines()
            for line in lines:
                fields = []
                if delimiter is not None:
                    fields = line.strip().split(delimiter)
                else:
                    fields = line.strip().split()
                if keepmultiples:
                    lig = fields[1]
                else:
                    lig = fields[1].split('_')[0]
                if lig not in self.ligands:
                    self.ligands.append(lig)
        print '%d ligands from %s' % (len(self.ligands), filename)

    def setGeneCluster(self, refgenes, allgenes):
        gc = gene_cluster(refgenes)
        gc.get_stats(allgenes)
        self.genecluster = gc

    def setUpregulatedGenes(self, ligand, upreg_genes, refgenes):
        ur = upregulated_genes(ligand, upreg_genes, refgenes)
        ur.get_stats()
        self.upregulated_genes = ur

    def setKnockoutPhenotypes(self, ligand, ko_genes, refgenes):
        ko = knockout_phenotypes(ligand, ko_genes, refgenes)
        ko.get_stats()
        self.knockouts = ko

    def setBackground(self, seamean, seastd):
        self.aggmean = seamean
        self.aggstd = seastd

    def setBackgroundFromFile(self, h5filename):
        if self.aggmean is None:
            with open_file(h5filename, 'r') as seadata_handle:
                attributes = seadata_handle.root.sea.block0_values.attrs
                self.aggmean = attributes.seamean
                self.aggstd = attributes.seastd

class knockout_phenotypes(object):
    def __init__(self, ligand, kogenes, refgenes):
        self.ligand = ligand
        self.kogenes = kogenes
        self.refgenes = refgenes
        self.mean = None
        self.stdev = None

    def get_stats(self):
        scores = []
        for i in range(3, len(self.refgenes)+1):
            combs = combinations(self.refgenes, i)
            for comb in combs:
                score = 0.0
                for gene in self.kogenes:
                    if gene in comb:
                        score += 1.0
                scores.append(score)
        self.mean = np.mean(scores)
        self.stdev = np.std(scores)

    def get_ko_score(self, ligands, pathgenes):
        score = 0.0
        if self.mean is None:
            print 'Distribution not calculated yet'
        else:
            if self.ligand in ligands:
                for gene in pathgenes:
                    if gene in self.kogenes:
                        score += 1.0
            return (score - self.mean)/self.stdev

class upregulated_genes(object):
    def __init__(self, ligand, upreg_genes, refgenes):
        self.ligand = ligand
        self.refgenes = refgenes
        self.upreg_genes = upreg_genes

    def get_stats(self):
        scores = []
        for i in range(3, len(self.refgenes)+1):
            combs = combinations(self.refgenes, i)
            for comb in combs:
                score = 0.0
                for gene in self.upreg_genes:
                    if gene in comb:
                        score += 1.0
                score = score/max(len(self.upreg_genes), i)
                scores.append(score)

        self.mean = np.mean(scores)
        self.stdev = np.std(scores)

    def get_upreg_score(self, ligands, pathgenes):
        score = 0.0
        if self.mean is None:
            print 'Distribution not calculated yet'
        else:
            if self.ligand in ligands:
                for gene in pathgenes:
                    if gene in self.upreg_genes:
                        score += 1.0
            score = score/max(len(self.upreg_genes), len(pathgenes))
            return (score - self.mean)/self.stdev

class gene_cluster(object):
    def __init__(self, refgenes):
        self.refpairs = set(combinations(refgenes, 2))
        self.mean = None
        self.stdev = None

    def get_stats(self, genes):
        scores = []
        for i in range(3, len(genes)+1):
            combs = combinations(genes, i)
            for comb in combs:
                pairs = set(combinations(comb, 2))
                length = len(pairs)
                counts = len(self.refpairs.intersection(pairs))
                rev = [(p[1], p[0]) for p in pairs]
                counts += len(self.refpairs.intersection(rev))
                score = float(counts)/max(length, len(self.refpairs))
                scores.append(score)
        self.mean = np.mean(scores)
        self.stdev = np.std(scores)

    def get_cluster_score(self, pathgenes):
        pairs = set(combinations(pathgenes, 2))
        length = len(pairs)
        reflength = len(self.refpairs)
        counts = len(self.refpairs.intersection(pairs))
        rev = [(p[1], p[0]) for p in pairs]
        counts += len(self.refpairs.intersection(rev))

        score = float(counts)/max(length, reflength)
        if self.mean is None:
            print 'Distribution not calculated yet'
        else:
            return (score - self.mean)/self.stdev


class background_sea_data(object):
    def __init__(self, dictmeans={}, dictstds={}, maxN=20, numsamples=100000):
        self.dictmeans = dictmeans
        self.dictstds = dictstds
        self.maxN = maxN
        self.numsamples = numsamples

    def calculate_background_sets(self, seafile):
        elist, dlist = get_sea_dict_from_file(seafile)
        seavals = pd.DataFrame(dlist, index=elist)

        dictmeans = {}
        dictstds = {}

        allproteins = seavals.index
        for num in range(3, self.maxN):
            scores = []
            for i in range(self.numsamples):
                subset = rand.sample(allproteins, num)
                s = get_score(subset, seavals)
                scores.append(s)
            dictmeans[num] = np.mean(scores)
            dictstds[num] = np.std(scores)
        print "Random sets calculated"
        self.dictmeans = dictmeans
        self.dictstds = dictstds

    def calculate_background_sets_from_df(self, df):
        dictmeans = {}
        dictstds = {}
        allproteins = df.index
        for num in range(3, self.maxN):
            scores = []
            for i in range(self.numsamples):
                subset = rand.sample(allproteins)
                s = get_score(subset, df)
                scores.append(s)
            dictmeans[num] = np.mean(scores)
            dictstds[num] = np.std(scores)
        print "Random sets calculated"
        self.dictmeans = dictmeans
        self.dictstds = dictstds


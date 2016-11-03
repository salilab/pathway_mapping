#!/bin/python

import sys
from optparse import OptionParser
from pathway_restraints import *
from sample_graph import Paths
from tables import open_file
from zscore_compute import set_up_data
from zscore_compute import ligandData
import pandas as pd

kinase=['[C:1][O:2][H]>>[C:1][O:2][P](O)(O)=O', 
        '[C:1][O:2][P](O)(O)=O>>[C:1][O:2][H]']
dehydratase=['[H][O:1][C:2]([H])[C:3]([H])[O][H]>>[O:1]=[C:2]([H])[C:3]([H])[H]',
             '[O:1]=[C:2]([H])[C:3]([H])[H]>>[H][O:1][C:2]([H])[C:3]([H])[O][H]']
aldolase=['[A:1][C:2][C:3]([O:4][H])[A:5]>>[A:1][C:2][H].[O:4]=[C:3][A:5]']
dehydrogenase=['[H][O:1][C:2]([H])([A:3])[A:4]>>[O:1]=[C:2]([A:3])([A:4])',
               '[O:1]=[C:2]([A:3])([A:4])>>[H][O:1][C:2]([H])([A:3])[A:4]']
reductase=['[O:1]=[C:2]([A:3])([A:4])>>[H][O:1][C@:2]([H])([A:3])[A:4]',
           '[H][O:1][C@:2]([H])([A:3])[A:4]>>[O:1]=[C:2]([A:3])([A:4])']

def get_ligands_gulonate():
    from compound_search import enzymeData, write_out_candidates
    starttime=time.time()
    tcmin = 0.75
    reactions = {'1':dehydrogenase,
                 '2':reductase,
                 '3':dehydratase,
                 '4':kinase,
                 '5':aldolase}
    
    reactionlist = [x for sublist in reactions.values() for x in sublist]

    keggsmifile='background_data/ZINC_all_neutral.smi'
    smilesfile='gulonate_test/top1000_smiles_cleaned.smi'

    enzs = [enzymeData(dockingfile=smilesfile, reactionlist=reactions['1']),
            enzymeData(dockingfile=smilesfile, reactionlist=reactions['2']),
            enzymeData(dockingfile=smilesfile, reactionlist=reactions['3']),
            enzymeData(dockingfile=smilesfile, reactionlist=reactions['4']),
            enzymeData(dockingfile=smilesfile, reactionlist=reactions['5'])]

    candidates=[]
    keggdata = ligandData(smilesfile=keggsmifile)
    for e in enzs:
        e.get_molids(column=1, cutoff=10000)
        e.molecules = ['%s' % m for m in e.molecules]
        newc = e.find_similar_compounds(keggdata, tcmin=tcmin)
        candidates.extend(newc)
        print '> Number of candidates: %d' % len(newc)
    print time.time() - starttime
    candidates = list(set(candidates))
    write_out_candidates('gulonate_test/top1000_similar_updated.smi', candidates, keggdata)


def set_up_gulonate(ligandfile, seafile, datafile, evidencefile, thermofluorfile, cmetabfile):
    reactions = {'1':dehydrogenase,
                 '2':reductase,
                 '3':dehydratase,
                 '4':kinase,
                 '5':aldolase}

    dockfiles = {'0':'gulonate_test/old_kegg/combine.scores.C9MHP2',
                 '1':'gulonate_test/old_kegg/combine.scores.C9MHP1',
                 '2':'gulonate_test/old_kegg/combine.scores.C9MHP6',
                 '3':'gulonate_test/old_kegg/combine.scores.C9MHN9',
                 '4':'gulonate_test/old_kegg/combine.scores.C9MHP5',
                 #'5':'gulonate_test/combine.scores.1FQ0'}
                 #'5':'gulonate_test/old_kegg/combine.scores.trim1_C9MHP7'}
                 '5':'gulonate_test/combine.scores.C9MHP7.rescore'}
                 #'5':'gulonate_test/old_kegg/combine.scores.C9MHP7'}

    enzymes=['0', '1', '2', '3', '4', '5']
    set_up_data(ligandfile, None, dockfiles, reactions, seafile, datafile, useMorganFP=True, proteins=enzymes)
    # evidencefile
    gulonate_evidence(ligandfile, datafile, ['3'], [evidencefile])
    #thermofluorfile
    gulonate_thermofluor(ligandfile, datafile, thermofluorfile)
    #cmetabfile
    gulonate_central_metab_comp(ligandfile, datafile, cmetabfile)

# Adds a dummy enzyme 
def set_up_gulonate_dummy(ligandfile, seafile, datafile, evidencefile, thermofluorfile, cmetabfile):
    set_up_gulonate(ligandfile, seafile, datafile, evidencefile, thermofluorfile, cmetabfile)
    pn = pd.read_hdf(datafile, 'rsim')
    df = pd.read_hdf(datafile, 'chsim')
    rsmean = np.mean(df.values)
    rsstd = np.std(df.values)
    zscorefxn = lambda x: (x - rsmean) / rsstd
    df = df.apply(zscorefxn)
    pn['D'] = df
    pn.to_hdf(datafile, 'rsim')
    dockdf = pd.read_hdf(datafile, 'dock')
    s = pd.DataFrame(index=['D'], columns=dockdf.columns)
    s.ix['D'] = 0.0
    dockdf = dockdf.append(s)
    dockdf.to_hdf(datafile, 'dock')

def gulonate_thermofluor(ligandfile, datafile, thermofluor):
    lig = ligandData(smilesfile=ligandfile, useMorganFP=True)
    tf_fingerprints = []
    with open(thermofluor, 'r') as tfhandle:
        lines = tfhandle.readlines()
        for line in lines:
            smiles = line.strip()
            fp = smilesToMorganFingerprint(smiles)
            tf_fingerprints.append(fp)
    tfdata = []
    tfdict = {}
    for fpkey in lig.fingerprints.keys():
        maxtc = 0.0
        for tffp in tf_fingerprints:
            tc = DataStructs.TanimotoSimilarity(tffp, lig.fingerprints[fpkey])
            maxtc = max(tc, maxtc)
        tfdict[fpkey] = maxtc
    df = pd.DataFrame(tfdict, index=['0'])
    meantc = np.mean(df.values)
    stdtc = np.std(df.values)
    zscorefxn = lambda x: (x - meantc) / stdtc
    df = df.apply(zscorefxn)
    df.to_hdf(datafile, 'tfluor')

def gulonate_evidence(ligandfile, datafile, enzymes, filenames):
    lig = ligandData(smilesfile=ligandfile, useMorganFP=True)
    enzdict = {}
    for i in range(len(enzymes)):
        tcdict = lig.compare_metabolites_to_list(filenames[i])
        series = pd.Series(tcdict)
        meantc = np.mean(series.values)
        stdtc = np.std(series.values)
        zscorefxn = lambda x : (x - meantc) / stdtc
        series = series.apply(zscorefxn)
        enzdict[enzymes[i]] = series
    df = pd.DataFrame(enzdict)
    df.to_hdf(datafile, 'evidence')

def gulonate_central_metab_comp(ligandfile, datafile, cmetabfile):
    lig = ligandData(smilesfile=ligandfile, useMorganFP=True)
    tcdict = lig.compare_metabolites_to_list(cmetabfile)
    series = pd.Series(tcdict)
    meantc = np.mean(series.values)
    stdtc = np.std(series.values)
    zscorefxn = lambda x : (x - meantc) / stdtc
    series = series.apply(zscorefxn)
    series.to_hdf(datafile, 'cmetab')

def run_gul(number, steps, proteins, ligandfile, datafile, fileprefix,
                          seed=None, kparameter_set=0, parameter_set=0):
    num_prots = len(proteins)
    num_ligands = num_prots + 1

    outfile = '%s%d.h5' % (fileprefix, int(number))

    d = dataTables(proteins=proteins[:])
    d.readInLigands(ligandfile)
    d.readInTables(datafile)
    d.readInThermofluor(datafile)
    d.readInCentralMetabComparison(datafile)
    d.readInEvidence(datafile)
    
    with open_file(outfile, 'w') as h5outhandle:
        ps = Paths(d, h5outhandle)
        res1 = reaction_restraint('RSIM')
        res2 = sea_restraint('SEA')
        res3 = dock_restraint('VS')
        res4 = evidence_restraint('EV')
        res5 = cmetab_restraint('CM')
        res6 = tfluor_restraint('TF')
        restraints = [res1, res2, res3, res4, res5, res6]
        ratio = ps.sample_labels_monte_carlo(steps, num_prots, num_ligands, restraints, 
                                             kparameter_set=kparameter_set, parameter_set=parameter_set,
                                             transporter='0')

def test_prediction(datafile):
    from rxngraphs import modelGraph

    proteins = ['1', '2', '3', '4', '5']
    ligands = ['ZINC03869787',
               'ZINC04095492',
               'ZINC02040884',
               'ZINC01532568',
               'ZINC01529165',
               'ZINC03869936']

    inters = []
    for i in range(len(proteins)):
        inters.append([ligands[i], proteins[i]])
        if i < len(ligands)-1:
            inters.append([proteins[i], ligands[i+1]])
    d = dataTables(proteins=proteins, ligands=ligands)
    d.readInTables(datafile)
    d.readInThermofluor(datafile)
    d.readInCentralMetabComparison(datafile)
    d.readInEvidence(datafile)
    g = modelGraph(data=d, proteins=proteins, ligands=ligands, interactions=inters)
    g.add_transporter('0')

    res1 = reaction_restraint('RSIM')
    res2 = sea_restraint('SEA')
    res3 = dock_restraint('VS')
    res4 = evidence_restraint('EV')
    res5 = cmetab_restraint('CM')
    res6 = tfluor_restraint('TF')

    restraints = [res1, res2, res3, res4, res5, res6]
    g.add_restraints(restraints)
    g.compute_scores()
    print g
    g.print_scores()

def gulonate_sample_random(steps, prots, ligandfile, datafile, randomh5, randompicklefile):
    print 'sampling %d random graphs' % steps
    with open_file(randomh5, mode='w', title='random') as h5outhandle:

        d = dataTables(proteins=prots)
        d.readInLigands(ligandfile)
        d.readInTables(datafile)
        d.readInThermofluor(datafile)
        d.readInCentralMetabComparison(datafile)
        d.readInEvidence(datafile)

        ps = Paths(d, h5outhandle)
        res1 = reaction_restraint('RSIM')
        res2 = sea_restraint('SEA')
        res3 = dock_restraint('VS')
        res4 = evidence_restraint('EV')
        res5 = cmetab_restraint('CM')
        res6 = tfluor_restraint('TF')
        restraints = [res1, res2, res3, res4, res5, res6]
 
        ps.random_graphs(steps, len(prots), len(prots)+1, restraints, transporter='0')

    scorelist = []
    print 'getting stats'
    with open_file(randomh5, 'r') as h5file:
        table = h5file.root.paths.pathTable
        for row in table.iterrows():
            scorelist.append(row['obj'])
    sl = np.array(scorelist)
    pickle.dump(sl, open(randompicklefile, 'w'))
    print 'Mean: %.3f' %  np.mean(sl)
    print 'Stdev: %.3f' % np.std(sl)


def analyze(fileprefix, picklefile, randompicklefile, clusteroutfile):
    from pathway_analysis import cluster_pathways
    from pathway_analysis import pathens
    from pathway_analysis import convergence_test

    #outfile = '/scrapp2/scalhoun/output/gulonate_test_' 
    #picklefile = '/netapp/home/scalhoun/bin/repo/output/gul_2sd.pickle'
    #sl = pickle.load(open('gul_random_scores.pickle', 'r'))

    numstd = 2
    totalruns = 1000

    best = None
    for i in range(1, totalruns+1):
        currfile = '%s%d.h5' % (fileprefix, i)
        with open_file(currfile, 'r') as h5:
            currtable = h5.root.paths.pathTable
            runbest = currtable.attrs.bestscore
            best = np.max([runbest, best])

    cutoff = best - 2*np.std(sl)
    sl = pickle.load(open(randompicklefile, 'r'))
    #cluster_pathways(cutoff=cutoff, outpickle=picklefile, fileprefix=outfile, totalruns=1000)
    convergence_test(cutoff=cutoff,
                     outpickle=picklefile, 
                     fileprefix=,
                     totalruns=totalruns, 
                     maxgroupsize=totalruns,
                     outfile=clusteroutfile)


def main():
    if len(sys.argv) < 2:
        print "usage %s: <run> <steps>" % sys.argv[0]
        sys.exit()
    n = sys.argv[1]
    s = sys.argv[2]
    parameter_set=1
    kparameter_set=4
    run_gul(n, int(s), kparameter_set=kparameter_set, parameter_set=parameter_set)


if __name__ == "__main__":

    #############################################################
    # FILENAMES
    ligandfile = 'gulonate_test/gul_candidate_ligands.txt'
    seafile = 'gulonate_test/sea_condensed.csv'
    thermofluor='gulonate_test/thermofluor_hits.txt'
    cmetabfile = 'background_data/central_metab.smi'
    evidencefile = 'gulonate_test/dehydratase_substrate.txt'
    
    datafile = 'gulonate_test/gulonate_data.h5'
    fileprefix = 'output/gulonate_'
    randomh5 = 'output/gulonate_random.h5'
    randompicklefile = 'output/gulonate_random.pickle'
    picklefile = 'output/gul.pickle'
    clusteroutfile = 'output/gul_clusters_out.txt'

    proteins = ['1', '2', '3', '4', '5']
    
    #############################################################

    parser = OptionParser()
    parser.add_option('-e', action='store', dest='function_call',
                      help='Options: setup, run, random, analyze')
    parser.add_option('-n', '--run_number', action='store', dest='run_number')
    parser.add_option('-s', '--steps', action='store', type='int', dest='num_steps')

    options, args = parser.parse_args()

    if options.function_call == 'setup':
        set_up_gulonate(ligandfile, seafile, datafile, evidencefile, thermofluorfile, cmetabfile)
    elif options.function_call == 'run':
        n = options.run_number
        s = options.num_steps
        if (n is None) or (s is None):
            if len(args) > 1:
                n = args[0]
                s = int(args[1])
        
        run_gul(n, int(s), proteins, ligandfile, datafile, fileprefix, kparameter_set=4,
                                                                        parameter_set=1)
    elif options.function_call == 'random':
        s = options.num_steps
        if s is None:
            if len(args) > 0:
             s = int(args[0])

        # Run random sampling to get background statistics on distribution
        gulonate_sample_random(int(s), proteins, ligandfile, datafile, randomh5, randompicklefile)

    elif options.function_call == 'analyze':
        analyze(fileprefix, picklefile, randompicklefile, clusteroutfile)

    else:
        test_prediction(datafile)


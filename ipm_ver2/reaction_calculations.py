from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from itertools import product
from itertools import permutations
from copy import copy
from tables import open_file
import numpy as np
import pandas as pd
#import pybel  # Remove dependency on pybel eventually?
import sys
import math

def spam(n):
    out = []
    for perm in getPerms(n):
        elem = [int(i) for i in list(perm)]
        out.append(elem)
    return out

def getPerms(n):
    for i in getCandidates(n):
        for perm in set(permutations(i)):
            yield ''.join(perm)

def getCandidates(n):
    for i in range(0, n+1):
        res = "1"*i+"0" *(n-i)
        yield res

# Adapted code from
# https://github.com/rdkit/rdkit/issues/626
def GetStereoIsomers(mol, maxNum=12):
    out = []

    chiralCenters = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

    # keep only unassigned chiral centers
    chiralCenters = [c for c in chiralCenters if c[1] == "?"]

    # return the molecule object if no unassigned centers were found
    if chiralCenters == []:
        return [mol]

    # All bit permutations with number of bits equals numbers of chiralCenters
    elements = spam(len(chiralCenters))

    for isoId, element in enumerate(elements):
        for centerId, i in enumerate(element):
            atomId = chiralCenters[centerId][0]
            if i == 0:
                mol.GetAtomWithIdx(atomId).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
            elif i == 1:
                mol.GetAtomWithIdx(atomId).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
        outmol = copy(mol)
        out.append(outmol)
        if len(out) >= maxNum:
            break
    return out


# Iterates over all values in list of lists
# From http://stackoverflow.com/questions/6340351/python-iterating-through-list-of-list
def traverse(o, tree_types=(list, tuple)):
    if isinstance(o, tree_types):
        for value in o:
            for subvalue in traverse(value, tree_types):
                yield subvalue
    else:
        yield o

def convert_sugar_forms(molecule):
    rxn1 = '[O:1]1[C:2]([C:8])[C:3][C:4][C:5][C:6]1[O:7][H]>>[H][O:1][C:2]([C:8])[C:3][C:4][C:5][C:6]=[O:7]'
    rxn2 = '[O:6]1[C:2]([O:7][H])([C:1])[C:3][C:4][C:5]1>>[C:1][C:2](=[O:7])[C:3][C:4][C:5][O:6][H]'
    rxn3 = '[O:1]1[C:2]([H:8])[C:3][C:4][C:5][C:6]1[O:7][H]>>[H][O:1][C:2]([H:8])[C:3][C:4][C:5][C:6]=[O:7]'
    rxn4 = '[O:6]1[C:2]([O:7][H])([H:1])[C:3][C:4][C:5]1>>[H:1][C:2](=[O:7])[C:3][C:4][C:5][O:6][H]'
    sugarrxns = [AllChem.ReactionFromSmarts(rxn1), 
                 AllChem.ReactionFromSmarts(rxn2),
                 AllChem.ReactionFromSmarts(rxn3),
                 AllChem.ReactionFromSmarts(rxn4)]
    rxnproducts = [molecule]
    seen = set()
    seen.add(Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True))    
    for sugarrxn in sugarrxns:
        prods = sugarrxn.RunReactants((molecule,))
        for p in traverse(prods):
            
            smilesprod = Chem.MolToSmiles(p, isomericSmiles=True)
            if smilesprod not in seen:
                pmol = Chem.MolFromSmiles(smilesprod)
                pmol = Chem.AddHs(pmol)
                seen.add(smilesprod)
                rxnproducts.append(pmol)
    return rxnproducts

def getMaxTC(fps1, fps2):
    if fps1 == [] or fps2 == []:
        return 0.
    tcs = []
    for pairfps in product(fps1, fps2):
        tc = DataStructs.TanimotoSimilarity(pairfps[0], pairfps[1])
        if tc == 1.0:
            return tc
        tcs.append(tc)
    return np.max(tcs)

def adjustment_factor_for_sea_data(seadf, proteins, maxcutoff=50.):
    seacopy = seadf.copy()
    for i in proteins:
        sea_ii = seadf.get_value(i, i)
        for j in proteins:
            sea_jj = seadf.get_value(j, j)
            normfactor = (1./maxcutoff)*min(sea_ii, sea_jj)
            seacopy[j][i] = seadf.get_value(i, j)*normfactor
    return seacopy

####################################################################
#
# Calculate the sea score for a pathway using the scores from the dataframe
# Input: pathway - a list of proteins in order of the pathway
#        df - a pandas dataframe with the sea scores
# Output: sea score
def calc_sea_score(pathway, df):
    scores = [df.get_value(pathway[p], pathway[p+1]) for p in range(len(pathway) - 1)]
    return sum(scores)


class ligandData(object):
    
    def __init__(self, datafile=None, smilesfile=None, smiles={},
                 molecules={}, delimiter=None, linearforms=True, keepmultiples=False):
        self.datafile = datafile
        self.smiles = smiles
        self.molecules = molecules
        self.fingerprints = {}
        self.linearforms = linearforms
        self.keepmultiples = keepmultiples
        self.linearmolecules = {}
        if smilesfile is not None:
            self.read_in_molecules(smilesfile=smilesfile,
                                   delimiter=delimiter,
                                   linearforms=linearforms)
            print 'Number of molecules: %d' % len(self.molecules)
            self.make_fingerprints()

    # Possibly deprecated, need to check    
    def read_in_smiles(self, smilesfilename, delimiter=None):
        with open(smilesfilename, 'r') as handle:
            smiles = {}
            lines = handle.readlines()
            for line in lines:
                if delimiter is not None:
                    fields = line.strip().split(delimiter)
                else:
                    fields = line.strip().split()
                fields = [f.strip() for f in fields]
                #smiles_string = pybel.readstring('smi', fields[0]).write('can')
                smiles_string = Chem.MolToSmiles(Chem.MolFromSmiles(fields[0]), True)
                smiles[fields[1]] = smiles_string
            self.smiles = smiles


    def read_in_molecules(self, smilesfile, delimiter=None, linearforms=True):
        with open(smilesfile, 'r') as handle:
            molecules = {}
            lines = handle.readlines()
            smilesdict = {}
            linearmolecules = {}
            for line in lines:
                if delimiter is not None:
                    fields = line.strip().split(delimiter)
                else:
                    fields = line.strip().split()
                if self.keepmultiples:
                    identifier = fields[1].strip()
                else:
                    identifier = fields[1].split('_')[0]
                
                smiles = fields[0].strip()
                mol = Chem.MolFromSmiles(smiles)

                #if linearforms:
                #    # Convert to linear forms
                #    identifier = fields[1].strip()
                #    linearmols = convert_sugar_conversions(mol)
                #else:
                #    # Preserve form from database
                #    identifier = fields[1].split('_')[0]
                #    linearmols = [mol]

                #if identifier not in molecules.keys():
                #    if len(linearmols) > 1:
                #        molecules[identifier] = linearmols[1]
                #        smiles = Chem.MolToSmiles(linearmols[1], isomericSmiles=True)
                #    else:
                #        molecules[identifier] = mol
                #    smiles = pybel.readstring('smi', smiles).write('can')
                #    smilesdict[identifier] = smiles
                if identifier not in molecules.keys():
                    molecules[identifier] = mol
                    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                    smilesdict[identifier] = smiles
                    m = Chem.AddHs(mol)
                    if self.linearforms:
                        ms = convert_sugar_forms(m)
                    else:
                        ms = [m]
                    linearmolecules[identifier] = ms
            self.molecules = molecules
            self.smiles = smilesdict
            self.linearmolecules = linearmolecules

    def make_fingerprints(self):
        fingerprints = {}
        for k in self.molecules.keys():
            m = self.molecules[k]
            m = Chem.AddHs(m)
            if self.linearforms:
                ms = convert_sugar_forms(m)
            else:
                ms = [m]

            refs = []
            for refmol in ms:
                refs.extend(GetStereoIsomers(refmol))
            refs = [Chem.RemoveHs(rdkmol) for rdkmol in refs]
            refs = [Chem.MolFromSmiles(Chem.MolToSmiles(rdkmol, isomericSmiles=True)) for rdkmol in refs]
            fps = [AllChem.GetMorganFingerprintAsBitVect(rdkmol, 2, useChirality=True) for rdkmol in refs]
            fingerprints[k] = fps
        self.fingerprints = fingerprints

    def compute_match(self):
        idxs = []
        dictlist = []
        for f1 in self.fingerprints.keys():
            ds = {}
            for f2 in self.fingerprints.keys():
                fpA = self.fingerprints[f1]
                fpB = self.fingerprints[f2]
                #sim = DataStructs.TanimotoSimilarity(fpA, fpB)
                sim = getMaxTC(fpA, fpB)
                ds[f2] = sim
            idxs.append(f1)
            dictlist.append(ds)
        return idxs, dictlist

    def compare_metabolites_to_list(self, reffile):
        ref_fps = []
        with open(reffile, 'r') as refhandle:
            lines = refhandle.readlines()
            for line in lines:

                smiles = line.strip()
                m = Chem.MolFromSmiles(smiles)
                m = Chem.AddHs(m)
                if self.linearforms:
                    ms = convert_sugar_forms(m)
                else:
                    ms = [m] 
                refs = []
                for refmol in ms: 
                    refs.extend(GetStereoIsomers(refmol))
                refs = [Chem.RemoveHs(rdkmol) for rdkmol in refs]
                refs = [Chem.MolFromSmiles(Chem.MolToSmiles(rdkmol, isomericSmiles=True)) for rdkmol in refs]
                fps = [AllChem.GetMorganFingerprintAsBitVect(rdkmol, 2, useChirality=True) for rdkmol in refs]
                ref_fps.append(fps)
        tcdict = {}
        for fpkey in self.fingerprints.keys():
            maxtc = 0.0 
            for ref_fp in ref_fps:
                tc = getMaxTC(ref_fp, self.fingerprints[fpkey])
                maxtc = max(tc, maxtc)
            tcdict[fpkey] = maxtc
        return tcdict

    def reaction_sets(self, reactions):
        idxs = []
        count = 0
        edges = []

        for k in self.molecules.keys():
            count += 1
            idxs.append(k)
            fps = []
            rs = self.linearmolecules[k]
            for reaction in reactions:
                rdk_rxn = AllChem.ReactionFromSmarts(reaction)
                
                for r in rs:
                    rxn_products = rdk_rxn.RunReactants((r,))
                    
                    for rxn_p in traverse(rxn_products):
                        try:
                            rxn_p_noH = Chem.RemoveHs(rxn_p)
                            pstereos = GetStereoIsomers(rxn_p_noH)
                            smis = [Chem.MolToSmiles(p, isomericSmiles=True) for p in pstereos]
                            pmols = [Chem.MolFromSmiles(s) for s in smis]
                            cfps = [AllChem.GetMorganFingerprintAsBitVect(pm, 2, useChirality=True) for pm in pmols]
                            fps.extend(cfps)
                
                        except Exception, e:
                            sys.stderr.write('ERROR: %s\n\t%s' % (str(e), k))

                            
            if len(fps) > 0:
                for j in self.fingerprints.keys():
                    score = getMaxTC(self.fingerprints[j], fps)
                    if score == 1.0:
                        edges.append([k, j])
        return edges

    # Input: list of reactions
    # Output: Pandas DataFrame with molecules as indices and columns
    #         and tanimoto coefficients as values
    def reaction_calculations(self, reactions):
        list_of_dicts = []
        idxs = []
        count = 0
        edges = []

        for k in self.molecules.keys():
            count += 1
            fps = []
            rs = self.linearmolecules[k]
            
            # Perform each reaction in list on the molecule
            for reaction in reactions:
                rdk_rxn = AllChem.ReactionFromSmarts(reaction)

                for r in rs:
                    rxn_products = rdk_rxn.RunReactants((r,))

                    for rxn_p in traverse(rxn_products):
                        try:
                            rxn_p_noH = Chem.RemoveHs(rxn_p)
                            pstereos = GetStereoIsomers(rxn_p_noH)
                            smis = [Chem.MolToSmiles(p, isomericSmiles=True) for p in pstereos]
                            pmols = [Chem.MolFromSmiles(s) for s in smis]
                            cfps = [AllChem.GetMorganFingerprintAsBitVect(pm, 2, useChirality=True) for pm in pmols]
                            fps.extend(cfps)

                        except Exception, e:
                            sys.stderr.write('ERROR: %s\n\t%s' % (str(e), k))

            # Calculate tanimoto coefficients between molecule and all other molecules
            dict2 = {}
            for j in self.fingerprints.keys():
                fpr = self.fingerprints[j]
                score = getMaxTC(fpr, fps)
                dict2[j] = score

            list_of_dicts.append(dict2)
            idxs.append(k)

        df = pd.DataFrame(list_of_dicts, index=idxs)
        return df

    def write_possible_reaction_sets(self, reactiondict, textfile):
        with open(textfile, 'w') as handle:
            for key in reactiondict.keys():
                print key
                reactions = reactiondict[key]
                edges = self.reaction_sets(reactions)
                for edge in edges:
                    handle.write('%s, %s, %s\n' % (key, edge[0], edge[1]))

    def set_dock_scores_pd(self, scoresfile, reactionkey, idprefix='ZINC', id_column=1, score_column=7,
                           delimiter=None):
        ligdict = {}

        with open(scoresfile, 'r') as handle:
            lines = handle.readlines()
            allscores = []
            for line in lines:
                if delimiter:
                    fields = line.split(delimiter)
                else:
                    fields = line.split()
                if len(fields) > 1:
                    val = float(fields[score_column-1])
                    allscores.append(val)
            nar = np.array(allscores)
            mean = np.mean(nar)
            std = np.std(nar)
            maxscore = np.max(nar)

            ms = []
            for line in lines:
                if delimiter:
                    fields = line.split(delimiter)
                else:
                    fields = line.split()
                if len(fields) > 1:
                    molid = fields[id_column-1].strip()
                    #molid = '%s%s' % (idprefix, molid)
                    
                    if not self.keepmultiples:
                        molid = molid.split('_')[0]

                    val = float(fields[score_column-1])
                    score = -(val-mean)/std
                    if molid in self.molecules.keys():
                        ms.append(molid)
                        if molid in ligdict.keys():
                            score = np.max([score, ligdict[molid]])
                        ligdict[molid] = score

        for i in self.molecules.keys():
            if i not in ligdict.keys():
                ligdict[i] = 0
        return ligdict

    # returns a dictionary {molecule id: docking score} for a given enzyme
    def set_scores_pd(self, scoresfile, reactionkey, delimiter=';', dmean=None, dstd=None):
        ligdict = {}
        scores = []
        with open(scoresfile, 'r') as handle:
            lines = handle.readlines()
        if not dmean:
            for line in lines:
                if not line.startswith('#'):
                    fields = line.split(delimiter)
                    scores.append(float(fields[1]))
            nar = np.array(scores)
            dmean = np.mean(nar)
            dstd = np.std(nar)
            maxscore = np.max(nar)
        else:
            maxscore = dmean

        for line in lines:
            fields = line.split(delimiter)
            if len(fields) > 1:
                #ligandid = fields[0].split('_')[0]
                if self.keepmultiples:
                    ligandid = fields[0]
                else:
                    ligandid = fields[0].split('_')[0]
                score = float(fields[1])
                if line.startswith('#'):
                    ligandid = ligandid.strip('#')
                    ligdict[ligandid] = -(maxscore - dmean)/dstd
                else:
                    zscore = -(score - dmean)/dstd
                    if ligandid in ligdict.keys():
                        zscore = np.max([zscore, ligdict[ligandid]])
                    ligdict[ligandid] = zscore
        for i in self.molecules.keys():
            if i not in ligdict.keys():
                ligdict[i] = 0.
        return ligdict

    def get_sea_dict_from_file(self, seafile):
        cutoff = 50
        enz_idxs = []
        allscores_nonmatch = []
        dictlist = []
        with open(seafile, 'r') as handle:
            lines = handle.readlines()

            # Remove header
            if lines[0].startswith('target'):
                lines.pop(0)

            # Get enzyme ids and best and worst scores
            for line in lines:
                fields = line.split(',')
                if len(fields) > 1 and fields[0] != fields[1]:
                    try:
                        num = float(fields[2])
                        val = -np.log10(float(fields[2]))
                        if math.isnan(val):
                            val = 0.
                        if val > cutoff:
                            val = cutoff
                        #if val < 0:
                        #    val = 0.
                        allscores_nonmatch.append(val)
                        enz_idxs.append(fields[0])
                    except Exception as e:
                        print "Error reading in e-value: %s" % (fields[2])
                        print "Errors: ", sys.exc_info()[0], e

            nar = np.array(allscores_nonmatch)
            best = nar.max()
            worst = nar.min()

            enz_idxs = list(set(enz_idxs))
            ds = {}
            for e in enz_idxs:
                ds[e] = {}

            for line in lines:
                fields = line.split(',')
                val = 1.
                if len(fields) > 1:
                    val = -np.log10(float(fields[2]))
                    if math.isnan(val):
                        val = 0.
                    if val > cutoff:
                        val = cutoff
                    #if val < 0.:
                    #    val = 0.
                    #val = (val-worst)/(best-worst)
                    #val = min(val, 1.0)
                ds[fields[0]][fields[1]] = val

            for e in enz_idxs:
                dictlist.append(ds[e])

        return enz_idxs, dictlist

    ######################################################################
    # Docking scores
    # 
    # Save a pandas dataframe associating docking scores between 
    # molecules and enzymes
    # Input: Dictionary of enzymes (keys) and filenames (values)
    #
    ######################################################################
    def save_docking_scores(self, dockfiles, idprefix='ZINC', id_column=2, score_column=7,
                            delimiter=None):
        ld = []
        enz = []
        for key in dockfiles.keys():
            ld.append(self.set_dock_scores_pd(dockfiles[key], key, 
                                              idprefix=idprefix, id_column=id_column,
                                              score_column=score_column, delimiter=delimiter))
            enz.append(key)
        dockscores = pd.DataFrame(ld, index=enz)
        dockscores.to_hdf(self.datafile, 'dock')
        print 'Docking scores saved'

    ######################################################################
    # Chemical transformation/reaction scores
    #
    # Save a pandas panel associating reaction scores with enzymes, 
    # and pairs of molecules
    # Input: Dictionary of enzymes (keys) and lists of SMARTs pattern
    #        chemical transformations (values)
    #
    ######################################################################
    def save_reaction_scores(self, reactiondict):        
        rxn_calc_dict = {}
        for key in reactiondict.keys():
            rxn_calc_dict[key] = self.reaction_calculations(reactiondict[key])
            
        rxn_panel = pd.Panel(rxn_calc_dict)
        rxn_tc_mean = np.mean(rxn_panel.values)
        rxn_tc_std = np.std(rxn_panel.values)

        zscore_fxn = lambda x: (x - rxn_tc_mean)/rxn_tc_std
        
        rxn_panel = rxn_panel.apply(zscore_fxn)
        rxn_panel.to_hdf(self.datafile, 'rsim')
        print 'Chemical transformation/reaction similiarity scores saved'


    ######################################################################
    # SEA scores
    # 
    # Save computed z-scores for SEA pathway scores
    # Input: File with pairwise SEA scores and
    #        list of pathway proteins
    #
    ######################################################################
    def save_sea_scores(self, seafile, proteinlist=[]):
        elist, dlist = self.get_sea_dict_from_file(seafile)
        seavals = pd.DataFrame(dlist, index=elist)
        seavals = adjustment_factor_for_sea_data(seavals, elist)
        seavals.to_hdf(self.datafile, 'sea')
        if len(proteinlist) == 0:
            proteinlist = elist[:]
        seascores = [calc_sea_score(pathway, seavals) for pathway in permutations(proteinlist, len(proteinlist))]
        smu = np.mean(seascores)
        ssd = np.std(seascores)

        with open_file(self.datafile, 'a') as h5data:
            h5data.root.sea.block0_values.attrs.seamean = smu
            h5data.root.sea.block0_values.attrs.seastd = ssd

        print 'Sea scores saved'

    def save_thermofluor_scores(self, tffile, transporter='0'):
        tcdict = self.compare_metabolites_to_list(tffile)
        df = pd.DataFrame(tcdict, index=[transporter])
        meantc = np.mean(df.values)
        stdtc = np.mean(df.values)
        zscorefxn = lambda x: (x - meantc)/stdtc
        df = df.apply(zscorefxn)
        df.to_hdf(self.datafile, 'tfluor')

        print 'Thermofluor scores saved'

    def save_evidence_scores(self, evidencedict):
        enzdict = {}
        for e in evidencedict.keys():
            tcdict = self.compare_metabolites_to_list(evidencedict[e])
            series = pd.Series(tcdict)
            meantc = np.mean(series.values)
            stdtc = np.std(series.values)
            zscorefxn = lambda x : (x - meantc)/stdtc
            series = series.apply(zscorefxn)
            enzdict[e] = series
        df = pd.DataFrame(enzdict)
        df.to_hdf(self.datafile, 'evidence')

        print 'Evidence scores saved'

    def save_central_metabolism_endpt_scores(self, cmetabfile):
        tcdict = self.compare_metabolites_to_list(cmetabfile)
        series = pd.Series(tcdict)
        meantc = np.mean(series.values)
        stdtc = np.std(series.values)
        zscorefxn = lambda x: (x - meantc)/stdtc
        series = series.apply(zscorefxn)
        series.to_hdf(self.datafile, 'cmetab')

        print 'Central metabolism scores saved'

    def save_smiles(self):
        smileseries = pd.Series(self.smiles)
        smileseries.to_hdf(self.datafile, 'smiles')

    def save_reaction_sets(self, reactiondict):
        rows = []
        for key in reactiondict.keys():
            rxns = reactiondict[key]
            print key
            for r in rxns:
                print '', r
            edges = self.reaction_sets(rxns)
            curr_row = [(key, edge[0], edge[1]) for edge in edges]
            rows.extend(curr_row)
        df = pd.DataFrame(rows, columns=['enzyme', 'substrate', 'product'])
        df.to_hdf(self.datafile, 'reactions')

        print 'Reaction sets saved'

        


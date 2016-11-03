#!/usr/bin/python

import math
import numpy as np
import pandas as pd
import pickle
import pybel
import re

from itertools import product
from itertools import permutations
from openeye.oechem import *
from openeye.oegraphsim import *
from openeye.oeomega import OEFlipper
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from tables import open_file
from rxngraphs import *



def convert_to_zscore_array(alist):
    narray = np.array(alist)
    mu = narray.mean()
    sigma = narray.std()
    zscore_np = [(x-mu)/sigma for x in narray]
    return zscore_np

def make_molecule(smi):
    mol = OEGraphMol()
    OEParseSmiles(mol, smi)
    #OEAssignImplicitHydrogens(mol)
    OEAssignFormalCharges(mol)
    return mol

def check_molecule(m, mol):
    match = False
    for other in mol:
        match = match or OEExactGraphMatch(m, other)
    return match

def smilesToRDKitMol(smiles):
    rdkmol = Chem.MolFromSmiles(smiles)
    Chem.SanitizeMol(rdkmol)
    AllChem.AssignStereochemistry(rdkmol)
    return rdkmol

def OEMolToRDKitMol(mol):
    smiles = OECreateIsoSmiString(mol) 
    return smilesToRDKitMol(smiles)

def smilesToMorganFingerprint(smiles):
    smiles = pybel.readstring('smi', smiles).write('can')
    rdkmol = smilesToRDKitMol(smiles)
    fp = AllChem.GetMorganFingerprintAsBitVect(rdkmol, 4, useChirality=True)
    return fp

def standardizeMol(mol):
    OEAssignAromaticFlags(mol)
    OEAssignFormalCharges(mol)
    return mol

def OEMolToMorganFingerprint(mol):
    smiles = OECreateIsoSmiString(mol)
    rdkmol = Chem.MolFromSmiles(smiles)
    if rdkmol:
        Chem.SanitizeMol(rdkmol)
        AllChem.AssignStereochemistry(rdkmol)
        fp = AllChem.GetMorganFingerprintAsBitVect(rdkmol, 4, useChirality=True)
        return fp
    else:
        return False

def useFlipper(mol):
    stereomols = [m for m in OEFlipper(mol, 4)]
    return stereomols

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

def do_sugar_conversions(p):    
    sugarconversion1 = '[#1,#6:1][C:2](=[O:8])[CH:3][CH:4][CH:5][CH:6][O:7][H]>>[O:7]1[C:2]([O:8][H])([#1,#6:1])[C:3][C:4][C:5][C:6]1'
    sugarconversion2 = '[#1,#6:1][C:2](=[O:7])[CH:3][CH:4][CH:5][O:6][H]>>[O:6]1[C:2]([O:7][H])([#1,#6:1])[C:3][C:4][C:5]1'
    revsugarconversion1 = '[O:7]1[C:2]([O:8][H])([#1,#6:1])[C:3][C:4][C:5][C:6]1>>[#1,#6:1][C:2](=[O:8])[CH:3][CH:4][CH:5][CH:6][O:7][H]'
    revsugarconversion2 = '[O:6]1[C:2]([O:7][H])([#1,#6:1])[C:3][C:4][C:5]1>>[#1,#6:1][C:2](=[O:7])[CH:3][CH:4][CH:5][O:6][H]'
    sixlinear= '[H][O:1][C:2]([#1,#6:8])[C:3][C:4][C:5][C:6]=[O:7]'
    sixcyclic = '[O:1]1[C:2]([#1,#6:8])[C:3][C:4][C:5][C:6]1[O:7][H]'
    revsugarconversion1 =  '%s>>%s' % (sixcyclic, sixlinear)
    #sugarrxns = [sugarconversion1, sugarconversion2, revsugarconversion1, revsugarconversion2]
    sugarrxns = [revsugarconversion1, revsugarconversion2]
    sugarLibs = [OELibraryGen(s) for s in sugarrxns]
    rxnproducts = [p]
    seen = set()
    seen.add(OECreateIsoSmiString(p))
    for slib in sugarLibs:
        slib.SetValenceCorrection(True)
        m = slib.AddStartingMaterial(p, 0)
        if m:
            for k in slib.GetProducts():
                smilesk = OECreateIsoSmiString(k)
                if smilesk not in seen:
                    mol = OEGraphMol()
                    OEParseSmiles(mol, smilesk)
                    OEPerceiveChiral(mol)
                    rxnproducts.append(mol)
                    seen.add(smilesk)
    return rxnproducts

def set_up_data(smilesfile, delimiter_in_smiles, dockfiles, reactions, 
                seafile, filename, useMorganFP=False, proteins=None, compute_chemsim=True):
    enz = []
    ld = []
    lig = ligandData(smilesfile=smilesfile, delimiter=delimiter_in_smiles,
                     reactions=reactions, useMorganFP=useMorganFP)

    for key in dockfiles.keys():
        if delimiter_in_smiles is not None:
            ld.append(lig.set_scores_pd(dockfiles[key], key))
        else:
            ld.append(lig.set_dock_scores_pd(dockfiles[key], key))
        enz.append(key)
    dockscores = pd.DataFrame(ld, index=enz)

    if seafile is not None:
        elist, dlist = lig.get_sea_dict_from_file(seafile)
        seavals = pd.DataFrame(dlist, index=elist)
        print 'Adjusting SEA Data'
        seavals = adjustment_factor_for_sea_data(seavals, elist)
        seavals.to_hdf(filename, 'sea')
        print "Getting SEA stats for normalization"
        # Normalize for sea pathway scores
        if proteins is None:
            proteins = elist
        
        seascores = [calc_sea_score(pathway, seavals) for pathway in permutations(proteins, len(proteins))]
        smu = np.mean(seascores)
        ssd = np.std(seascores)
        
        with open_file(filename, 'a') as h5data:
            h5data.root.sea.block0_values.attrs.seamean = smu
            h5data.root.sea.block0_values.attrs.seastd = ssd

    if reactions is not None:
        reactdict = {}
        for key in reactions.keys():
            reactdict[key], edges = lig.react_pd(key)

        reactscores = pd.Panel(reactdict)
        rsmean = np.mean(reactscores.values)
        rsstd = np.std(reactscores.values)
        zscorefxn = lambda x: (x - rsmean) / rsstd
        reactscores = reactscores.apply(zscorefxn)
        reactscores.to_hdf(filename, 'rsim')
        #elist, dlist = lig.get_sea_zscores_dict_from_file(seafile)
    
    smileseries = pd.Series(lig.smiles)
    smileseries.to_hdf(filename, 'smiles')
    if compute_chemsim:
        llist, ldict = lig.compute_match()
        ligscores = pd.DataFrame(ldict, index=llist)
        ligscores.to_hdf(filename, 'chsim')
    dockscores.to_hdf(filename, 'dock')
       
    return lig


def find_possible_reaction_sets(smilesfile, delimiter_in_smiles, reactions, textfile, useMorganFP=False):
    enz = []
    ld = []
    lig = ligandData(smilesfile, delimiter=delimiter_in_smiles,
                     reactions=reactions, useMorganFP=useMorganFP)
    with open(textfile, 'w') as handle:
        for key in reactions.keys():
            r, edges = lig.react_pd(key)
            for edge in edges:
                handle.write('%s, %s, %s\n' % (key, edge[0], edge[1]))

# Calculates tanimoto coefficient between all pairs of ligands
# saves smiles and tanimoto coefficients to a hdf5 file in pandas
# dataframe or series and if provided, writes out to an output text file
def calculate_pairwise_tc(inputsmiles, delimiter, outputfile, outputtextfile=None):
    ld = ligandData(smilesfile=inputsmiles, useMorganFP=True, delimiter=delimiter)
    llist, ldict = ld.compute_match()
    ligscores = pd.DataFrame(ldict, index=llist)
    smileseries = pd.Series(ld.smiles)
    smileseries.to_hdf(outputfile, 'smiles')
    ligscores.to_hdf(outputfile, 'chsim')
    
    if outputtextfile is not None:
        with open(outputtextfile, 'w') as handle:
            for i, j in ligscores.iterrows():
                for h, k in j.iteritems():
                    handle.write('%s\t%s\t%.10f\n' % (i, h, k))

def adjustment_factor_for_sea_data(seadf, proteins, maxcutoff=50.):
    seacopy = seadf.copy()
    for i in proteins:
        sea_ii = seadf.get_value(i, i)
        for j in proteins:
            sea_jj = seadf.get_value(j, j)
            normfactor = (1./maxcutoff)*min(sea_ii, sea_jj)
            seacopy[j][i] = seadf.get_value(i, j)*normfactor
    return seacopy

# Calculate the sea score for a pathway using the scores from the dataframe
# Input: pathway - a list of proteins in order of the pathway
#        df - a pandas dataframe with the sea scores
# Output: sea score
def calc_sea_score(pathway, df):
    scores = [df.get_value(pathway[p], pathway[p+1]) for p in range(len(pathway) - 1)]
    return sum(scores)

class ligandData(object):
    def __init__(self, smilesfile=None, smirksfile=None, reactions={}, smiles={}, 
                 molecules={}, useMorganFP=False, delimiter=None):
        self.smiles = smiles
        self.molecules = molecules
        self.reactions = reactions
        self.useMorganFP = useMorganFP
        
        if smilesfile is not None:
            self.read_in_molecules(smilesfilename=smilesfile, delimiter=delimiter)
            print 'Number of molecules: %d' % len(self.molecules)
        if smirksfile is not None:
            self.read_in_smirks(smirksfile)
        self.make_fingerprints()
        

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
                smiles_string = pybel.readstring('smi', fields[0]).write('can')
                smiles[fields[1]] = smiles_string
            self.smiles = smiles    

    def read_in_molecules(self, smilesfilename, delimiter=None):
        with open(smilesfilename, 'r') as handle:
            handle = open(smilesfilename, 'r')
            molecules = {}
            lines = handle.readlines()
            smilesdict = {}
            for line in lines:
                if delimiter is not None:
                    fields = line.strip().split(delimiter)
                else:
                    fields = line.strip().split()
                smiles = fields[0].strip()
                #print fields[1]
                mol = OEGraphMol()
                OEParseSmiles(mol, smiles)
                #OEAssignImplicitHydrogens(mol)
                #OEAssignFormalCharges(mol)
                mol = standardizeMol(mol)

                # Convert to linear forms
                identifier = fields[1].strip()
                linearmols = do_sugar_conversions(mol)

                # Preserve form from database
                #identifier = fields[1].split('_')[0]
                #linearmols = [mol]

                if identifier not in molecules.keys():
                    if len(linearmols) > 1:
                        molecules[identifier] = linearmols[1]
                        smiles = OECreateIsoSmiString(linearmols[1])
                    else:
                        molecules[identifier] = mol
                    smiles = pybel.readstring('smi', smiles).write('can')
                    smilesdict[identifier] = smiles
            self.molecules = molecules
            self.smiles = smilesdict
    

    def read_in_smirks(self, smirksfilename):
        with open(smirksfilename, 'r') as handle:
            reactions = {}
            lines = handle.readlines()
            for line in lines:
                fields = lines.strip().split()
                reactions[fields[1]] = fields[2]
            self.reactions = reactions

    def make_fingerprints(self):
        fingerprints = {}
        if self.useMorganFP:
            for k in self.smiles.keys():
                s = self.smiles[k]
                fp = smilesToMorganFingerprint(s)
                fingerprints[k] = fp
        else:
            for k in self.molecules.keys():
                m = self.molecules[k]
                fp = OEFingerPrint()
                OEMakeFP(fp, m, OEFPType_Circular)
                fingerprints[k] = fp
        self.fingerprints = fingerprints

    def compute_match(self, table):
        row = table.row
        for f1 in self.fingerprints.keys():
            for f2 in self.fingerprints.keys():
                fpA = self.fingerprints[f1]
                fpB = self.fingerprints[f2]
                row['ligand1'] = f1
                row['ligand2'] = f2
                
                if self.useMorganFP:
                    row['similarity'] = DataStructs.TanimotoSimilarity(fpA, fpB)
                else:
                    row['similarity'] = OETanimoto(fpA, fpB)
                row.append()
        table.flush()

    def compute_match(self):
        idxs = []
        dictlist = []
        for f1 in self.fingerprints.keys():
            ds = {}
            for f2 in self.fingerprints.keys():
                fpA = self.fingerprints[f1]
                fpB = self.fingerprints[f2]
                if self.useMorganFP:
                    sim = DataStructs.TanimotoSimilarity(fpA, fpB)
                else:
                    sim = OETanimoto(fpA, fpB)
                ds[f2] = sim
            idxs.append(f1)
            dictlist.append(ds)
        return idxs, dictlist

    def react_pd(self, reactionkey):
       
        list_of_dicts = []
        idxs = []
        count = 0
        edges = []

        if self.useMorganFP:
            isofpdict = {}
            for j in self.molecules.keys():
                isomers = useFlipper(self.molecules[j])
                isofpdict[j] = [OEMolToMorganFingerprint(isomol) for isomol in isomers]

        for k in self.molecules.keys():
            count += 1
            idxs.append(k)
            dict2 = {}
            fps = []
            m = self.molecules[k]
            rs = self.reactions[reactionkey]
            for r in rs:
                lg = OELibraryGen(r)
                lg.AddStartingMaterial(m, 0)
                #lg.SetExplicitHydrogens(True)
                lg.SetValenceCorrection(True)
                for p in lg.GetProducts():
                    numparts, partlist = OEDetermineComponents(p)
                    if numparts == 1:
                        rxnproducts = do_sugar_conversions(p)
                        #rxnproducts = [p]
                        for rxnproduct in rxnproducts:
                            #OEAssignFormalCharges(p)
                            if self.useMorganFP:
                                isomers = useFlipper(rxnproduct)
                                for molisomer in isomers:
                                    fp = OEMolToMorganFingerprint(molisomer)
                                    if fp:
                                        fps.append(fp)
                            else:
                                fp = OEFingerPrint()
                                OEMakeFP(fp, rxnproduct, OEFPType_Circular)
                                fps.append(fp)
                    elif numparts > 1:
                        pred = OEPartPredAtom(partlist)
                        for i in range(1, numparts + 1):
                            pred.SelectPart(i)
                            partmol = OEGraphMol()
                            OESubsetMol(partmol, p, pred)
                            #OEAssignFormalCharges(partmol)
                            rxnproducts = do_sugar_conversions(partmol)
                            #rxnproducts = [partmol]
                            for rxnproduct in rxnproducts:
                                if self.useMorganFP:
                                    isomers = useFlipper(rxnproduct)
                                    for molisomer in isomers:
                                        fp = OEMolToMorganFingerprint(molisomer)
                                        if fp:
                                            fps.append(fp)
                                else:
                                    fp = OEFingerPrint()
                                    OEMakeFP(fp, rxnproduct, OEFPType_Circular)
                                    fps.append(fp)


            if self.useMorganFP:
                for j in isofpdict.keys():
                    score = getMaxTC(isofpdict[j], fps)
                    dict2[j] = score
                    if score == 1.0:
                        edges.append([k, j])
            else:
                for j in self.fingerprints.keys():
                    fpr = self.fingerprints[j]
                    tc = [OETanimoto(fpr, fpp) for fpp in fps]
                    tc.append(0.0)
                    score = max(tc)
                    dict2[j] = score

            list_of_dicts.append(dict2)

        x = pd.DataFrame(list_of_dicts, index=idxs)
        return x, edges


    def set_dock_scores_pd(self, scoresfile, reactionkey):
        ligdict = {}
        with open(scoresfile, 'r') as handle:
            lines = handle.readlines()
            allscores = []
            for line in lines:
                fields = line.split()
                if len(fields) > 1:
                    val = float(fields[6])
                    if val <= 0.0:
                        allscores.append(val)
            nar = np.array(allscores)
            mean = np.mean(nar)
            std = np.std(nar)
            max = np.max(nar)

            for i in self.molecules.keys():
                score = 0
                pattern = i.lstrip('ZINC')
                match = re.compile(pattern)
                for line in lines:
                    if re.search(match, line):
                        fields = line.split()
                        if len(fields) > 1:
                            val = float(fields[6])
                            if val <=0.0:
                                score = (val-mean)/std
                            else:
                                score = (max-mean)/std
                ligdict[i] = -score
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

    def get_sea_zscores_dict_from_file(self, seafile, zcol=4):
        enz_idxs = []
        allscores_non_match = []
        dictlist = []
        with open(seafile, 'r') as handle:
            lines = handle.readlines()
            lines.pop(0)
            for line in lines:
                fields = line.split(',')
                if len(fields) > 1 and fields[0] != fields[1]:
                    try:
                        enz_idxs.append(fields[0])
                    except Exception as e:
                        print "Error occured while  trying to convert pe value: %s to double " % (fields[2])
                        print "Errors :\n\t\t", sys.exec_info()[0], ":\t", e
            enz_idxs = list(set(enz_idxs))
            ds = {}
            for e in enz_idxs:
                ds[e] = {}
            for line in lines:
                fields = line.split(',')
                val = 1.
                if len(fields) > 1:
                    val = float(fields[zcol])
                    if math.isnan(val):
                        val = 0.
                ds[fields[0]][fields[1]] = val
            for e in enz_idxs:
                dictlist.append(ds[e])
        return enz_idxs, dictlist

    def compare_metabolites_to_list(self, reffile):
        ref_fps = []
        with open(reffile, 'r') as refhandle:
            lines = refhandle.readlines()
            for line in lines:
                smiles = line.strip()
                fp = smilesToMorganFingerprint(smiles)
                ref_fps.append(fp)
        tcdict = {}
        for fpkey in self.fingerprints.keys():
            maxtc = 0.0
            for ref_fp in ref_fps:
                tc = DataStructs.TanimotoSimilarity(ref_fp, self.fingerprints[fpkey])
                maxtc = max(tc, maxtc)
            tcdict[fpkey] = maxtc
        return tcdict

def main():
    if len(sys.argv) < 2:
        print 'usage %s: <run> <steps>' % sys.argv[0]
        sys.exit()
    n = sys.argv[1] # run number
    s = sys.argv[2] # number of steps
    kdo_sample_zscores(n, int(s))

if __name__ == "__main__":

    #kdo_set_up()
    #kdo_sample_random(1000000)
    #main()
    #get_kdo_score()
    pass

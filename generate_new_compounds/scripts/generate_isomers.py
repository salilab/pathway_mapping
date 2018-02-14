from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import numpy as np
import pandas as pd
import os,sys,string,math
import pprint
from rdkit.Chem.Fingerprints import FingerprintMols
from itertools import permutations
from copy import copy


#######file to write
tmxl=open("smi_all_isomers.csv",'w')
print >>tmxl,"#smiles strings"

#########file to read
filename=sys.argv[1]
sf=open(filename,'r')
smilist=[]
smidict={}
for i,ln in enumerate(sf.readlines()):
    line =ln.strip().split()
    smidict[line[1]]=line[2]

########list of smiles
for i in smidict.values():
    smilist.append(i) 

############define functions
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
    chiralCenters = [c for c in chiralCenters if c[1] == "?"]
    if chiralCenters == []:
        return [mol]
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

####################write out isomers
for i in smilist:
    print i,smilist.index(i),smidict.keys()[smidict.values().index(i)]
    keyn=smidict.keys()[smidict.values().index(i)]
    m=Chem.MolFromSmiles(i)
    refs=GetStereoIsomers(m)
    smis = [Chem.MolToSmiles(p, isomericSmiles=True) for p in refs]
    for j in range(0, len(smis)):
        print >> tmxl, str(keyn)+'_iso_'+str(j),smis[j]


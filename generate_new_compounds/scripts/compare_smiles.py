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
from collections import Counter


#######file to write
tmxl=open("common_smiles.txt",'w')
print >>tmxl,"#smiles ID1, smiles ID2,  smiles of ID1, smiles of ID2, Tanimoto coefficient"
tmxl1=open("cleaned_smilesa.txt",'w')
#########file to read
filename1=sys.argv[1]
filename2=sys.argv[2]
sf1=open(filename1,'r')
sf2=open(filename2,'r')
smidict={};smidict2={}
for i,ln in enumerate(sf1.readlines()):
    line =ln.strip().split()
    print line
    smidict[line[1]]=line[2]

for i,ln in enumerate(sf2.readlines()):
    line =ln.strip().split()
    smidict2[line[1]]=line[2]


############define functions
def getTC(fps1, fps2):
    tc = DataStructs.TanimotoSimilarity(fps1, fps2)
    return tc

def traverse(o, tree_types=(list, tuple)):
    if isinstance(o, tree_types):
        for value in o:
            for subvalue in traverse(value, tree_types):
                yield subvalue
    else:
        yield o

def getFP(smis):
    m=Chem.MolFromSmiles(smis,sanitize=True)
    rxn_p_noH = Chem.RemoveHs(m)
    smis = Chem.MolToSmiles(rxn_p_noH, isomericSmiles=True)
    pmols = Chem.MolFromSmiles(smis)
    return AllChem.GetMorganFingerprintAsBitVect(m, 15, useChirality=False)

def writedictofile(dictn,filen):
    for i,j in dictn.iteritems():
        print >> filen, dictn.keys().index(i), i, j


def compare_dicts(dict1,dict2):
    dict_return=dict1.copy();same_keys=[]
    for k,v in dict1.items():
        for k1,v1 in dict2.items():
            if (getTC(getFP(v),getFP(v1))==1 and k != k1):
                print k,k1,v,v1,getTC(getFP(v),getFP(v1))
                print >> tmxl, k,k1,v,v1,getTC(getFP(v),getFP(v1))
                same_keys.append(k)
                if k1 not in same_keys:
                    dict_return.pop(k1,None)
                break
    return dict_return

smidict=compare_dicts(smidict,smidict2)
writedictofile(smidict, tmxl1)


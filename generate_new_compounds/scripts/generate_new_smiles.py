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
tmxl=open("smi_all_rxn.csv",'w')
print >>tmxl,"#smiles strings"
tmxl1=open("productfile.csv",'w')
tmxl2=open("reactantfile.csv",'w')
tmxl3=open('smi_clean.csv','w')
#########file to read
filename=sys.argv[1]
number_of_runs=int(sys.argv[2])
sf=open(filename,'r')
smilist=[]
smidict={};smidict2={}
for i,ln in enumerate(sf.readlines()):
    line =ln.strip().split(';')
    smidict[line[1]]=line[0]
    
########maintain two dictionaries to cross check answers
for k,i in smidict.iteritems():
    smilist.append(i) 
    smidict2[smidict.keys().index(k)]=k
    
######### reactions list
###all reaction smirks checked by John Irwin
#########
phenyllactatedehydratase1=["n1cnc(N)c(c12)ncn2[C@H]3[C@H](O)[C@H](OP(=O)([O-])[O-])[C@H](O3)COP(=O)([O-])OP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCS-[C:1](=O)-[C:2](-[OD1])-[Ch:3]>>n1cnc(N)c(c12)ncn2[C@H]3[C@H](O)[C@H](OP(=O)([O-])[O-])[C@H](O3)COP(=O)([O-])OP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCS-[C:1](=O)-[C:2]=[C:3]"]
phenyllactatedehydratase2=["n1cnc(N)c(c12)ncn2[C@H]3[C@H](O)[C@H](OP(=O)([O-])[O-])[C@H](O3)COP(=O)([O-])OP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCS-[C:1](=O)-[C:2]=[C:3]>>n1cnc(N)c(c12)ncn2[C@H]3[C@H](O)[C@H](OP(=O)([O-])[O-])[C@H](O3)COP(=O)([O-])OP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCS-[C:1](=O)-[C:2](-[OD1])-[Ch:3]"]
betatoalpha=["n1cnc(N)c(c12)ncn2[C@H]3[C@H](O)[C@H](OP(=O)([O-])[O-])[C@H](O3)COP(=O)([O-])OP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCS-[C:1](=O)-[Ch2:2]-[C:3](-[OD1])>>n1cnc(N)c(c12)ncn2[C@H]3[C@H](O)[C@H](OP(=O)([O-])[O-])[C@H](O3)COP(=O)([O-])OP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCS-[C:1](=O)-[C:2](-[OD1])=[C:3]"]
alphatobeta=["n1cnc(N)c(c12)ncn2[C@H]3[C@H](O)[C@H](OP(=O)([O-])[O-])[C@H](O3)COP(=O)([O-])OP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCS-[C:1](=O)-[C:2](-[OD1])-[Ch2:3]>>n1cnc(N)c(c12)ncn2[C@H]3[C@H](O)[C@H](OP(=O)([O-])[O-])[C@H](O3)COP(=O)([O-])OP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCS-[C:1](=O)-[C:2]=[C:3](-[OD1])"]
dehydrogenase1=["n1cnc(N)c(c12)ncn2[C@H]3[C@H](O)[C@H](OP(=O)([O-])[O-])[C@H](O3)COP(=O)([O-])OP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCS-[C:1](=O)-[C:2]-[Ch:3](-[OD1])>>n1cnc(N)c(c12)ncn2[C@H]3[C@H](O)[C@H](OP(=O)([O-])[O-])[C@H](O3)COP(=O)([O-])OP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCS-[C:1](=O)-[C:2]-[CX3:3](=O)"]
dehydrogenase2=["n1cnc(N)c(c12)ncn2[C@H]3[C@H](O)[C@H](OP(=O)([O-])[O-])[C@H](O3)COP(=O)([O-])OP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCS-[C:1](=O)-[C:2]-[C:3](=O)>>n1cnc(N)c(c12)ncn2[C@H]3[C@H](O)[C@H](OP(=O)([O-])[O-])[C@H](O3)COP(=O)([O-])OP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCS-[C:1](=O)-[C:2]-[Ch:3](-[OD1])"]
dehydratase1=["n1cnc(N)c(c12)ncn2[C@H]3[C@H](O)[C@H](OP(=O)([O-])[O-])[C@H](O3)COP(=O)([O-])OP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCS-[C:1](=O)-[C:2]=[C:3]>>n1cnc(N)c(c12)ncn2[C@H]3[C@H](O)[C@H](OP(=O)([O-])[O-])[C@H](O3)COP(=O)([O-])OP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCS-[C:1](=O)-[Ch:2]-[C:3](-[OD1])"]
dehydratase2=["n1cnc(N)c(c12)ncn2[C@H]3[C@H](O)[C@H](OP(=O)([O-])[O-])[C@H](O3)COP(=O)([O-])OP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCS-[C:1](=O)-[Ch:2]-[C:3](-[OD1])>>n1cnc(N)c(c12)ncn2[C@H]3[C@H](O)[C@H](OP(=O)([O-])[O-])[C@H](O3)COP(=O)([O-])OP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCS-[C:1](=O)-[C:2]=[C:3]"]
acetyltransferase1=["n1cnc(N)c(c12)ncn2[C@H]3[C@H](O)[C@H](OP(=O)([O-])[O-])[C@H](O3)COP(=O)([O-])OP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCS-[C:1](=O)-[CH2:2]-[C:3](=O)-[C:4]>>n1cnc(N)c(c12)ncn2[C@H]3[C@H](O)[C@H](OP(=O)([O-])[O-])[C@H](O3)COP(=O)([O-])OP(=O)([O-])OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)NCCS-[C:3](=O)-[C:4]"]

###############reacion dictionary
reactions = {1:phenyllactatedehydratase1, 2:phenyllactatedehydratase2,3:betatoalpha,4:alphatobeta, 5:dehydrogenase1,6:dehydrogenase2, 7:dehydratase1, 8:dehydratase2,9:acetyltransferase1}

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

def writelisttofile(listn,filen):
    for i in listn:
        if i in smidict.values():
            print >> filen, i, smidict.keys()[smidict.values().index(i)]
        else:
            print >> filen, i, "UNK"

def writetofile(line, num, filen):
    if line in smidict.values():
        print >> filen,num, line, smidict.keys()[smidict.values().index(line)]
    else: 
        print >> filen, num,line, "UNK"

def writesmitofile(num,line,keyn,filen):
    print >> filen,num, line, keyn, smidict.keys()[smidict.values().index(line)]


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

def writedictofile(dictn,filen):
    for i,j in dictn.iteritems():
        print >> filen, dictn.keys().index(i), i, j

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

def compare_dicts(dict1,dict2):
    print "comparing dictionaries now.."
    dict_return=dict1.copy();same_keys=[]
    for k,v in dict1.items():
        for k1,v1 in dict2.items():
            if (getTC(getFP(v),getFP(v1))==1 and k != k1):
                print >> tmxl5, k,k1,v,v1,getTC(getFP(v),getFP(v1))
                same_keys.append(k1)
                if k not in same_keys:
                    dict_return.pop(k,None)
                    print "removed " + k+ " from your dictionary, because a similar compound "+ k1 + " exists"
                break
    return dict_return



####################get finger prints of original smiles

print "Total number of distinct smiles and hence FPs is:",  len(smilist)

##############transform a smiles string to a product based on reactions
for k in range(0,number_of_runs):
    print "This is the " + str(k) + "th run"
    for i in range(1,10):
        rdk_rxn = AllChem.ReactionFromSmarts(reactions[i][0])
        for r in smilist:
            m=Chem.MolFromSmiles(r,sanitize=True)
            m1 = Chem.RemoveHs(m)
            r2=Chem.MolToSmiles(m1)
            m2=Chem.MolFromSmiles(r2,sanitize=True)
            rxn_products = rdk_rxn.RunReactants((m2,))
            if len(rxn_products)>0:
                writetofile(Chem.MolToSmiles(rxn_products[0][0]),i,tmxl1)
                writetofile(r,i,tmxl2)
            for rxn_p in traverse(rxn_products):
                rxn_p_noH = Chem.RemoveHs(rxn_p)
                smis = Chem.MolToSmiles(rxn_p_noH, isomericSmiles=True)
                pmols = Chem.MolFromSmiles(smis)
                cfps = AllChem.GetMorganFingerprintAsBitVect(pmols, 15, useChirality=False)
                for j in range(0,len(smilist)+1):
                    if getTC(getFP(smilist[j]),cfps)==1.0:
                        break
                    elif j==(len(smilist)-1) and getTC(getFP(smilist[j]),cfps)<1.0:
                        keyc=smidict2[smilist.index(r)]+'_'+str(i)
                        smidict2[len(smilist)]=keyc
                        print len(smilist),keyc, smis
                        smidict[keyc]=smis
                        smilist.append(smis)
                        print >> tmxl, len(smilist),keyc, smis
                        
writedictofile(smidict,tmxl3)

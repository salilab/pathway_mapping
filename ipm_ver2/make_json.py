from pathway_restraints import *
from rxngraphs import modelGraph
import pandas as pd
import json
import pickle

def pathstring_to_components(pathstring, offset=0):
    pstart = 1
    lstart = 0
    if offset:
        pstart = 0
        lstart = 1
    nodes = pathstring.split(' -> ')
    proteins = nodes[pstart::2]
    ligands = nodes[lstart::2]
    return proteins, ligands


# Takes rxngraph object and returns a json format
def make_cygraph(mygraph, smilesdf):
    cyjs = {'edges':[], 'nodes':[], 'restraints':[], 'scores':[]}

    #edges = mygraph.rxngraph.interactions[:]
    edges = mygraph.interactions[:]

    for p in mygraph.proteins:
        node = {'label':p,
                'node type':'enzyme'}
        cyjs['nodes'].append(node)

    for l in mygraph.ligands:
        smiles = smilesdf.get_value(l).strip()
        node = {'label':l,
                'node type':'metabolite',
                'smiles':smiles}
        cyjs['nodes'].append(node)

    if mygraph.rxngraph.transporter:
        node = {'label':mygraph.rxngraph.transporter,
                'node type':'enzyme'} # transporter
        cyjs['nodes'].append(node)

        first_substrate = mygraph.rxngraph.get_start()
        transporter_substrate = mygraph.rxngraph.graph.node[first_substrate]['label']
        tedge = [mygraph.rxngraph.transporter, transporter_substrate]
        edges.append(tedge)

    for e in edges:
        edge = {'source':e[0], 'target':e[1]}
        cyjs['edges'].append(edge)

    scores = {}
    obj = mygraph.compute_scores()
    for res in mygraph.restraints:
        newrestraints = res.to_json(mygraph.rxngraph)
        if len(newrestraints) > 0:
            cyjs['restraints'].extend(newrestraints)
        resname = '%s score' % res.name
        scores[resname] = res.score
    scores['score'] = obj
    cyjs['scores'] = scores

    return cyjs


# Write out pathway models to a json file
# Input:  jsonfile - Name of json file
#         data - object of data tables
#         smiledf - Dataframe with ligand ids as indices and smiles strings
#         restraints - list of pathway restraint objects
#         np_path_array - list of paths string
#         transporter - (Default: None), protein id
def write_out_models_to_json(jsonfile, data, smilesdf, restraints, np_path_array, transporter=None):
    offset = 0
    if transporter is not None:
        offset = 1

    model_array = []
    for pathstring in np_path_array['strrepr']:
        proteins, ligands = pathstring_to_components(pathstring, offset=offset)

        if transporter is not None:
            proteins = proteins[offset::]

        edges = []
        for i in range(len(proteins)):
            edges.append([ligands[i], proteins[i]])
            if i < len(ligands) - 1:
                edges.append([proteins[i], ligands[i+1]])

        mgraph = modelGraph(data=data, proteins=proteins, ligands=ligands, interactions=edges)
        if transporter is not None:
            mgraph.add_transporter(transporter)
        mgraph.add_restraints(restraints)
        model = make_cygraph(mgraph, smilesdf)
        model_array.append(model)

    models = {'models': model_array}

    with open(jsonfile, 'w') as handle:
        json.dump(models, handle, indent=2)

def serine_models_to_json():
    jsonfile = 'json_graph_files/prototype_serine_biosyn_models.json'
    nummodels = 25

    SEAMEAN=-0.680521232682
    SEASTD=5.93823741616
    
    proteins = ['1', '2', '3', '4', '5']
    
    datafile = '/salilab/park2/scalhoun/cluster_output/ser/ser_data_sea.h5'
    d = dataTables(proteins=proteins)
    d.setBackground(seamean=SEAMEAN, seastd=SEASTD)
    d.readInSelectedTables(datafile, ['rsim', 'dock', 'sea'])
    smilesdf = pd.read_hdf(datafile, 'smiles')

    restraint_list = [dock_restraint(),
                      sea_restraint(),
                      reaction_restraint()]

    picklefile = 'ser_sea_2sd.pickle'
    with open(picklefile, 'r') as phandle:
        pens = pickle.load(phandle)
        np_path_array = pens.np_path_array[0:nummodels]

    write_out_models_to_json(jsonfile, d, smilesdf, restraint_list,
                             np_path_array)

def gulonate_models_to_json():
    jsonfile = 'json_graph_files/gulonate_models_test.json'
    protein_candidates = ['0', '1', '2', '3', '4', '5']
    datafile = '/salilab/park2/scalhoun/cluster_output/gulonate/gulonate_sea.h5'
    picklefile = 'gul_2sd.pickle'
    restraints = ['rsim', 'dock', 'sea', 'cmetab', 'tfluor', 'evid']

    d = dataTables(proteins=protein_candidates)
    d.readInSelectedTables(datafile, restraints)
    smilesdf = pd.read_hdf(datafile, 'smiles')
    with open(picklefile, 'r') as phandle:
        pens = pickle.load(phandle)
   
    res1 = reaction_restraint()
    res2 = sea_restraint()
    res3 = dock_restraint()
    res4 = evidence_restraint()
    res5 = cmetab_restraint()
    res6 = tfluor_restraint()
    restraint_list = [res1, res2, res3, res4, res5, res6] 

    write_out_models_to_json(jsonfile, d, smilesdf, restraint_list,
                             pens.np_path_array, transporter='0')
    
if __name__ == "__main__":
    gulonate_models_to_json()
    #serine_models_to_json()

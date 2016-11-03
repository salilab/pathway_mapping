#!/bin/python
import networkx as nx
from networkx import DiGraph
from networkx.algorithms import bipartite
import math
import numpy as np
import random as rand
import pandas as pd
from itertools import combinations

class rxnGraph(object):
    def __init__(self, datatables, proteins=None, ligands=None, interactions=None): #, dockdf=None, seadf=None, rsimpn=None):
        # Construct graph
        self.graph = DiGraph()
        self.data = datatables
        self.proteins = proteins  # list of nodes
        self.ligands = ligands    # list of nodes
        self.protein_members = [] # subset of possible proteins, members in graph
        self.protein_choices = [] # subset of possible proteins, not members in graph
        self.ligand_members = []  # subset of possible ligands, members in graph
        self.ligand_choices = []  # subset of possible ligands, not members in graph

        self.transporter = None
       
        self.protein_nodes = [] # moveable nodes
        self.ligand_nodes = []
        self.protein_allnodes = []
        self.ligand_allnodes = []
        self.fixed = []
        self.keep = []
        
        if proteins is not None and ligands is not None and interactions is not None:
            self.construct_graph(proteins, ligands, interactions)
            self.protein_members = proteins[:]
            self.ligand_members = ligands[:]

    def add_choices(self, protein_choices=[], ligand_choices=[]):
        self.ligand_choices = ligand_choices
        self.protein_choices = protein_choices

    def add_transporter(self, transporter):
        self.transporter = transporter
        #n = self.get_start()
        #self.transporter_substrate = self.graph.node[n]['label']

    def copy(self):
        newGraph = rxnGraph(self.data, proteins=self.proteins[:], ligands = self.ligands[:])
                            #dockdf=self.dockdf, seadf=self.seadf, rsimpn=self.rsimpn)
        newGraph.graph = self.graph.copy()
        #newGraph.constraints = self.constraints
        newGraph.ligand_members = self.ligand_members[:]
        newGraph.protein_members = self.protein_members[:]
        newGraph.ligand_choices = self.ligand_choices[:]
        newGraph.protein_choices = self.protein_choices[:]

        newGraph.protein_nodes = self.protein_nodes
        newGraph.ligand_nodes = self.ligand_nodes
        newGraph.protein_allnodes = self.protein_allnodes
        newGraph.ligand_allnodes = self.ligand_allnodes
        newGraph.fixed = self.fixed
        newGraph.keep = self.keep
        newGraph.transporter = self.transporter
        return newGraph

    def construct_graph_from_path(self, proteins, ligands):
        #self.graph.add_nodes_from(ligands, bipartite=0)
        #self.graph.add_nodes_from(proteins, bipartite=1)
        interactions = []
        for i in range(len(proteins)):
            interactions.append([ligands[i], proteins[i]])
            if i > 0:
                interactions.append([proteins[i-1], ligands[i]])
        if len(ligands) == len(proteins) + 1:
            interactions.append([proteins[-1], ligands[-1]])
        self.construct_graph(proteins, ligands, interactions)

    def construct_graph(self, proteins, ligands, interactions):
        # interactions: [[protein, ligand], [protein, ligand]]
        self.graph.clear()
        labelposdict = {}
        for i, p in enumerate(proteins):
            self.graph.add_node(i, label=p, bipartite=1)
            labelposdict[p] = i
            self.protein_allnodes.append(i)
        for j, l in enumerate(ligands):
            self.graph.add_node(j + len(proteins), label=l, bipartite=0)
            labelposdict[l] = j + len(proteins)
            self.ligand_allnodes.append(j + len(proteins))

        #for i in range(len(interactions)):
        #    node1 = labelposdict[interactions[i][0]]
        #    node2 = labelposdict[interactions[i][1]]
        #    self.graph.add_edge(node1, node2)
        [self.graph.add_edge(labelposdict[interactions[i][0]], labelposdict[interactions[i][1]]) for i in range(len(interactions))]
        
        self.proteins = proteins
        self.ligands = ligands
        self.protein_members = proteins[:]
        self.ligand_members = ligands[:]
        self.ligand_nodes = self.ligand_allnodes[:]
        self.protein_nodes = self.protein_allnodes[:]

        for i in self.fixed:
            if i in self.protein_nodes:
                self.protein_nodes.remove(i)
            if i in self.ligand_nodes:
                self.ligand_nodes.remove(i)

    # updates the proteins and ligands in order
    # maybe not necessary
    def update_proteins_ligands(self):
        n = self.get_start()
        all_ligands = self.ligand_members[:]
        all_ligands.extend(self.ligand_choices[:])
        all_proteins = self.protein_members[:]
        all_proteins.extend(self.protein_choices[:])
        self.ligands = []
        self.proteins = []
        if n in all_proteins:
            self.proteins.append(n)
        elif n in all_ligands:
            self.ligands.append(n)
        while len(self.graph.successors(n)) > 0:
            n = self.graph.successors(n)[0]
            if n in all_proteins: 
                self.proteins.append(n)
            elif n in all_ligands:
                self.ligands.append(n)

    # also, maybe not necessary
    def get_start(self):
        start = None
        for n in self.graph.nodes_iter():
            if len(self.graph.predecessors(n)) < 1:
                start = n
        if start is None:
            print self.graph.edges()
        return start

    def get_single_predecessor_node(self, node):
        pred = self.graph.predecessors(node)
        if len(pred) > 0:
            return pred[0]
        else:
            return None

    def get_single_successor_node(self, node):
        succ = self.graph.successors(node)
        if len(succ) > 0:
            return succ[0]
        else:
            return None

    def move_remove_single(self):
        all = self.protein_members[:]
        all.extend(self.ligand_members[:])
        node = self.choose_by_degree(1, all, 1)
        if node is not None:
            self.graph.remove_node(node)
            if node in self.ligand_members:
                self.ligand_members.remove(node)
                self.ligand_choices.append(node)
            elif node in self.protein_members:
                self.protein_members.remove(node)
                self.protein_choices.append(node)

    # Removes a pair of consecutive nodes, a protein node and a ligand node
    # Adds edge between predecessors to the successors of the removed nodes
    def move_remove_pair(self, pnode):
        lnodes = None
        if self.graph.in_edges(pnode) is not None and self.graph.out_edges(pnode) is not None:
            ppred = self.graph.predecessors(pnode)
            psucc = self.graph.successors(pnode)
            node_choices = []
            if len(ppred) > 0:
                node_choices.append(ppred)
            if len(psucc) > 0:
                node_choices.append(psucc)
            if len(node_choices) > 0:
                #lnodes = rand.choice(node_choices)
                lnode = rand.choice(node_choices)
                if len(ppred) > 0 and len(psucc) > 0:
                    for pp in ppred:
                        for ps in ppred:
                            self.graph.add_edge(pp, ps)
                self.graph.remove_node(pnode)
                self.protein_members.remove(pnode)
                self.protein_choices.append(pnode)
                
                if self.graph.in_edges(lnode) is not None and self.graph.out_edges(lnode) is not None:
                    lpred = self.graph.predecessors(lnode)
                    lsucc = self.graph.successors(lnode)
                    if len(lpred) > 0 and len(lsucc) > 0:
                        for lp in lpred:
                            for ls in lsucc:
                                self.graph.add_edge(lp, ls)
                self.graph.remove_node(lnode)
                self.ligand_members.remove(lnode)
                self.ligand_choices.append(lnode)

    # choose k nodes from a nodeset with a degree less than the limit
    def choose_by_degree_in(self, upperlim_deg, nodeset, k):
        choices = []
        for node in nodeset:
            degree = self.graph.in_degree(node)
            if degree < upperlim_deg:
                choices.append(node)
        if len(choices) >= k:
            return rand.sample(choices, k)
        else:
            return None

    def choose_by_degree_out(self, upperlim_deg, nodeset, k):
        choices = []
        for node in nodeset:
            degree = self.graph.out_degree(node)
            if degree < upperlim_deg:
                choices.append(node)
        if len(choices) >= k:
            return rand.sample(choices, k)
        else:
            return None

    def choose_by_degree(self, upperlim_deg, nodeset, k):
        choices = []
        for node in nodeset:
            degree = self.graph.degree(node)
            if degree < upperlim_deg:
                choices.append(node)
        if len(choices) >= k:
            return rand.sample(choices, k)
        else:
            return None

    # Adds a pair of consecutive nodes to the end of the linear graph
    def move_add_pair(self, pnode, lnode):
        last = self.choose_by_degree_in(2, self.ligand_members, 1)[0]
        if last is not None:
            self.graph.add_node(pnode)
            self.graph.add_node(lnode)
            self.graph.add_edge(last, pnode)
            self.graph.add_edge(pnode, lnode)

            # Updates lists to reflect the current state of the full system
            self.ligand_members.append(lnode)
            self.ligand_choices.remove(lnode)
            self.protein_members.append(pnode)
            self.protein_choices.remove(pnode)


    def move_add_edge(self):
        if rand.choice([0, 1]):
            pnode = self.choose_by_degree_out(2, self.protein_members, 1)
            lnode = self.choose_by_degree_in(2, self.ligand_members, 1)
            if pnode is not None and lnode is not None:
                self.graph.add_edge(pnode[0], lnode[0])
        else:
            pnode = self.choose_by_degree_in(2, self.protein_members, 1)
            lnode = self.choose_by_degree_out(2, self.ligand_members, 1)
            if pnode is not None and lnode is not None:
                self.graph.add_edge(lnode[0], pnode[0])

    def move_swap_pair(self, swapLigands=True):
        if swapLigands:
            node1, node2 = rand.sample(self.ligand_nodes, 2)
        else:
            node1, node2 = rand.sample(self.protein_nodes, 2)

        if (node1 not in self.fixed) and (node2 not in self.fixed):
            label1 = self.graph.node[node1]['label']
            label2 = self.graph.node[node2]['label']
            self.graph.node[node1]['label'] = label2
            self.graph.node[node2]['label'] = label1

    def move_swap_neighboring_pair(self, swapLigands=True):
        
        if swapLigands:
            node = rand.choice(self.ligand_nodes)
        else:
            node = rand.choice(self.protein_nodes)
        predprot = self.get_single_predecessor_node(node)
        succprot = self.get_single_successor_node(node)
        options = []
        if predprot is not None:
            predpred = self.get_single_predecessor_node(predprot)
            if predpred is not None:
                options.append(predpred)
        if succprot is not None:
            succsucc = self.get_single_successor_node(succprot)
            if succsucc is not None:
                options.append(succsucc)
        if self.fixed in options:
            options.remove(self.fixed)
        if len(options) > 0:
            swapnode = rand.choice(options)
            label1 = self.graph.node[node]['label']
            label2 = self.graph.node[swapnode]['label']
            self.graph.node[node]['label'] = label2
            self.graph.node[swapnode]['label'] = label1

    # assumes linear
    def move_lp_swap_pairs(self, pnode1, pnode2):
        lpred = self.graph.predecessors(pnode1)
        if not len(lpred) > 0:
            tempnode = pnode1
            pnode1 = pnode2
            pnode2 = tempnode
            lpred = self.graph.predecessors(pnode1)
        ppred = self.graph.predecessors(lpred[0])
        lsucc = self.graph.successors(pnode1)
        if len(ppred) > 0:
            self.graph.remove_edge(ppred[0], lpred[0])
        if len(lsucc) > 0:
            self.graph.remove_edge(pnode1, lsucc[0])
        if len(ppred) > 0 and len(lsucc) > 0:
            self.graph.add_edge(ppred[0], lsucc[0])
        lsucc2 = self.graph.successors(pnode2)
        if len(lsucc2) > 0:
            self.graph.remove_edge(pnode2, lsucc2[0])
            self.graph.add_edge(pnode1, lsucc2[0])
        self.graph.add_edge(pnode2, lpred[0])

    def move_swap_in(self, swapLigands=True):
        if swapLigands:
            select = True
            while select or onodelabel in self.keep:
                onode = rand.choice(self.ligand_nodes)
                onodelabel = self.graph.node[onode]['label']
                select = False
            newnodelabel = rand.choice(self.ligand_choices)
        else:
            onode = rand.choice(self.protein_nodes)
            onodelabel = self.graph.node[onode]['label']
            newnodelabel = rand.choice(self.protein_choices)
        self.graph.node[onode]['label'] = newnodelabel

        if onode in self.ligand_allnodes:
            self.ligand_choices.remove(newnodelabel)
            self.ligand_choices.append(onodelabel)
        elif onode in self.protein_allnodes:
            self.protein_choices.append(onodelabel)
            self.protein_choices.remove(newnodelabel)
            self.protein_members.append(newnodelabel)
            self.protein_members.remove(onodelabel)

    # lnode is the new node to swap in
    # tested, doesn't help at all
    def move_swap_in_similar(self, lnode):
        npdt = np.dtype([('ligand', np.str_, 35), ('tc', np.float32)])
        tcarray = []
        tcsum = 0.0
        for lig in self.ligand_members:
            tc = self.data.chsimdf.get_value(lig, lnode)
            tcarray.append((lig, tc))
        nptc = np.array(tcarray, dtype = npdt)

        cutoff = rand.uniform(0, 0.5)
        bool = nptc['tc'] > cutoff
        if bool.any():
            options = nptc['ligand'][bool]
        else:
            options = nptc['ligand']
        outnode = rand.choice(options) 

        pred = self.get_single_predecessor_node(outnode)
        succ = self.get_single_successor_node(outnode)

        self.graph.add_node(lnode)
        if pred is not None:
            self.graph.add_edge(pred, lnode)
            self.graph.remove_edge(pred, outnode)
        if succ is not None:
            self.graph.add_edge(lnode, succ)
            self.graph.remove_edge(outnode, succ)
        self.graph.remove_node(outnode)

        self.ligand_members.remove(outnode)
        self.ligand_choices.append(outnode)
        self.ligand_members.append(lnode)
        self.ligand_members.remove(lnode)

    def get_subgraph_successors(self, pnode, numnodes):
        currnode = pnode
        subgraph_nodes = [pnode]

        for n in range(numnodes - 1):
            if currnode is not None:
                nextnode = self.get_single_successor_node(currnode)
                if nextnode is not None:
                    currnode = self.get_single_successor_node(nextnode)
                    if currnode is self.fixed:
                        break
                    if currnode is not None:
                        subgraph_nodes.append(currnode)
            else:
                return False
        return subgraph_nodes

    def get_subgraph_predecessors(self, pnode, numnodes):
        currnode = pnode
        subgraph_nodes = [pnode]

        for n in range(numnodes - 1):
            if currnode is not None:
                nextnode = self.get_single_predecessor_node(currnode)
                if nextnode is not None:
                    currnode = self.get_single_successor_node(nextnode)
                    if currnode is self.fixed:
                        break
                    if currnode is not None:
                        subgraph_nodes.append(currnode)
            else:
                return False
        return subgraph_nodes.reverse()

    def move_subgraph(self):
        pnode = rand.choice(self.protein_nodes)
        numnodes = rand.choice(range(2, len(self.protein_nodes)/2 + 1))

        subgraph_nodes = self.get_subgraph_successors(pnode, numnodes)
        if not subgraph_nodes:
            subgraph_nodes = self.get_subgraph_predecessors(pnode, numnodes)
        if subgraph_nodes:
            pnode = subgraph_nodes[0]
            qnode = subgraph_nodes[-1]
            lnode = self.get_single_predecessor_node(pnode)

            rchoices = self.protein_nodes[:]
            for s in subgraph_nodes:
                if s in rchoices:
                    rchoices.remove(s)
            rnode = rand.choice(rchoices)

            # Remove edges
            if lnode is not None:
                lpred = self.get_single_predecessor_node(lnode)
                qsucc = self.get_single_successor_node(qnode)
                rsucc = self.get_single_successor_node(rnode)
                if lnode == rsucc:
                    rchoices.remove(rnode)
                    rnode = rand.choice(rchoices)
                    rsucc = self.get_single_successor_node(rnode)

                if lpred is not None:
                    self.graph.remove_edge(lpred, lnode)
                    if qsucc is not None:
                        self.graph.add_edge(lpred, qsucc)
                if qsucc is not None:
                    self.graph.remove_edge(qnode, qsucc)
                if rsucc is not None:
                    self.graph.remove_edge(rnode, rsucc)
                    self.graph.add_edge(qnode, rsucc)
                self.graph.add_edge(rnode, lnode)

    def string_repr(self):
        string = ''
        for e in self.graph.edges():
            string += "('%s', '%s'), " % (e[0], e[1])
        return string.rstrip(', ')

    def string_interactions(self):
        string = ''
        for e in self.graph.edges():
            string += "'%s', '%s', " % (e[0], e[1])
        return string.rstrip(', ')

    def linear_to_string(self):
        string = ''
        if self.transporter is not None and self.transporter in self.data.dockdf.index:
            string += '%s -> ' % self.transporter
        n = self.get_start()
        string += self.graph.node[n]['label']
        while len(self.graph.successors(n)) > 0: 
            n = self.graph.successors(n)[0]
            label = self.graph.node[n]['label']
            string += ' -> '
            string += label
        return string

    def add_constraint(self, fun):
        self.constraints.append(fun)

    def test_constraints(self):
        for fun in self.constraints:
            if not fun(self):
                return False
        return True


class proteinGraph(object):
    def __init__(self, proteins, datatables):
        graph = DiGraph()
        graph.add_nodes_from(proteins)
        self.data = datatables
        for i in range(len(proteins)-1):
            graph.add_edge(proteins[i], proteins[i+1])
        self.graph = graph
        self.proteins = proteins

    def compute_score_sea(self):
        seascore = 0
        count = 0
        for p in self.proteins:
            prot1 = self.graph.predecessors(p)
            for i in prot1:
                seascore += self.data.seadf.get_value(p, i)
                count += 1
            #prot2 = self.graph.successors(p)
            #for i in prot1:
            #    for j in prot2:
            #        seascore += self.data.seadf.get_value(i, j)
            #        count += 1.
        if count == 0:
            self.seascore=0
            return 0., count
        self.seascore = seascore/count
        return seascore/count #, count

    def linear_to_string(self):
        string = self.proteins[0]
        n = self.proteins[0]
        while len(self.graph.successors(n)) > 0:
            n = self.graph.successors(n)[0]
            string += ' -> '
            string += n
        return string

class proteinLigandPairs(object):
    # same order of lists [P1, P2], [L1, L2]
    # [(P1, L1), (P2, L2)]
    def __init__(self, proteins, ligands):
        pairs = []
        for i, p in enumerate(proteins):
            pairs.append([p, ligands[i]])
        self.pairs = pairs

    def linear_to_string(self):
        stringlist = []
        for i, p in enumerate(self.pairs):
            stringlist.append('%s -> %s' % (p[1], p[0]))
        return ', '.join(stringlist)

    def compute_score_dock(self, data):
        dockscore = 0
        count = 0
        for p in self.pairs:
            substrate = p[1]
            enzyme = p[0]
            dockscore += self.data.dockdf.get_value(p, i)
            count += 1.
        if count == 0:
            return np.inf
        self.dockscore = dockscore/count
        return dockscore/count

class ligandGraph(object):
    def __init__(self, ligands):
        graph = DiGraph()
        graph.add_nodes_from(ligands)
        for i in range(len(ligands)-1):
            graph.add_edge(ligands[i], ligands[i+1])
        self.graph = graph
        self.ligands = ligands

    # again, does not currently compute tc score
    def compute_score_tc(self, data):
        tcscore = 0
        count = 0
        for lig1 in self.ligands:
            lig2 = self.graph.successors(lig1)
            if len(lig2) == 1:
                tcscore += tc[0]
                count += 1.
        if count == 0:
            return 0
        return tcscore/count

    def linear_to_string(self):
        string = self.ligands[0]
        n = self.ligands[0]
        while len(self.graph.successors(n)) > 0:
            n = self.graph.successors(n)[0]
            string += ' -> '
            string += n
        return string

class modelGraph(object):
    def __init__(self, data, proteins=None, ligands=None, interactions=None):
        self.data = data
        self.proteins = proteins
        self.ligands = ligands
        self.interactions = interactions
        self.rxngraph = rxnGraph(data, proteins=proteins, ligands=ligands, interactions=interactions)
        self.restraints = []
        self.objscore = None

    def copy(self):
        newgraph = modelGraph(self.data)
        newgraph.rxngraph = self.rxngraph.copy()
        newgraph.restraints = self.restraints[:]
        newgraph.objscore = self.objscore
        return newgraph

    def __str__(self):
        return self.rxngraph.linear_to_string()

    def add_restraint(self, restraint):
        self.restraints.append(restraint)

    def add_restraints(self, restraints):
        self.restraints += restraints

    def add_transporter(self, transporter):
        self.rxngraph.add_transporter(transporter)

    def compute_scores(self):
        score = 0.0
        for res in self.restraints:
            score += res.compute_score(self.rxngraph)
        self.objscore = score
        return self.objscore

    def print_scores(self):
        #score = self.rxngraph.compute_scores()
        for res in self.restraints:
            res.compute_score(self.rxngraph)
            print res
        print 'OBJ: %.3f' % self.objscore

    def get_scores_strings(self):
        score = self.rxngraph.compute_scores()
        string = ''
        for res in self.restraints:
            string += '%s, ' % str(res)
        string += 'OBJ: %.3f'
        return string

    def construct_graph_from_path(self, proteins, ligands):
        self.proteins = proteins
        self.ligands = ligands
        self.rxngraph.construct_graph_from_path(proteins, ligands)

    def add_choices(self, protein_choices, ligand_choices):
        self.rxngraph.add_choices(protein_choices, ligand_choices)

    def move_subgraph(self):
        self.rxngraph.move_subgraph()

    def move_swap_in(self, swapLigands=True):
        self.rxngraph.move_swap_in(swapLigands=swapLigands)

    def move_swap_neighboring_pair(self, swapLigands=True):
        self.rxngraph.move_swap_neighboring_pair(swapLigands=swapLigands)

    def move_swap_pair(self, swapLigands=True):
        self.rxngraph.move_swap_pair(swapLigands=swapLigands)



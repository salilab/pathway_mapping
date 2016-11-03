#!/bin/python

import numpy as np
from itertools import permutations,combinations,product
import sys
import random as rand
from rxngraphs import *
from pathway_tables import pathTable
import networkx as nx
import time

# splits a list, xs, into two sublists, where the first sublist
# contains n elements
def rand_split(xs, n):
    rand.shuffle(xs)
    return xs[0:n], xs[n:]

class Paths(object):

    # Initialize a Paths object
    # Arguments:
    #     datatables - chem_datatables object
    #     tblhandle - writeable handle for h5tables for pathways
    #     seed - (default: random int) seed for random number generator
    def __init__(self, datatables, tblhandle, seed=None):
        if seed is None:
            seed = rand.randint(0, sys.maxint)
        rand.seed(seed)
        print 'Random seed %d' % seed
        self.data = datatables
        g = tblhandle.create_group('/', 'paths', 'Pathway Tables')
        self.tbl = tblhandle.create_table(g, 'pathTable', pathTable, 'pathway data')
        
    def protein_permutations(self, num):
        proteins = self.data.proteins
        protein_paths = permutations(proteins, num)
        return protein_paths

    def protein_combinations(self, num):
        proteins = self.data.proteins
        protein_paths = combinations(proteins, num)
        return protein_paths

    def ligand_permutations(self, num):
        ligands = self.data.ligands
        ligand_paths = permutations(ligands, num)
        return ligand_paths

    def store_graph(self, graph, step, T):
        row = self.tbl.row
        row['obj'] = graph.objscore
        row['step'] = step
        row['strrepr'] = str(graph)
        row.append()
    
    def sample_single_iteration(self, currstate, steps, mcparam, iteration=1, iterationlength=50000):
            
        def accept(newstate, currstate, step_norm):
            diff = currstate.objscore - newstate.objscore
            if diff <= 0:
                return True
            else:
                T = mcparam.get_temperature(step_norm)
                p = np.exp(-(diff)/T)
                return rand.uniform(0, 1) < p 
    
        def choose_move():
            move_choice = []
            if len(newstate.rxngraph.protein_members) > 3:
                for j in range(mcparam.m1):
                    move_choice.append(0)
            for j in range(mcparam.m2):
                move_choice.append(1)
            for j in range(mcparam.m3):
                move_choice.append(2)
            if len(newstate.rxngraph.protein_choices) > 0:
                for j in range(mcparam.m4):
                    move_choice.append(3)

            #move_choice.append(4)
            return rand.choice(move_choice) 
        step = 1
        acc_count = 0
    
        currscore = currstate.compute_scores()
        self.store_graph(currstate, step+(iterationlength*(iteration-1)), mcparam.get_temperature(0))
        currbest = currstate.objscore
        bestgraph = currstate.copy()
    
        while step < steps:
            newstate = currstate.copy()
            move_int = choose_move()
    
            if move_int == 0:
                if len(newstate.rxngraph.protein_nodes) > 3 and rand.choice([0, 1]):
                    newstate.move_subgraph()
                else:
                    newstate.move_swap_pair(swapLigands=False)
            elif move_int == 1:
                newstate.move_swap_neighboring_pair(swapLigands=True)
            elif move_int == 2:
                newstate.move_swap_in(swapLigands=True)
            elif move_int == 3:
                newstate.move_swap_in(swapLigands=False)
            # TESTING TO INCLUDE NOT NEIGHBORING
            # PROTEIN NODE SWAPS
            #elif move_int == 4: 
            #    newstate.move_swap_pair(swapLigands=False)
    
            newscore = newstate.compute_scores()
            step += 1
            
            if accept(newstate, currstate, float(step)/steps):
                currstate = newstate
                currscore = currstate.objscore
                acc_count += 1
            self.store_graph(currstate, step+(iterationlength*(iteration-1)), mcparam.get_temperature(float(step)/steps))
            if currstate.objscore > currbest:
                currbest = currstate.objscore
                bestgraph = currstate.copy()
    
        self.tbl.flush()
        return bestgraph, acc_count


    def sample_labels_monte_carlo(self, steps, protnum, lignum, restraints, random=False, fixed=[], startprots=[], 
                                  keep=[], transporter=None, numiterations=1, kparameter_set=0, parameter_set=0):
        #numiterations = 5000000/steps
        
        #numiterations = 5
        iterationlength = steps
        mcparam = montecarlo_parameters(parameter_set, kparameter_set)
        print mcparam
        print 'Steps: %d' % steps
        print 'Num Iterations: %d' % numiterations
        start_time = time.time()
        acc_count = 0
        step = 0
        proteins = self.data.proteins
        ligands = self.data.ligands
        protein_members, protein_choices = rand_split(proteins, protnum)
        ligand_members, ligand_choices = rand_split(ligands, lignum)

        currstate = modelGraph(self.data)
        currstate.add_restraints(restraints)
        currstate.rxngraph.keep = keep

        for i in range(len(keep)):
            ligand_members[i] = keep[i]
            if keep[i] in ligand_choices:
                ligand_choices.remove(keep[i])
        rand.shuffle(ligand_members) 

        if startprots:
            currstate.construct_graph_from_path(startprots, ligand_members)
        else:
            currstate.construct_graph_from_path(protein_members, ligand_members)


        currstate.add_choices(protein_choices=protein_choices, ligand_choices=ligand_choices)
        if transporter is not None:
            currstate.add_transporter(transporter)


        currscore = currstate.compute_scores()
        self.store_graph(currstate, step, mcparam.get_temperature(float(step)/steps))
        currbest = currstate.objscore
        bestpath = ''
        print 'Initial Path: %s, O: %.3f' % (str(currstate), currstate.objscore)


        bestgraph = currstate
        for j in range(1, 1+numiterations):
            bestgraph, single_acc_count = self.sample_single_iteration(bestgraph, steps, mcparam,
                                                     iteration=j, iterationlength=iterationlength)
            print 'Iteration %d: %s, O: %.3f' % (j, str(bestgraph), bestgraph.objscore)
            acc_count += single_acc_count
        print '---------------------------------------------'
        self.tbl.attrs.bestscore = bestgraph.objscore

        t = time.time() - start_time
        print 'Best Observed:', str(bestgraph)
        bestgraph.print_scores()

        print "Time: %.2f Seconds" % t
        print "Time: %.2f Minutes" % (t/60.)
        totalsteps = iterationlength*numiterations
        v = float(acc_count)/totalsteps
        print "Acceptance Ratio: %.5f" % v
        return float(acc_count)/totalsteps


    def random_graphs_stats(self, steps, protnum, lignum, transporter=None, keep=[]):
        g = modelGraph(self.data)
        g.add_restraints(restraints)
        prots = rand.sample(self.data.proteins, protnum)
        ligs = rand.sample(self.data.ligands, lignum)
        g.construct_graph_from_path(prots, ligs)

        g.rxngraph.keep = keep
        if transporter is not None:
            g.add_transporter(transporter)

        scores = []
        for i in range(steps):
            prots = rand.sample(self.data.proteins, protnum)
            ligs = rand.sample(self.data.ligands, lignum)
            for k in keep:
                if k in self.data.ligands and k not in ligs:
                    ligs.pop()
                    ligs.append(k)
                    rand.shuffle(ligs)
                if k in self.data.proteins and k not in prots:
                    prots.pop()
                    prots.append(k)
                    rand.shuffle(prots)
            g.proteins = prots
            g.ligands = ligs
            g.rxngraph.ligand_members = ligs
            g.rxngraph.protein_members = prots
            pdict = dict((n, prots[n]) for n in range(len(prots)))
            ldict = dict((n+len(prots), ligs[n]) for n in range(len(ligs)))
            nx.set_node_attributes(g.rxngraph.graph, 'label', pdict)
            nx.set_node_attributes(g.rxngraph.graph, 'label', ldict)
            g.compute_scores()
            scores.append(g.objscore)

    def random_graphs(self, steps, protnum, lignum, restraints, transporter=None, keep=[]):
        row = self.tbl.row
        g = modelGraph(self.data)
        g.add_restraints(restraints)
        prots = rand.sample(self.data.proteins, protnum)
        ligs = rand.sample(self.data.ligands, lignum)
        g.construct_graph_from_path(prots, ligs)


        g.rxngraph.keep = keep
        if transporter is not None:
            g.add_transporter(transporter)

        for i in range(steps):
            prots = rand.sample(self.data.proteins, protnum)
            ligs = rand.sample(self.data.ligands, lignum)
            for k in keep:
                if k in self.data.ligands and k not in ligs:
                    ligs.pop()
                    ligs.append(k)
                    rand.shuffle(ligs)
                if k in self.data.proteins and k not in prots:
                    prots.pop()
                    prots.append(k)
                    rand.shuffle(prots)

            g.proteins = prots
            g.ligands = ligs 
            g.rxngraph.ligand_members = ligs
            g.rxngraph.protein_members = prots
            pdict = dict((n, prots[n]) for n in range(len(prots)))
            ldict = dict((n+len(prots), ligs[n]) for n in range(len(ligs)))
            nx.set_node_attributes(g.rxngraph.graph, 'label', pdict)
            nx.set_node_attributes(g.rxngraph.graph, 'label', ldict)
            g.compute_scores()
            row['obj'] = g.objscore
            row.append()
        self.tbl.flush()

    def random_stats(self, steps, protnum, lignum, restraints, transporter=None, keep=[]):

        g = modelGraph(self.data)
        g.add_restraints(restraints)
        prots = rand.sample(self.data.proteins, protnum)
        ligs = rand.sample(self.data.ligands, lignum)
        g.construct_graph_from_path(prots, ligs)

        g.rxngraph.keep = keep
        if transporter is not None:
            g.add_transporter(transporter)
        scores = []
        for i in range(steps):
            prots = rand.sample(self.data.proteins, protnum)
            ligs = rand.sample(self.data.ligands, lignum)
            for k in keep:
                if k in self.data.ligands and k not in ligs:
                    ligs.pop()
                    ligs.append(k)
                    rand.shuffle(ligs)
                if k in self.data.proteins and k not in prots:
                    prots.pop()
                    prots.append(k)
                    rand.shuffle(prots)

            g.proteins = prots
            g.ligands = ligs
            g.rxngraph.ligand_members = ligs
            g.rxngraph.protein_members = prots
            pdict = dict((n, prots[n]) for n in range(len(prots)))
            ldict = dict((n+len(prots), ligs[n]) for n in range(len(ligs)))
            nx.set_node_attributes(g.rxngraph.graph, 'label', pdict)
            nx.set_node_attributes(g.rxngraph.graph, 'label', ldict)
            g.compute_scores()
            scores.append(g.objscore)

        sl = np.array(scores)
        return np.mean(sl), np.std(sl)

    def score_solutions(self, solutions, restraints, transporter=None):
        row = self.tbl.row

        best = None
        bestgraph = None

        pathlengths = np.dtype([('strrepr', np.str_, 200), ('length', np.int)])
        sollengths = []

        for sol in solutions:
            sollist = sol.split(' -> ')
            sollengths.append((sol, len(sollist)))

        solution_arr = np.array(sollengths, dtype=pathlengths)
        lengths = set(solution_arr['length'])
        #solution_arr.sort(order='length')

        for length in lengths:
            subset_solutions = solution_arr[np.where(solution_arr['length']==length)]['strrepr']

            path = subset_solutions[0].split(' -> ')
            ligs = [path[2*i] for i in range(len(path)/2 + 1)] 
            prots = [path[2*i+1] for i in range(len(path)/2)]
    
            g = modelGraph(self.data)
            g.add_restraints(restraints)
            if transporter is not None:
                g.add_transporter(transporter)
            g.construct_graph_from_path(prots, ligs)

            for solution in subset_solutions:
                path = solution.split(' -> ')
                ligs = [path[2*i] for i in range(len(path)/2 + 1)] 
                prots = [path[2*i+1] for i in range(len(path)/2)]

                g.proteins = prots
                g.ligands = ligs
                g.rxngraph.ligand_members = ligs
                g.rxngraph.protein_members = prots

                pdict = dict((n, prots[n]) for n in range(len(prots)))
                ldict = dict((n+len(prots), ligs[n]) for n in range(len(ligs)))

                nx.set_node_attributes(g.rxngraph.graph, 'label', pdict)
                nx.set_node_attributes(g.rxngraph.graph, 'label', ldict)

                g.compute_scores()
                row['strrepr'] = str(g)
                row['obj'] = g.objscore
                row.append()
                if g.objscore > best:
                    best = g.objscore
                    bestgraph = g.copy()

        self.tbl.attrs.bestscore = best
        print 'Best Observed:', str(bestgraph)
        bestgraph.print_scores()

        self.tbl.flush()

    def filter_solutions_by_stats(self, patharray, steps, restraints, num_stdevs=1.0, transporter=None, keep=[]): 
        pathlengths = np.dtype([('strrepr', np.str_, 200), ('obj', np.float64), ('length', np.int)])
        sollengths = []

        for sol in patharray:
            sollist = sol['strrepr'].split(' -> ')
            sollengths.append((sol['strrepr'], sol['obj'], len(sollist)))

        solution_arr = np.array(sollengths, dtype=pathlengths)
        lengths = set(solution_arr['length'])
        #solution_arr.sort(order='length')

        good_solutions_array = np.array([], dtype=pathlengths)
        for length in lengths:
            subset_solutions = solution_arr[np.where(solution_arr['length']==length)]

            path = subset_solutions['strrepr'][0].split(' -> ')
            ligs = [path[2*i] for i in range(len(path)/2 + 1)] 
            prots = [path[2*i+1] for i in range(len(path)/2)]

            mu, std = self.random_stats(steps, len(prots), len(ligs), restraints, transporter=None, keep=[])
            maxscore = np.max(subset_solutions['obj'])
            cu = maxscore - num_stdevs*std

            good_solutions = subset_solutions[subset_solutions['obj'] >= cu] 
            good_solutions_array = np.concatenate([good_solutions_array, good_solutions])

        # Change back to two data columns
        pathdt = np.dtype([('obj', np.float64), ('strrepr', np.str_, 200)])
        templist = [(x['obj'], x['strrepr']) for x in good_solutions_array]
        good_solutions_array = np.array(templist, dtype=pathdt)

        print '%d paths filtered down to %d within %.1f standard deviation(s)' % (len(patharray), len(good_solutions_array), num_stdevs)
        return good_solutions_array

    def sample_all_ligands(self, numligands):
        row = self.tbl.row
        g = rxnGraph(self.data)
        proteins = self.data.proteins
        ligands = self.data.ligands
        numprots = len(proteins)
        count = 0
        for L in permutations(ligands, numligands):
            dockscore = 0.0
            simscore = 0.0
            tfscore = 0.0

            for j, p in enumerate(proteins):
                if j >= len(L):
                    sub = L[j-1]
                else:
                    sub = L[j]
                dock = g.data.dockdf.get_value(p, sub)
                if dock > 0.0 and not math.isnan(dock):
                    dockscore += dock
                if p in g.data.tfluordf.index:
                    tfscore += g.data.tfluordf.get_value(p, sub)
                if (p in g.data.rsimpn.items): #and j < len(L)-1:
                    sub = L[j-1]
                    prod = L[j]
                    #simscore += g.data.chsimdf.get_value(sub, prod)
                    #simscore += g.data.rsimpn.get_value(p, sub, prod)
                    simscore += g.data.rsimpn.get_value(p, prod, sub)
            
            obj = 0.0
            obj += dockscore
            obj += simscore
            obj += tfscore

            if dockscore > 0. and simscore > 8.0:
                row['obj'] = obj
                strrepr = [None]*(len(proteins)+len(L))
                strrepr[::2] = proteins
                strrepr[1::2] = L
                print ' -> '.join(strrepr), simscore, tfscore, dockscore, obj
                row['strrepr'] = ' -> '.join(strrepr)
                row.append()
            count += 1
        self.tbl.flush()
        print count

    def sample_ligands(self, numligands, proteins, restraints, transporter=None):
        row = self.tbl.row
        g = modelGraph(self.data)
        g.add_restraints(restraints)


        ligands = self.data.ligands
        initialligs = rand.sample(self.data.ligands, numligands)

        g.construct_graph_from_path(proteins, initialligs)
        
        if transporter is not None:
            g.add_transporter(transporter)

        g.compute_scores()
        bestscore = g.objscore
        bestgraphstr = str(g)
        count = 0
        for L in permutations(ligands, numligands):
            g.ligands = L[:]
            g.rxngraph.ligand_members = L[:]
            ldict = dict((n+len(proteins), L[n]) for n in range(len(L)))
            nx.set_node_attributes(g.rxngraph.graph, 'label', ldict)
            g.compute_scores()

            row['obj'] = g.objscore
            row['strrepr'] = str(g)
            row.append()

            if bestscore < g.objscore:
                bestscore = g.objscore
                bestgraphstr = str(g)
            count += 1

        self.tbl.attrs.bestscore = bestscore
        self.tbl.flush()
        
        
        print 'Best Observed:', bestgraphstr
        print 'Best Score: %.3f' % bestscore
        print 'Number of Permutations: %d' % count

class montecarlo_parameters:

    def __init__(self, parameter_set=0, kparameter_set=0):
        if parameter_set:
            if parameter_set == 1:
                self.m1 = 1
                self.m2 = 1
                self.m3 = 1
                self.m4 = 1
            elif parameter_set == 2:
                self.m1 = 1
                self.m2 = 1
                self.m3 = 10
                self.m4 = 1

            elif parameter_set == 3:
                self.m1 = 1
                self.m2 = 1
                self.m3 = 50
                self.m4 = 1
        else:
            self.m1 = 15  # swap protein positions
            self.m2 = 1  # swap ligand positions
            self.m3 = 30
            self.m4 = 1  # swap new protein

        if kparameter_set:
            if kparameter_set == 1:
                self.k1 = 0.5
                self.k2 = 0.1
                self.k3 = 0.05
            elif kparameter_set == 2:
                self.k1 = 0.5
                self.k2 = 0.2
                self.k3 = 0.1
            elif kparameter_set == 3:
                self.k1  = 0.1
                self.k2 = 0.01
                self.k3 = 0.01
            elif kparameter_set == 4:
                self.k1 = 0.3
                self.k2 = 0.2
                self.k3 = 0.1
            elif kparameter_set == 5:
                self.k1 = 0.2
                self.k2 = 0.1
                self.k3 = 0.1
            elif kparameter_set == 6:
                self.k1 = 1.0
                self.k2 = 0.1
                self.k3 = 0.1
        else:
            self.k1 = .3
            self.k2 = .3
            self.k3 = 0.05
        # Good parameters for serine, kdo pathways k1=0.5, k2=0.1, k3=0.05
        # m1 = 15, m2 = 1, m3 = 50, m4 = 1

    def __repr__(self):
        string = 'Monte Carlo sampling'
        string += '\nSimulated annealing parameters: %.4f x %.4f ** step + %.4f' % (self.k1, self.k2, self.k3)
        string += '\nMove set parameters: %d %d %d %d' % (self.m1, self.m2, self.m3, self.m4)
        return string

    def get_temperature(self, step_norm, linear=False):
        if linear:
            temperature = self.k1 - step_norm*(self.k1 - self.k2)
        else:
            temperature = self.k1*self.k2**step_norm + self.k3
        return temperature



import pickle

class node():
    def __init__(self, variables, mappings):
        self.variables = variables
        self.mappings = mappings

    def __repr__(self):
        vkeys = list(self.variables.keys())
        vkeys.sort()
        string = ''
        for i in range(len(self.mappings)):
            oneset = [self.mappings[i][self.variables[v]] for v in vkeys]
            onestring = ' -> '.join(oneset)
            string += '%s\n' % onestring
        return string.strip()

    def write_pickle(self, picklefile):
        import pickle
        list_of_strings = str(self).split('\n')
        with open(picklefile, 'w') as handle:
            pickle.dump(list_of_strings, handle)


    def get_solutions(self):
        vkeys = list(self.variables.keys())
        vkeys.sort()
        strings = []
        for i in range(len(self.mappings)):
            oneset = [self.mappings[i][self.variables[v]] for v in vkeys]
            onestring = ' -> '.join(oneset)
            strings.append(onestring)
        return strings

def merge(node1, node2):
    vs = []
    nodevars1 = node1.variables.keys()
    nodevars2 = node2.variables.keys()
    varlist = list(set(nodevars1 + nodevars2))
    varinter = set(nodevars2).intersection(set(nodevars1))
    varunion = set(varlist)

    combinedvars = {}
    for i in range(len(varlist)):
        combinedvars[varlist[i]] = i
    storemappings = []

    for i in range(len(node1.mappings)):
        currmapping1 = node1.mappings[i]

        for j in range(len(node2.mappings)):
            currmapping2 = node2.mappings[j]
            valid = True

            for vi in varinter:
                vidx1 = node1.variables[vi]
                vidx2 = node2.variables[vi]

                if currmapping1[vidx1] != currmapping2[vidx2]:
                    valid = False
                    break

                if valid:
                    for k in nodevars1:
                        vidx1 = node1.variables[k]
                        if k not in varinter:
                            for m in nodevars2:
                                vidx2 = node2.variables[m]
                                if currmapping1[vidx1] == currmapping2[vidx2]:
                                    valid = False
                                    break

                if valid:
                    newmapping = [0 for v in combinedvars.keys()]
                    for v in nodevars1:
                        oldidx = node1.variables[v]
                        newidx = combinedvars[v]
                        newmapping[newidx] = currmapping1[oldidx]

                    for v in nodevars2:
                        if v not in nodevars1:
                            oldidx = node2.variables[v]
                            newidx = combinedvars[v]
                            newmapping[newidx] = currmapping2[oldidx]
                    storemappings.append(newmapping)

    return node(combinedvars, storemappings)

def make_merge(nodelist):
    if len(nodelist) <= 2:
        first = merge(nodelist[0], nodelist[1])
        print 'Number of solutions: %d' % len(first.mappings)
        return first
    else:
        last = nodelist.pop()
        intermed = merge(last, make_merge(nodelist))
        print 'Number of solutions: %d' % len(intermed.mappings)
        return intermed

def read_in_reactions(filename):
    rxns = []
    with open(filename, 'r') as handle:
        lines = handle.readlines()
        for line in lines:
            myline = line.split(',')
            e = myline[0].strip()
            s = myline[1].strip()
            p = myline[2].strip()
            rxns.append([s, e, p])
    return rxns

def read_in_reactions_df(datafile):
    import pandas as pd
    df = pd.read_hdf(datafile, 'reactions')
    rxns = [(df['substrate'][i], df['enzyme'][i], df['product'][i]) for i in df.index]
    return rxns

def enumerate_pathways(rxns, length, picklefile, minlength=None):
    if minlength is None:
        minlength = length

    variables = [{2*n:0, 2*n+1:1, 2*n+2:2} for n in range(length)]
    nodelist = [node(v, rxns) for v in variables]
    currsol = make_merge(nodelist[0:minlength])

    solutions = []
    nodelist_partial = nodelist[minlength::]
    solutions = currsol.get_solutions()

    for y in nodelist_partial:
        mergelist = [currsol, y]
        currsol = make_merge(mergelist)
        solutions.extend(currsol.get_solutions())

    if len(solutions) > 0:
        with open(picklefile, 'w') as handle:
            print 'Writing %d pathway solutions to %s' % (len(solutions), picklefile)
            pickle.dump(solutions, handle)


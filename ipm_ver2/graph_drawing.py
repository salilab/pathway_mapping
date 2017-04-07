import pygraphviz as pgv
import colorsys

# Returns a list of hexcolor strings of length N
def get_N_HexCol(N=5):
    if N < 15:
        HSV_tuples = get_small_HexCol(N)
    else:
        HSV_tuples = []
        denom = (N+2)/3
        for x in xrange(N):
            rescalex = x % denom
            if x < N/3:
                HSV_tuples.append((rescalex*1.0/denom, 0.6, 1.0))
            elif x < 2*N/3:
                HSV_tuples.append((rescalex*1.0/denom, 1.0, 1.0))
            else:
                HSV_tuples.append((rescalex*1.0/denom, 0.85, 0.6))

    hex_out = []
    for rgb in HSV_tuples:
        rgb = map(lambda x: int(x*255), colorsys.hsv_to_rgb(*rgb))
        hex_out.append("".join(map(lambda x: chr(x).encode('hex'), rgb)))
    return hex_out

# Returns a list of hexcolor strings of length N
# that is linearly distributed in Red values of RGB colors
def get_small_HexCol(N):
    HSV_tuples = [(x*1.0/N, 1.0, 0.9) for x in xrange(N)]
    return HSV_tuples

# Builds a pygraphviz graph
def build_pvgraph(nodes, edges, nodecolordict, edgecolordict, index_true,
                  truepathway, fontsize, circles=None, circlelabels=False,
                  diamondlabels=True, directed=False, labeldict={}):
    
    G = pgv.AGraph(strict=False, directed=directed, rankdir="LR", 
                   center=1, outputorder="edgesfirst", nodesep=0.05)
    if circles is None:
        circles = nodecolordict.keys()
    if labeldict == {}:
        for node in nodes:
            labeldict[node] = node
    for node in nodes:
        la = node
        if la in circles:
            color = nodecolordict[la]
            fillcolor = color
            penwidth = 1
            if node in truepathway:
                color = 'black'
                penwidth = 3
            if circlelabels:
                xlabel = labeldict[node]
            else:
                xlabel = ''
            G.add_node(node, width=0.3, height=0.3, shape="circle", label='',
                       xlabel=xlabel, fontname="helvetica", color=color, 
                       fillcolor=fillcolor, style="filled", fontsize=fontsize,
                       penwidth=penwidth)
        else:
            if la in nodecolordict.keys():
                color = nodecolordict[la]
                fillcolor = color
            else:
                fillcolor = 'white'
                color = 'black'
            fontcolor = 'black'
            if node in truepathway:
                color = 'black'
                fontcolor = 'white'
            if diamondlabels:
                label = labeldict[node]
            else:
                label = ''
            G.add_node(node, width=0.6, height=0.6, shape="diamond",
                       label=label, fontsize=fontsize, margin=0.001, 
                       fontname="helvetica", fontcolor=fontcolor, 
                       fillcolor=fillcolor, style="filled", color=color)
   
    if index_true is None:
        defaultcolor='black'
    else:
        defaultcolor='#008B8B'

    for edge in edges.keys():
        edgevals = edges[edge]
        colors = ''
        penwidth = 2
        for val in edgevals:
            if edgecolordict is not None:
                colors += '%s:' % edgecolordict[val]
            elif colors == '' and val != index_true:
                colors += '%s:' % defaultcolor

        for val in edgevals:
            if val == index_true:
                penwidth = 5
                colors += 'black:'
        colors.rstrip(':')
        G.add_edge(edge, color=colors, penwidth=penwidth)
    return G

def make_graph_img(paths, outfilename, nodecolordict=None, index_true=None, skeleton=False, fontsize=22):
    
    pathdict = {paths[i]: i for i in range(len(paths))}
    
    subtract = 0
    indexes = range(len(paths))
    if index_true is not None:
        indexes.remove(index_true)
        subtract = 1
    cscale = get_N_HexCol(len(paths)-subtract)
    edgecolordict = {i: '#' + cscale[i-subtract] for i in indexes}

    nodes = []
    posnodecolordict = {}
    edges = {}
    truep = []
    cnodes = []
    labeldict = {}
    
    if index_true is not None:
        edgecolordict[index_true] = 'black'
        truep = paths[index_true].split(' -> ')
        truep = ['%s_%d' % (truep[i], i) for i in range(len(truep))]
    
    for path in pathdict.keys():
        positions = path.split(' -> ')
        currnodes = []
        
        for pos in range(len(positions)):
            name = '%s_%d' % (positions[pos], pos)
            currnodes.append(name)
            labeldict[name] = positions[pos]
            if pos % 2 == 0:
                posnodecolordict[positions[pos]] = 0
                cnodes.append(name)
        nodes.extend(currnodes)
        
        for n in range(len(currnodes)-1):
            edge = (currnodes[n], currnodes[n+1])
            if edge not in edges.keys():
                edges[edge] = [pathdict[path]]
            else:
                edges[edge].append(pathdict[path])
    
    if nodecolordict is None:
        nodecolordict = posnodecolordict
        cscale = get_N_HexCol(len(nodecolordict))
        for i in range(len(nodecolordict)):
            k = nodecolordict.keys()[i]
            nodecolordict[k] = '#' + cscale[i]
        nodes = list(set(nodes))

    #allnodes_colordict = {}
    for j in cnodes:
        key = labeldict[j]
        ncolor = nodecolordict[key]
        #allnodes_colordict[j] = ncolor
        nodecolordict[j] = ncolor
    #nodecolordict = allnodes_colordict
    
    G = build_pvgraph(nodes, edges, nodecolordict, edgecolordict, index_true,
                      truep, fontsize, circlelabels=False, labeldict=labeldict)

    if skeleton:
        length = len(paths[0].split(' -> '))/2 + 1 
        G = add_skeleton(G, length)
        
    G.layout(prog="dot")
    G.draw(outfilename)
    return nodecolordict


def make_graph_img_net(paths, outfilename, nodecolordict=None, index_true=None, fontsize=22,
                       circlelabels=False, diamondlabels=True, directed=False):
    
    pathdict = {paths[i]: i for i in range(len(paths))}
    
    indexes = range(len(paths))
    subtract = 0
    if index_true is not None:
        indexes.remove(index_true)
        subtract = 1
    cscale = get_N_HexCol(len(paths)-subtract)
    edgecolordict = {i: '#' + cscale[i-subtract] for i in indexes}
    
    nodes = []
    circles = []
    edges = {}
    truep = []
    
    if index_true is not None:
        edgecolordict[index_true] = 'black'
        truep = paths[index_true].split(' -> ')
        
    for path in pathdict.keys():
        currnodes = path.split(' -> ')
        nodes.extend(currnodes)
        for pos in range(len(currnodes)):
            if pos % 2 == 0:
                circles.append(currnodes[pos])
        
        for n in range(len(currnodes)-1):
            edge = (currnodes[n], currnodes[n+1])
            if edge not in edges.keys():
                edges[edge] = [pathdict[path]]
            else:
                edges[edge].append(pathdict[path])

    nodes = list(set(nodes))
    circles = list(set(circles))

    if nodecolordict is None:
        cscale = get_N_HexCol(len(circles))
        new_nodecolordict = {}
        for i in range(len(circles)):
            new_nodecolordict[circles[i]] = '#' + cscale[i]
        nodecolordict = new_nodecolordict

    G = build_pvgraph(nodes, edges, nodecolordict, edgecolordict, index_true,
                      truep, fontsize, circles=circles, circlelabels=circlelabels,
                      diamondlabels=diamondlabels, directed=directed, labeldict={})

    G.layout(prog="dot")
    G.draw(outfilename)
    return nodecolordict


# Adds an outline, or "skeleton", of linear pathway of
# a given length to the graph
# Input: 
#  G - pygraphviz graph object
#  length - number of nodes in skeleton pathway
# Output:
#  Edited pygraphviz graph object with the skeleton
#  pathway added
def add_skeleton(G, length):
    for i in range(length):
        G.add_node(2*i, width=0.3, height=0.3, shape="circle", label="")
        if i > 0:
            G.add_edge(2*i-1, 2*i)
        if i < length - 1:
            G.add_node(2*i+1, width=0.6, height=0.6, shape="diamond", label="")
            G.add_edge((2*i, 2*i+1))
    return G


def make_graph_img_net_single(paths, outfilename, nodecolordict=None, index_true=None, fontsize=22,
                       circlelabels=False, diamondlabels=True, directed=False,
                       dotfile=None, true_edges_subenz=[], true_edges_enzprod=[]):
       
    cscale = get_N_HexCol(len(paths))
       
    nodes = []
    edges = {}
    
    circles = []

    idxpathdict = {i:paths[i] for i in range(len(paths))}
    for pathkey in idxpathdict.keys():
        currpath = idxpathdict[pathkey]
        currnodes = currpath.split(' -> ')
        nodes.extend(currnodes)
        for pos in range(len(currnodes)):
            if pos % 2 == 0:
                circles.append(currnodes[pos])
        for n in range(len(currnodes)-1):
            edge = (currnodes[n], currnodes[n+1])
            if edge not in edges.keys():
                edges[edge] = [pathkey]
            else:
                edges[edge].append(pathkey)

    if len(true_edges_subenz) > 0 and len(true_edges_enzprod) > 0:
        index_true = -1
        for edge_subenz in true_edges_subenz:
            nodes.append(edge_subenz[0])
            nodes.append(edge_subenz[1])
            circles.append(edge_subenz[0])
            if edge_subenz not in edges.keys():
                edges[edge_subenz] = [-1]
            else:
                edges[edge_subenz].append(-1)
        for edge_enzprod in true_edges_enzprod:
            nodes.append(edge_enzprod[0])
            nodes.append(edge_enzprod[1])
            circles.append(edge_enzprod[1])
            if edge_enzprod not in edges.keys():
                edges[edge_enzprod] = [-1]
            else:
                edges[edge_enzprod].append(-1)

    nodes = list(set(nodes))
    circles = list(set(circles))

    if nodecolordict is None:
        cscale = get_N_HexCol(len(circles))
        new_nodecolordict = {}
        for i in range(len(circles)):
            new_nodecolordict[circles[i]] = '#' + cscale[i]
        nodecolordict = new_nodecolordict

    G = build_pvgraph(nodes, edges, nodecolordict, None, index_true,
                      [], fontsize, circles=circles, circlelabels=circlelabels, 
                      diamondlabels=diamondlabels, directed=directed, labeldict={})
       
    G.layout(prog="dot")
    G.draw(outfilename)
    if dotfile is not None:
        G.write(dotfile)
    return nodecolordict

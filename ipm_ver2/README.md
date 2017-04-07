# Integrative Pathway Mapping

## Dependencies:
- pandas
- numpy
- RDKit
- PyTables
- NetworkX
- alignment
- scipy
- pygraphviz
- graphviz

## Python modules:

#### reaction_calculations.py
Module for processing molecules and performing calculations for restraints/scores.
Dependencies:
- pandas
- numpy
- RDKit
- PyTables

#### rxngraphs.py
Module containing classes for model representations of pathway graphs.
Dependencies:
- NetworkX
- pandas
- numpy

#### sample_graph.py
Module for pathway sampling (Monte Carlo methods) and scoring pathways
Dependencies:
- networkx
- numpy
- pathway_tables
    -PyTables

#### pathway_analysis.py
Module for clustering and filtering solutions
Dependencies:
- numpy

For clustering:
- alignment
- scipy

For loading pathways as a numpy array from a hdf5 file:
- PyTables

#### backtracking.py
Module for enumerating pathways

#### graph_drawing.py
Functions for network visualization; Abstract visualizations for networks and pathways using standard graph drawing libraries
Dependencies:
- pygraphviz
- graphviz

#### make_json.py
Functions for creating json files for pathway graphs that can be imported into Cytoscape NetIMP
Dependencies:
- rxngraphs
- pathway_restraints


#!/bin/python

#from tables import *
from tables import IsDescription
from tables import StringCol
from tables import Float32Col
from tables import Int32Col

class pathTable(IsDescription):
    strrepr = StringCol(200)
    obj = Float32Col()
    step = Int32Col()                    
    #dock = Float32Col()
    #tc = Float32Col()
    #mcorr = Float32Col() # metabolite-metabolite correlations
    #terms = Int32Col() # number of terms used in the objective function
    #                   # would like to track to see if less terms lead to lower
    #                   # objective fxns just because there are less to compute
    #rxn = Float32Col()
    #obj2 = Float32Col()

    #tf = Float32Col()

    #temperature = Float32Col()

class scoreTable(IsDescription):
    obj = Float32Col()
    step = Int32Col()

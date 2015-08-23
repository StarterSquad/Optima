# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 18:36:39 2015

@author: cliffk
"""

import add_optima_paths # analysis:ignore
from region import Region
import os
from numpy import sort
from sim import Sim


# List all the json files in the regions sub-directory.
templist = sort([x for x in os.listdir('./regions/') if x.endswith('.json')])
r1 = Region.load('./regions/' + templist[0])

s = Sim('test-sim',r1)
s.initialise()
S = s.run()

from makeccocs import plotallcurves
D = dict()
D['G'] = r1.metadata
D['data'] = r1.data
D['opt'] = r1.options
D['programs'] = r1.metadata['programs']
D['P'] = s.parsdata
D['S'] = S
D['M'] = s.parsmodel
for program in D['programs']:
    progname = program['name']
    if len(program['effects']):
        print(progname)
        plotdata_cco, plotdata_co, plotdata_cc, effects, D = plotallcurves(D, unicode(progname))
    else:
        print(progname+' not used')
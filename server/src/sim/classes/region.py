# -*- coding: utf-8 -*-
"""
Created on Fri May 29 23:16:12 2015

@author: David Kedziora
"""

import defaults
from simbox import SimBox


class Region:
    def __init__(self, name,populations,programs,datastart,dataend):
        self.D = dict()                 # Data structure for saving everything. Will hopefully be broken down eventually.
              
        self.metadata = {}            # Loosely analogous to D['G']
        self.metadata['datastart'] = datastart
        self.metadata['dataend'] = dataend
        self.metadata['populations'] = populations
        self.metadata['programs'] = programs
        self.metadata['name'] = name

        self.options = None             # This used to be D['opt']. Is it constant? Or should it be tagged 'default'?
        
        self.simboxlist = []            # Container for simbox objects (e.g. optimisations, grouped scenarios, etc.)
        
    def createsimbox(self, simboxname):
        self.simboxlist.append(SimBox(simboxname))
        
        # Makes sure that a new simbox has at least one sim object, ready to run. It is passed regional data, metadata and options.
        self.simboxlist[-1].createsim(simboxname+'-initial', self.data, self.metadata, self.options)
        
    # Runs through every simulation in simbox (if not processed) and processes them.
    def runsimbox(self, simbox):
        simbox.runallsims(self.data, self.metadata, self.options, self.programs, forcerun = False)
    
    # Runs through every simulation in simbox (if not processed) and optimises them.
    # Currently uses default settings.
    def optsimbox(self, simbox):
        simbox.optallsims(self.data, self.metadata, self.options, self.programs, forcerun = True)
        
    # Runs through every simulation in simbox (if processed) and plots them.
    def plotsimbox(self, simbox):
        simbox.plotallsims()
        
    def printdata(self):
        print(self.data)
    
    def printmetadata(self):
        print(self.metadata)
        
    def printoptions(self):
        print(self.options)
        
    def printprograms(self):
        print(self.programs)
        
    def printsimboxlist(self, assubset = False):
        # Prints with nice arrow formats if assubset is true. Otherwise numbers the list.        
        if assubset:
            if len(self.simboxlist) > 0:
                for simbox in self.simboxlist:
                    print(' --> %s' % simbox.getsimboxname())
                    simbox.printsimlist(assubsubset = True)
        else:
            if len(self.simboxlist) == 0:
                print('No simulations are currently associated with region %s.' % self.getregionname())
            else:
                print('Collections of simulations associated with this project...')
                fid = 0
                for simbox in self.simboxlist:
                    fid += 1
                    print('%i: %s' % (fid, simbox.getsimboxname()))
                    simbox.printsimlist(assubsubset = False)
                
    def setdata(self, data):
        self.data = data
        
    def getdata(self):
        return self.data
        
    def setmetadata(self, metadata):
        self.metadata = metadata
        
    def getmetadata(self):
        return self.metadata
        
    def setoptions(self, options):
        self.options = options
        
    def getoptions(self):
        return self.options
    
    def setprograms(self, programs):
        self.programs = programs
        
    def getprograms(self):
        return self.programs
        
    def setregionname(self, regionname):
        self.metadata['name'] = regionname
        
    def getregionname(self):
        return self.metadata['name']
        
    ### Refers to legacy D.
        
    def loadDfrom(self, path):
        from dataio import loaddata
        tempD = loaddata(path)
        self.setD(tempD)                # It would be great to get rid of setD one day. But only when data is fully decomposed.
        self.setdata(tempD['data'])
        self.setmetadata(tempD['G'])
        self.setoptions(tempD['opt'])
        self.setprograms(tempD['programs'])
        
    def setD(self, D):
        self.D = D
        
    def getD(self):
        return self.D
        
    def makeproject(self, projectname='example', pops = defaults.default_pops, progs = defaults.default_progs, datastart = defaults.default_datastart, \
        dataend = defaults.default_dataend, nsims = defaults.default_nsims, verbose=2, savetofile = True, domakeworkbook=True):
        """
        Initializes the empty project. Only the "Global" and "Fitted" parameters are added on this step.
        The rest of the parameters are calculated after the model is updated with the data from the workbook.
        
        Version: 2015jan27 by cliffk
        """
        
        from dataio import savedata, projectpath
        from printv import printv
        from numpy import arange
        from copy import deepcopy
    
        printv('Making project...', 1, verbose)
    
        self.D['plot'] = dict() # Initialize plotting data
        
        # Initialize options
        from setoptions import setoptions
        self.D['opt'] = setoptions(nsims=nsims)
        
        # Set up "G" -- general parameters structure
        self.D['G'] = dict()
        self.D['G']['version'] = defaults.current_version # so that we know the version of new project with regard to data structure
        self.D['G']['projectname'] = projectname  
        self.D['G']['projectfilename'] = projectpath(projectname+'.prj')
        self.D['G']['workbookname'] = self.D['G']['projectname'] + '.xlsx'
        self.D['G']['npops'] = len(pops)
        self.D['G']['nprogs'] = len(progs)
        self.D['G']['datastart'] = datastart
        self.D['G']['dataend'] = dataend
        self.D['G']['datayears'] = arange(self.D['G']['datastart'], self.D['G']['dataend']+1)
        self.D['G']['inputprograms'] = deepcopy(progs) # remember input programs with their possible deviations from standard parameter set (if entered from GUI). 
        # Hate duplicating the data, but can't think of a cleaner way of export/import.
        self.D['G']['inputpopulations'] = deepcopy(pops) # should be there as well, otherwise we cannot export project without data
        # Health states
        self.D['G']['healthstates'] = ['acute', 'gt500', 'gt350', 'gt200', 'gt50', 'aids']
        self.D['G']['ncd4'] = len(self.D['G']['healthstates'])
        self.D['G']['nstates'] = 1+self.D['G']['ncd4']*5 # Five are undiagnosed, diagnosed, 1st line, failure, 2nd line, plus susceptible
        self.D['G']['sus']  = arange(0,1)
        self.D['G']['undx'] = arange(0*self.D['G']['ncd4']+1, 1*self.D['G']['ncd4']+1)
        self.D['G']['dx']   = arange(1*self.D['G']['ncd4']+1, 2*self.D['G']['ncd4']+1)
        self.D['G']['tx1']  = arange(2*self.D['G']['ncd4']+1, 3*self.D['G']['ncd4']+1)
        self.D['G']['fail'] = arange(3*self.D['G']['ncd4']+1, 4*self.D['G']['ncd4']+1)
        self.D['G']['tx2']  = arange(4*self.D['G']['ncd4']+1, 5*self.D['G']['ncd4']+1)
        for i,h in enumerate(self.D['G']['healthstates']): self.D['G'][h] = [self.D['G'][state][i] for state in ['undx', 'dx', 'tx1', 'fail', 'tx2']]
        
        if savetofile: #False if we are using database
            savedata(self.D['G']['projectfilename'], self.D, verbose=verbose) # Create project -- #TODO: check if an existing project exists and don't overwrite it
        
        # Make an Excel template and then prompt the user to save it
        if projectname == 'example': # Don't make a new workbook, but just use the existing one, if the project name is "example"
            print('WARNING, Project name set to "example", not creating a new workbook!')
        else: # Make a new workbook
            if domakeworkbook:
                self.makeworkbook(self.D['G']['workbookname'], pops, progs, datastart, dataend, verbose=verbose)
        
        printv('  ...done making project.', 2, verbose)
    
    
    def makeworkbook(self, name, pops, progs, datastart=defaults.default_datastart, dataend=defaults.default_dataend, verbose=2):
        """ Generate the Optima workbook -- the hard work is done by makeworkbook.py """
        from printv import printv
        from dataio import templatepath
        from makeworkbook import OptimaWorkbook
    
        printv("""Generating workbook with parameters:
                 name = %s, pops = %s, progs = %s, datastart = %s, dataend = %s""" \
                 % (name, pops, progs, datastart, dataend), 1, verbose)
        path = templatepath(name)
        book = OptimaWorkbook(name, pops, progs, datastart, dataend)
        book.create(path)
        
        printv('  ...done making workbook %s.' % path, 2, verbose)
        return path
"""
This module defines the Program and Programset classes, which are
used to define a single program/modality (e.g., FSW programs) and a
set of programs, respectively.

Version: 2017feb15
"""

from optima import OptimaException, Link, odict, objrepr, promotetoarray, promotetolist, defaultrepr, checktype
from numpy.random import uniform, seed, get_state


class Programset(object):
    """
    Object to store all programs. Coverage-outcome data and functions belong to the program set, 
    while cost-coverage data/functions belong to the individual programs.
    """

    def __init__(self, name='default', parsetname=-1, project=None, programs=None):
        """ Initialize """
        if project is None:
            errormsg = 'To create a program set, you must supply a project as an argument'
            raise OptimaException(errormsg)
        
        self.name = name
        self.programs = odict()
        self.covout = odict()
        self.parsetname = parsetname # Store the parset name
        self.projectref = Link(project) # Store pointer for the project, if available
#        self.denominators = self.setdenominators() # Calculate the denominators for different coverage values
        if programs is not None: self.addprograms(programs)


    def __repr__(self):
        """ Print out useful information"""
        output = objrepr(self)
        output += '    Program set name: %s\n'    % self.name
        output += '            Programs: %s\n'    % self.programs.keys()
        output += '============================================================\n'
        return output
    
    
    def addprograms(self, progs=None, replace=False):
        ''' Add a list of programs '''
        if progs is not None:
            progs = promotetolist(progs)
        else:
            errormsg = 'Programs to add should not be None'
            raise OptimaException(errormsg)
        if replace:
            self.programs = odict()
        for prog in progs:
            if isinstance(prog, dict):
                prog = Program(**prog)
            if type(prog)!=Program:
                errormsg = 'Programs to add must be either dicts or program objects, not %s' % type(prog)
                raise OptimaException(errormsg)
            self.programs[prog.short] = prog
        return None


    def rmprograms(self, progs=None, die=True):
        ''' Remove one or more programs '''
        if progs is None:
            self.programs = odict() # Remove them all
        progs = promotetolist(progs)
        for prog in progs:
            try:
                self.programs.pop[prog]
            except:
                errormsg = 'Could not remove program named %s' % prog
                if die: raise OptimaException(errormsg)
                else: print(errormsg)
            for co in self.covout.values(): # Remove from coverage-outcome functions too
                co.progs.pop(prog, None)
        
        
        return None
    
    
    def addcovout(self):
        ''' add coverage-outcome parameter '''
        pass

    def defaultbudget(self, total=True):
        ''' Get default budget -- either per program or total '''
        if total:      budget = 0
        else: budget = odict()
        for prog in self.programs.values():
            if total: budget += prog.spend
            else:     budget[prog.short] = prog.spend
        return budget

    def coverage2budget(self):
        ''' get budget from coverage '''
        pass

    def budget2coverage(self):
        ''' get coverage from budget '''
        pass

    def reconcile(self):
        ''' reconcile with parset '''
        pass

    def compareoutcomes(self):
        ''' compare textually '''
        pass
    
    def check(self):
        ''' checks that all costcov and covout data is entered '''
        pass



class Program(object):
    ''' Defines a single program. '''
    
    def __init__(self, short=None, name=None, category=None, spend=None, basespend=None, coverage=None, unitcost=None, saturation=None, targetpops=None, targetpars=None):
        self.short      = None # short name
        self.name       = None # full name
        self.category   = None # spending category
        self.spend      = None # latest or estimated expenditure
        self.basespend  = None # non-optimizable spending
        self.coverage   = None # latest or estimated coverage (? -- not used)
        self.unitcost   = None # dataframe of [t, best, low, high]
        self.saturation = None # saturation coverage value
        self.targetpops = None # key(s) for targeted populations
        self.targetpars = None # which parameters are targeted
        self.update(short=short, name=name, category=category, spend=spend, basespend=basespend, coverage=coverage, unitcost=unitcost, saturation=saturation, targetpops=targetpops, targetpars=targetpars)
#        self.denominators = self.setdenominators()
    
    def __repr__(self):
        output = defaultrepr(self)
        return output
    
    
    def update(self, short=None, name=None, category=None, spend=None, basespend=None, coverage=None, unitcost=None, saturation=None, targetpops=None, targetpars=None):
        ''' Add data to a program '''
        if short      is not None: self.short      = checktype(short,    'string') # short name
        if name       is not None: self.name       = checktype(name,     'string') # full name
        if category   is not None: self.category   = checktype(category, 'string') # spending category
        if spend      is not None: self.spend      = None # latest or estimated expenditure
        if basespend  is not None: self.basespend  = None # non-optimizable spending
        if coverage   is not None: self.coverage   = None # latest or estimated coverage (? -- not used)
        if unitcost   is not None: self.unitcost   = None # dataframe of [t, best, low, high]
        if saturation is not None: self.saturation = None # saturation coverage value
        if targetpops is not None: self.targetpops = promotetolist(targetpops, 'string') # key(s) for targeted populations
        if targetpars is not None: self.targetpars = None # which parameters are targeted



class Covout(object):
    ''' A coverage-outcome object -- cost-outcome objects are incorporated in programs '''
    
    def __init__(self, par=None, pop=None, lowerlim=None, upperlim=None, progs=None):
        self.par = par
        self.pop = pop
        self.lowerlim = Val(lowerlim)
        self.upperlim = Val(upperlim)
        self.progs = odict()
        if isinstance(progs, dict):
            for key,val in progs.items():
                self.progs[key] = Val(val)
        return None
    
    def add(self, prog=None, val=None):
        ''' 
        Accepts either
        self.add({'FSW':[0.3,0.1,0.4]})
        or
        self.add('FSW', 0.3)
        '''
        if isinstance(prog, dict):
            if isinstance(prog, dict):
                for key,val in prog.items():
                    self.progs[key] = Val(val)
        elif isinstance(prog, basestring) and val is not None:
            self.progs[prog] = Val(val)
        else:
            errormsg = 'Could not understand prog=%s and val=%s' % (prog, val)
            raise OptimaException(errormsg)
        return None
            




class Val(object):
    ''' A single value for a coverage-outcome curve -- seems insanely complicated, I know! '''
    
    def __init__(self, best=None, low=None, high=None, dist=None):
        ''' Allow the object to be initialized, but keep the same infrastructure for updating '''
        self.best = None
        self.low = None
        self.high = None
        self.dist = None
        self.update(best=best, low=low, high=high, dist=dist)
        return None
    
    
    def __call__(self, *args, **kwargs):
        ''' Convenience function for both update and get '''
        
        # If it's None or if the key is a string (e.g. 'best'), get the values:
        if len(args)+len(kwargs)==0 or 'what' in kwargs or (len(args) and type(args[0])==str):
            return self.get(*args, **kwargs)
        else: # Otherwise, try to set the values
            self.update(*args, **kwargs)
    
    
    def update(self, best=None, low=None, high=None, dist=None):
        ''' Actually set the values -- very convoluted, but should be flexible and work :)'''
        
        # Reset these values if already supplied
        if best is None and self.best is not None: best = self.best
        if low  is None and self.low  is not None: low  = self.low 
        if high is None and self.high is not None: high = self.high 
        if dist is None and self.dist is not None: dist = self.dist
        
        # Handle values
        if best is None: # Best is not supplied, so use high and low, e.g. Val(low=0.2, high=0.4)
            if low is None or high is None:
                errormsg = 'If not supplying a best value, you must supply both low and high values'
                raise OptimaException(errormsg)
            else:
                best = (low+high)/2. # Take the average
        elif isinstance(best, dict):
            self.update(**best) # Assume it's a dict of args, e.g. Val({'best':0.3, 'low':0.2, 'high':0.4})
        else: # Best is supplied
            best = promotetoarray(best)
            if len(best)==1: # Only a single value supplied, e.g. Val(0.3)
                best = best[0] # Convert back to number
                if low is None: low = best # If these are missing, just replace them with best
                if high is None: high = best
            elif len(best)==2: # If length 2, assume high-low supplied, e.g. Val([0.2, 0.4])
                if low is not None and high is not None:
                    errormsg = 'If first argument has length 2, you cannot supply high and low values'
                    raise OptimaException(errormsg)
                low = best[0]
                high = best[1]
                best = (low+high)/2.
            elif len(best)==3: # Assume it's called like Val([0.3, 0.2, 0.4])
                low, best, high = sorted(best) # Allows values to be provided in any order
            else:
                errormsg = 'Could not understand input of best=%s, low=%s, high=%s' % (best, low, high)
                raise OptimaException(errormsg)
        
        # Handle distributions
        validdists = ['uniform']
        if dist is None: dist = 'uniform'
        if dist not in validdists:
            errormsg = 'Distribution "%s" not valid; choices are: %s' % (dist, validdists)
            raise OptimaException(errormsg) 
        
        # Store values
        self.best = best
        self.low  = low
        self.high = high
        self.dist = dist
        
        return None
    
    
    def get(self, what=None, n=1):
        ''' Get the value from this distribution '''
        if what is None or what is 'best': # Haha this is funny but works
            val = self.best
        elif what is 'low':
            val = self.low
        elif what is 'high':
            val = self.high
        elif what in ['rand','random']:
            seed(get_state()[1][0]) # Reset the random seed, if specified -- WARNING, will affect rest of Optima!
            val = uniform(low=self.low, high=self.high, size=n)
        return val
    
    

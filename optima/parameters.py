"""
This module defines the Constant, Metapar, Timepar, and Popsizepar classes, which are 
used to define a single parameter (e.g., hivtest) and the full set of
parameters, the Parameterset class.

Version: 2.1 (2017apr04)
"""

from numpy import array, nan, isnan, isfinite, zeros, argmax, mean, log, polyfit, exp, maximum, minimum, Inf, linspace, median, shape
from numpy.random import uniform, normal, seed
from optima import OptimaException, Link, odict, dataframe, printv, sanitize, uuid, today, getdate, makefilepath, smoothinterp, dcp, defaultrepr, isnumber, findinds, getvaliddata, promotetoarray, promotetolist, inclusiverange # Utilities 
from optima import Settings, getresults, convertlimits, gettvecdt, loadpartable, loadtranstable # Heftier functions
import optima as op

defaultsmoothness = 1.0 # The number of years of smoothing to do by default
generalkeys = ['male', 'female', 'popkeys', 'injects', 'fromto', 'transmatrix'] # General parameter keys that are just copied
staticmatrixkeys = ['birthtransit','agetransit','risktransit'] # Static keys that are also copied, but differently :)


#################################################################################################################################
### Define the parameter set
#################################################################################################################################



class Parameterset(object):
    ''' Class to hold all parameters and information on how they were generated, and perform operations on them'''
    
    def __init__(self, name='default', project=None, progsetname=None, budget=None, start=None, end=None):
        self.name = name # Name of the parameter set, e.g. 'default'
        self.uid = uuid() # ID
        self.projectref = Link(project) # Store pointer for the project, if available
        self.created = today() # Date created
        self.modified = today() # Date modified
        self.pars = None
        self.popkeys = [] # List of populations
        self.posterior = odict() # an odict, comparable to pars, for storing posterior values of m -- WARNING, not used yet
        self.resultsref = None # Store pointer to results
        self.progsetname = progsetname # Store the name of the progset that generated the parset, if any
        self.budget = budget # Store the budget that generated the parset, if any
        self.start = start # Store the startyear of the parset
        self.end = end # Store the endyear of the parset
        self.isfixed = None # Store whether props are fixed or not
        
    
    def __repr__(self):
        ''' Print out useful information when called'''
        output  = defaultrepr(self)
        output += 'Parameter set name: %s\n'    % self.name
        output += '      Date created: %s\n'    % getdate(self.created)
        output += '     Date modified: %s\n'    % getdate(self.modified)
        output += '               UID: %s\n'    % self.uid
        output += '============================================================\n'
        return output
    
    
    def getresults(self, die=True):
        ''' Method for getting the results '''
        if self.resultsref is not None and self.projectref() is not None:
            results = getresults(project=self.projectref(), pointer=self.resultsref, die=die)
            return results
        else:
            raise OptimaException('No results associated with this parameter set')
    
    
    def getcovpars(self):
        '''Method for getting a list of coverage-only parameters'''
        coveragepars = [par.short for par in self.pars.values() if isinstance(par, Par) and par.iscoveragepar()]
        return coveragepars
        

    def parkeys(self):
        ''' Return a list of the keys in pars that are actually parameter objects '''
        parslist = []
        for key,par in self.pars.items():
            if issubclass(type(par), Par):
                parslist.append(key)
        return parslist
    
    
    def makepars(self, data=None, fix=True, verbose=2, start=None, end=None):
        self.pars = makepars(data=data, verbose=verbose) # Initialize as list with single entry
        self.fixprops(fix=fix)
        self.popkeys = dcp(self.pars['popkeys']) # Store population keys more accessibly
        if start is None: self.start = data['years'][0] # Store the start year -- if not supplied, use beginning of data
        else:             self.start = start
        if end is None:   self.end   = Settings().endyear # Store the end year -- if not supplied, use default
        else:             self.end   = end
        return None


    def interp(self, keys=None, start=None, end=2030, dt=0.2, tvec=None, smoothness=20, asarray=True, samples=None, verbose=2):
        """ Prepares model parameters to run the simulation. """
        printv('Making model parameters...', 1, verbose),
        
        if start is None: start = self.start

        simparslist = []
        if isnumber(tvec): tvec = array([tvec]) # Convert to 1-element array -- WARNING, not sure if this is necessary or should be handled lower down
        if samples is None: samples = [None]
        for sample in samples:
            simpars = makesimpars(pars=self.pars, name=self.name, keys=keys, start=start, end=end, dt=dt, tvec=tvec, smoothness=smoothness, asarray=asarray, sample=sample, verbose=verbose)
            simparslist.append(simpars) # Wrap up
        
        return simparslist
    
    
    def updateprior(self):
        ''' Update the prior for all of the variables '''
        for key in self.parkeys():
            self.pars[key].updateprior()
        return None
    
    
    def fixprops(self, fix=None, which=None, startyear=None):
        '''
        Fix or unfix the proportions of people on ART and suppressed.
        
        To fix:   P.parset().fixprops()
        To unfix: P.parset().fixprops(False)
        
        You can also specify a start year. "which" can be a string
        or a list of strings, to specify which of ['dx', 'tx', 'supp']
        you want to fix.
        '''
        if fix is None:
            fix = True # By default, do fix
        self.isfixed = fix # Store fixed status
        if   which is None:  which = ['tx','supp']
        elif which is 'all': which = ['dx','tx','supp']
        else:                which = promotetolist(which)
        if startyear is None:
            if fix:  startyear = self.pars['numtx'].t['tot'][-1]
            else:    startyear = 2100
        if 'dx'   in which: self.pars['fixpropdx'].t   = startyear
        if 'tx'   in which: self.pars['fixproptx'].t   = startyear
        if 'supp' in which: self.pars['fixpropsupp'].t = startyear # Doesn't make sense to assume proportion on treatment without assuming proportion suppressed....also, crashes otherwise :)
        return None
        
    
    def usedataprops(self, use=None, which=None):
        '''
        Use the data in the Optional indicators tab to create parameters for proportions of people in cascade stages
        
        To fix:   P.parset().usedataprops()
        To unfix: P.parset().usedataprops(False)
        To fix specific props:   P.parset().usedataprops('supp')

        '''
        if use is None: use = True # By default, use the data
        if   which is None:  which = ['supp'] # By default, only use the indicator on the proportion suppressed
        elif which is 'all': which = ['dx','tx','supp']
        else:                which = promotetolist(which)
        
        data = self.projectref().data
        if use:
            for key in which:
                tmp = data2timepar(data=data['optprop'+key], years=data['years'], keys=self.pars['prop'+key].t.keys(), name='tmp', short='tmp')
                self.pars['prop'+key].y = tmp.y
                self.pars['prop'+key].t = tmp.t
        else:
            for key in which:
                self.pars['prop'+key].y[0] = array([nan])
                self.pars['prop'+key].t[0] = array([0.0])

        return None


    def printpars(self, output=False):
        outstr = ''
        count = 0
        for par in self.pars.values():
            if hasattr(par,'p'): print('WARNING, population size not implemented!')
            if hasattr(par,'y'):
                if hasattr(par.y, 'keys'):
                    count += 1
                    if len(par.keys())>1:
                        outstr += '%3i: %s\n' % (count, par.name)
                        for key in par.keys():
                            outstr += '     %s = %s\n' % (key, par.y[key])
                    elif len(par.keys())==1:
                        outstr += '%3i: %s = %s\n\n' % (count, par.name, par.y[0])
                    elif len(par.keys())==0:
                        outstr += '%3i: %s = (empty)' % (count, par.name)
                    else:
                        print('WARNING, not sure what to do with %s: %s' % (par.name, par.y))
                else:
                    count += 1
                    outstr += '%3i: %s = %s\n\n' % (count, par.name, par.y)
        print(outstr)
        if output: return outstr
        else: return None


    def listattributes(self):
        ''' Go through all the parameters and make a list of their possible attributes '''
        
        maxlen = 20
        pars = self.pars
        
        print('\n\n\n')
        print('CONTENTS OF PARS, BY TYPE:')
        partypes = []
        for key in pars: partypes.append(type(pars[key]))
        partypes = set(partypes)
        count1 = 0
        count2 = 0
        for partype in set(partypes): 
            count1 += 1
            print('  %i..%s' % (count1, str(partype)))
            for key in pars:
                if type(pars[key])==partype:
                    count2 += 1
                    print('      %i.... %s' % (count2, str(key)))
        
        print('\n\n\n')
        print('ATTRIBUTES:')
        attributes = {}
        for key in self.parkeys():
            theseattr = list(pars[key].__dict__.keys())
            for attr in theseattr:
                if attr not in attributes.keys(): attributes[attr] = []
                attributes[attr].append(getattr(pars[key], attr))
        for key in attributes:
            print('  ..%s' % key)
        print('\n\n')
        for key in attributes:
            count = 0
            print('  ..%s' % key)
            items = []
            for item in attributes[key]:
                try: 
                    string = str(item)
                    if string not in items: 
                        if len(string)>maxlen: string = string[:maxlen]
                        items.append(string) 
                except: 
                    items.append('Failed to append item')
            for item in items:
                count += 1
                print('      %i....%s' % (count, str(item)))
        return None


    def manualfitlists(self, parsubset=None, advanced=None):
        ''' WARNING -- not sure if this function is needed; if it is needed, it should be combined with manualgui,py '''
        if not self.pars:
            raise OptimaException("No parameters available!")
    
        # Check parname subset is valid
        if parsubset is None:
            tmppars = self.pars
        else:
            if type(parsubset)==str: parsubset=[parsubset]
            if parsubset and type(parsubset) not in (list, str):
                raise OptimaException("Expecting parsubset to be a list or a string!")
            for item in parsubset:
                if item not in [par.short for par in self.pars.values() if hasattr(par,'manual') and par.manual!='no']:
                    raise OptimaException("Parameter %s is not a manual parameter.")
            tmppars = {par.short:par for par in self.pars.values() if hasattr(par,'manual') and par.manual!='no' and par.short in parsubset}
            
        mflists = {'keys': [], 'subkeys': [], 'types': [], 'values': [], 'labels': []}
        keylist = mflists['keys']
        subkeylist = mflists['subkeys']
        typelist = mflists['types']
        valuelist = mflists['values']
        labellist = mflists['labels']
        
        for key in tmppars.keys():
            par = tmppars[key]
            if hasattr(par, 'manual') and par.manual != 'no':  # Don't worry if it doesn't work, not everything in tmppars is actually a parameter
                if par.manual=='meta':
                    if advanced: # By default, don't include these
                        keylist.append(key)
                        subkeylist.append(None)
                        typelist.append('meta')
                        valuelist.append(par.m)
                        labellist.append('%s: meta' % par.name)
                elif par.manual=='const':
                    keylist.append(key)
                    subkeylist.append(None)
                    typelist.append('const')
                    valuelist.append(par.y)
                    labellist.append(par.name)
                elif par.manual=='advanced': # These are also constants, but skip by default
                    if advanced:
                        keylist.append(key)
                        subkeylist.append(None)
                        typelist.append('const')
                        valuelist.append(par.y)
                        labellist.append(par.name)
                elif par.manual=='year':
                    keylist.append(key)
                    subkeylist.append(None)
                    typelist.append('year')
                    valuelist.append(par.t)
                    labellist.append(par.name)
                elif par.manual=='pop':
                    for subkey in par.keys():
                        keylist.append(key)
                        subkeylist.append(subkey)
                        typelist.append('pop')
                        valuelist.append(par.y[subkey])
                        labellist.append('%s: %s' % (par.name, str(subkey)))
                elif par.manual=='exp':
                    for subkey in par.keys():
                        keylist.append(key)
                        subkeylist.append(subkey)
                        typelist.append('exp')
                        valuelist.append(par.i[subkey])
                        labellist.append('%s: %s' % (par.name, str(subkey)))
                else:
                    print('Parameter type "%s" not implemented!' % par.manual)
    
        return mflists
    
    
    ## Define update step
    def update(self, mflists, verbose=2):
        ''' Update Parameterset with new results -- WARNING, duplicates the function in gui.py!!!! '''
        if not self.pars:
            raise OptimaException("No parameters available!")
    
        tmppars = self.pars
    
        keylist    = mflists['keys']
        subkeylist = mflists['subkeys']
        typelist   = mflists['types']
        valuelist  = mflists['values']
        
        ## Loop over all parameters and update them
        for (key, subkey, ptype, value) in zip(keylist, subkeylist, typelist, valuelist):
            if ptype=='meta':  # Metaparameters
                tmppars[key].m = float(value)
                printv('%s.m = %s' % (key, value), 4, verbose)
            elif ptype=='pop':  # Populations or partnerships
                tmppars[key].y[subkey] = float(value)
                printv('%s.y[%s] = %s' % (key, subkey, value), 4, verbose)
            elif ptype=='exp':  # Population growth
                tmppars[key].i[subkey] = float(value)
                printv('%s.i[%s] = %s' % (key, subkey, value), 4, verbose)
            elif ptype in ['const', 'advanced']:  # Constants
                tmppars[key].y = float(value)
                printv('%s.y = %s' % (key, value), 4, verbose)
            elif ptype=='year':  # Year parameters
                tmppars[key].t = float(value)
                printv('%s.t = %s' % (key, value), 4, verbose)
            else:
                errormsg = 'Parameter type "%s" not implemented!' % ptype
                raise OptimaException(errormsg)
    
                # parset.interp() and calculate results are supposed to be called from the outside
    
    def export(self, filename=None, folder=None, compare=None):
        '''
        Little function to export code for the current parameter set. To use, do something like:
        
        pars = P.pars()
        
        and then paste in the output of this function.
        
        If compare is not None, then only print out parameter values that differ. Most useful for
        comparing to default, e.g.
        P.parsets[-1].export(compare='default')
        '''
        cpars, cvalues = None, None
        if compare is not None:
            try: 
                cpars = self.projectref().parsets[compare].pars
            except: 
                print('Could not compare parset %s to parset %s; printing all parameters' % (self.name, compare))
                compare = None
        
        def oneline(values): return str(values).replace('\n',' ') 
        
        output = ''
        for parname,par in self.pars.items():
            prefix2 = None # Handle the fact that some parameters need more than one line to print
            values2 = None
            cvalues2 = None
            if hasattr(par,'manual'):
                if par.manual=='pop': 
                    values = par.y[:].tolist()
                    prefix = "pars['%s'].y[:] = " % parname
                    if cpars is not None: cvalues = cpars[parname].y[:].tolist()
                elif par.manual in ['const', 'advanced']: 
                    values = par.y
                    prefix = "pars['%s'].y = " % parname
                    if cpars is not None: cvalues = cpars[parname].y
                elif par.manual=='year': 
                    values = par.t
                    prefix = "pars['%s'].t = " % parname
                    if cpars is not None: cvalues = cpars[parname].t
                elif par.manual=='meta':
                    values = par.m
                    prefix = "pars['%s'].m = " % parname
                    if cpars is not None: cvalues = cpars[parname].m
                elif par.manual=='exp':
                    values  = par.i[:].tolist()
                    values2 = par.e[:].tolist()
                    prefix  = "pars['%s'].i[:] = " % parname
                    prefix2 = "pars['%s'].e[:] = " % parname
                    if cpars is not None: 
                        cvalues  = cpars[parname].i[:].tolist()
                        cvalues2 = cpars[parname].e[:].tolist()
                elif par.manual=='no':
                    values = None
                else: 
                    print('Parameter manual type "%s" not implemented' % par.manual)
                    values = None
                if values is not None:
                    if compare is None or (values!=cvalues) or (values2!=cvalues2):
                        output += prefix+oneline(values)+'\n'
                        if prefix2 is not None:
                            output += prefix2+oneline(values2)+'\n'
        
        if filename is not None:
            fullpath = makefilepath(filename=filename, folder=folder, default=self.name, ext='par')
            with open(fullpath, 'w') as f:
                f.write(output)
            return fullpath
        else:
            return output






#################################################################################################################################
### Define the other classes
#################################################################################################################################

class Par(object):
    '''
    The base class for epidemiological model parameters.
    
    There are four subclasses:
        * Constant objects store a single scalar value in y and an uncertainty sample in ysample -- e.g., transmfi
        * Metapar objects store an odict of y values, have a single metaparameter m, and an odict of ysample -- e.g., force
        * Timepar objects store an odict of y values, have a single metaparameter m, and uncertainty scalar msample -- e.g., condcas
        * Popsizepar objects are like Timepar objects except have odicts of i (intercept) and e (exponent) values
        * Yearpar objects store a single time value -- e.g., fixpropdx
    
    These four thus have different structures (where [] denotes dict):
        * Constants   have y, ysample
        * Metapars    have y[], ysample[], m, msample
        * Timepars    have y[], m, msample
        * Popsizepars have i[], e[], m, msample
        * Yearpars    have t
    
    Consequently, some of them have different sample(), updateprior(), and interp() methods; in brief:
        * Constants have sample() = ysample, interp() = y
        * Metapars have sample() = ysample[], interp() = m*y[]
        * Timepars have sample() = msample, interp() = m*y[]
        * Popsizepars have sample() = msample, interp() = m*i[]*exp(e[])
        * Yearpars have no sampling methods, and interp() = t
    
    Version: 2016nov06 
    '''
    def __init__(self, short=None, name=None, limits=(0.,1.), by=None, manual='', fromdata=None, m=1.0, prior=None, verbose=None, **defaultargs): # "type" data needed for parameter table, but doesn't need to be stored
        ''' To initialize with a prior, prior should be a dict with keys 'dist' and 'pars' '''
        self.short = short # The short name, e.g. "hivtest"
        self.name = name # The full name, e.g. "HIV testing rate"
        self.limits = limits # The limits, e.g. (0,1) -- a tuple since immutable
        self.by = by # Whether it's by population, partnership, or total
        self.manual = manual # Whether or not this parameter can be manually fitted: options are '', 'meta', 'pop', 'exp', etc...
        self.fromdata = fromdata # Whether or not the parameter is made from data
        self.m = m # Multiplicative metaparameter, e.g. 1
        self.msample = None # The latest sampled version of the metaparameter -- None unless uncertainty has been run, and only used for uncertainty runs 
        if prior is None:             self.prior = Dist() # Not supplied, create default distribution
        elif isinstance(prior, dict): self.prior = Dist(**prior) # Supplied as a dict, use it to create a distribution
        elif isinstance(prior, Dist): self.prior = prior # Supplied as a distribution, use directly
        else:
            errormsg = 'Prior must either be None, a Dist, or a dict with keys "dist" and "pars", not %s' % type(prior)
            raise OptimaException(errormsg)
    
    def __repr__(self):
        ''' Print out useful information when called'''
        output = defaultrepr(self)
        return output
    
    def iscoveragepar(self):
        ''' Determine whether it's a coverage parameter'''
        return True if self.limits[1] == 'maxpopsize' else False


class Constant(Par):
    ''' The definition of a single constant parameter, which may or may not vary by population '''
    
    def __init__(self, y=None, **defaultargs):
        Par.__init__(self, **defaultargs)
        del self.m # These don't exist for the Constant class
        del self.msample 
        self.y = y # y-value data, e.g. 0.3
        self.ysample = None # y-value data generated from the prior, e.g. 0.24353
    
    def keys(self):
        ''' Constants don't have any keys '''
        return None 
    
    def sample(self, randseed=None):
        ''' Recalculate ysample '''
        self.ysample = self.prior.sample(n=1, randseed=randseed)[0]
        return None
    
    def updateprior(self, verbose=2):
        ''' Update the prior parameters to match the metaparameter, so e.g. can recalibrate and then do uncertainty '''
        if self.prior.dist=='uniform':
            tmppars = array(self.prior.pars) # Convert to array for numerical magic
            self.prior.pars = sanitizeprior(self.prior.pars, self.y)
            printv('Priors updated for %s' % self.short, 3, verbose)
        else:
            errormsg = 'Distribution "%s" not defined; available choices are: uniform or bust, bro!' % self.dist
            raise OptimaException(errormsg)
        return None
    
    def interp(self, tvec=None, dt=None, smoothness=None, asarray=True, sample=False, randseed=None, popkeys=None): # Keyword arguments are for consistency but not actually used
        """
        Take parameters and turn them into model parameters -- here, just return a constant value at every time point
        
        There are however 3 options with the interpolation:
            * False -- use existing y value
            * 'old' -- use existing ysample value
            * 'new' -- recalculate ysample value
        """
        # Figure out sample
        if not sample: 
            y = self.y
        else:
            if sample=='new' or self.ysample is None: self.sample(randseed=randseed) # msample doesn't exist, make it
            y = self.ysample
            
        # Do interpolation
        dt = gettvecdt(tvec=tvec, dt=dt, justdt=True) # Method for getting dt
        output = applylimits(par=self, y=y, limits=self.limits, dt=dt)
        if not asarray: output = odict([('tot',output)])
        return output



class Metapar(Par):
    ''' The definition of a single metaparameter, such as force of infection, which usually does vary by population '''

    def __init__(self, y=None, prior=None, **defaultargs):
        Par.__init__(self, **defaultargs)
        self.y = y # y-value data, e.g. {'FSW:'0.3, 'MSM':0.7}
        self.ysample = None
        if isinstance(prior, dict):
            self.prior = prior
        elif prior is None:
            self.prior = odict()
            for key in self.keys():
                self.prior[key] = Dist() # Initialize with defaults
        else:
            errormsg = 'Prior for metaparameters must be an odict, not %s' % type(prior)
            raise OptimaException(errormsg)
            
    def keys(self):
        ''' Return the valid keys for using with this parameter '''
        return self.y.keys()
    
    def sample(self, randseed=None, verbose=False):
        ''' Recalculate ysample '''
        self.ysample = odict()
        for key in self.keys():
            self.ysample[key] = self.prior[key].sample(randseed=randseed)[0]
            if verbose: print(self.name, key, self.prior[key].pars, self.ysample[key]) # For debugging only
        return None
    
    def updateprior(self, verbose=2):
        ''' Update the prior parameters to match the y values, so e.g. can recalibrate and then do uncertainty '''
        for key in self.keys():
            if self.prior[key].dist=='uniform':
                self.prior[key].pars = sanitizeprior(self.prior[key].pars, self.y[key])
                printv('Priors updated for %s' % self.short, 3, verbose)
            else:
                errormsg = 'Distribution "%s" not defined; available choices are: uniform or bust, bro!' % self.dist
                raise OptimaException(errormsg)
        return None
    
    def interp(self, tvec=None, dt=None, smoothness=None, asarray=True, sample=None, randseed=None, popkeys=None): # Keyword arguments are for consistency but not actually used
        """ Take parameters and turn them into model parameters -- here, just return a constant value at every time point """
        
        # Figure out sample
        if not sample: 
            y = self.y
        else:
            if sample=='new' or self.ysample is None: self.sample(randseed=randseed) # msample doesn't exist, make it
            y = self.ysample
                
        dt = gettvecdt(tvec=tvec, dt=dt, justdt=True) # Method for getting dt
        outkeys = getoutkeys(self, popkeys) # Get the list of keys for the output
        if asarray: output = zeros(len(outkeys))
        else: output = odict()
        
        for pop,key in enumerate(outkeys): # Loop over each population, always returning an [npops x npts] array
            if key in self.keys(): yval = y[key]*self.m
            else:                  yval = 0. # Population not present, set to zero
            yinterp = applylimits(par=self, y=yval, limits=self.limits, dt=dt) 
            if asarray: output[pop] = yinterp
            else:       output[key] = yinterp
        return output
    


class Timepar(Par):
    ''' The definition of a single time-varying parameter, which may or may not vary by population '''
    
    def __init__(self, t=None, y=None, **defaultargs):
        Par.__init__(self, **defaultargs)
        if t is None: t = odict()
        if y is None: y = odict()
        self.t = t # Time data, e.g. [2002, 2008]
        self.y = y # Value data, e.g. [0.3, 0.7]
    
    def keys(self):
        ''' Return the valid keys for using with this parameter '''
        return self.y.keys()
    
    def df(self, key=None, data=None):
        '''
        Return t,y data as a data frame; or if data is supplied, replace current t,y values.
        Example: use df() to export data, work with it as a dataframe, and then import it back in:
        aidstest = P.pars()['aidstest'].df()
        aidstest.addrow([2005, 0.3])
        P.pars()['aidstest'].df(data=aidstest)
        '''
        if key is None: key = self.keys()[0] # Pull out first key if not specified -- e.g., 'tot'
        output = dataframe(['t','y'], [self.t[key], self.y[key]])
        if data is not None:
            if isinstance(data, dataframe):
                self.t[key] = array(data['t'],dtype=float)
                self.y[key] = array(data['y'],dtype=float)
            else:
                errormsg = 'Data argument must be a dataframe, not "%s"' % type(data)
                raise OptimaException(errormsg)
            return None
        else:
            return output
    
    def sample(self, randseed=None):
        ''' Recalculate msample '''
        self.msample = self.prior.sample(n=1, randseed=randseed)[0]
        return None
    
    def updateprior(self, verbose=2):
        ''' Update the prior parameters to match the metaparameter, so e.g. can recalibrate and then do uncertainty '''
        if self.prior.dist=='uniform':
            self.prior.pars = sanitizeprior(self.prior.pars, self.m)
            printv('Priors updated for %s' % self.short, 3, verbose)
        else:
            errormsg = 'Distribution "%s" not defined; available choices are: uniform or bust, bro!' % self.dist
            raise OptimaException(errormsg)
        return None
    
    def interp(self, tvec=None, dt=None, smoothness=None, asarray=True, sample=None, randseed=None, popkeys=None):
        """ Take parameters and turn them into model parameters """
        
        # Validate input
        if tvec is None: 
            errormsg = 'Cannot interpolate parameter "%s" with no time vector specified' % self.name
            raise OptimaException(errormsg)
        tvec, dt = gettvecdt(tvec=tvec, dt=dt) # Method for getting these as best possible
        if smoothness is None: smoothness = int(defaultsmoothness/dt) # Handle smoothness
        outkeys = getoutkeys(self, popkeys) # Get the list of keys for the output
        
        # Figure out sample
        if not sample:
            m = self.m
        else:
            if sample=='new' or self.msample is None: self.sample(randseed=randseed) # msample doesn't exist, make it
            m = self.msample
        
        # Set things up and do the interpolation
        npops = len(outkeys)
        if self.by=='pship': asarray= False # Force odict since too dangerous otherwise
        if asarray: output = zeros((npops,len(tvec)))
        else:       output = odict()
        for pop,key in enumerate(outkeys): # Loop over each population, always returning an [npops x npts] array
            if key in self.keys():
                yinterp = m * smoothinterp(tvec, self.t[key], self.y[key], smoothness=smoothness) # Use interpolation
                yinterp = applylimits(par=self, y=yinterp, limits=self.limits, dt=dt)
            else:
                yinterp = zeros(len(tvec)) # Population not present, just set to zero
            if asarray: output[pop,:] = yinterp
            else:       output[key]   = yinterp
        if npops==1 and self.by=='tot' and asarray: return output[0,:] # npops should always be 1 if by==tot, but just be doubly sure
        else: return output



class Popsizepar(Par):
    ''' The definition of the population size parameter '''
    
    def __init__(self, i=None, e=None, m=1.0, start=2000., **defaultargs):
        Par.__init__(self, **defaultargs)
        if i is None: i = odict()
        if e is None: e = odict()
        self.i = i # Exponential fit intercept, e.g. 3.4e6
        self.e = e # Exponential fit exponent, e.g. 0.03
        self.m = m # Multiplicative metaparameter, e.g. 1
        self.start = start # Year for which population growth start is calibrated to
    
    def keys(self):
        ''' Return the valid keys for using with this parameter '''
        return self.i.keys()
    
    def sample(self, randseed=None):
        ''' Recalculate msample -- same as Timepar'''
        self.msample = self.prior.sample(n=1, randseed=randseed)[0]
        return None
    
    def updateprior(self, verbose=2):
        ''' Update the prior parameters to match the metaparameter -- same as Timepar '''
        if self.prior.dist=='uniform':
            self.prior.pars = sanitizeprior(self.prior.pars, self.m)
            printv('Priors updated for %s' % self.short, 3, verbose)
        else:
            errormsg = 'Distribution "%s" not defined; available choices are: uniform or bust, bro!' % self.dist
            raise OptimaException(errormsg)
        return None

    def interp(self, tvec=None, dt=None, smoothness=None, asarray=True, sample=None, randseed=None, popkeys=None): # WARNING: smoothness isn't used, but kept for consistency with other methods...
        """ Take population size parameter and turn it into a model parameters """
        
        # Validate input
        if tvec is None: 
            errormsg = 'Cannot interpolate parameter "%s" with no time vector specified' % self.name
            raise OptimaException(errormsg)
        tvec, dt = gettvecdt(tvec=tvec, dt=dt) # Method for getting these as best possible
        outkeys = getoutkeys(self, popkeys) # Get the list of keys for the output
        
        # Figure out sample
        if not sample:
            m = self.m
        else:
            if sample=='new' or self.msample is None: self.sample(randseed=randseed) # msample doesn't exist, make it
            m = self.msample
        
        # Do interpolation
        npops = len(outkeys)
        if asarray: output = zeros((npops,len(tvec)))
        else: output = odict()
        for pop,key in enumerate(outkeys):
            if key in self.keys():
                yinterp = m * self.i[key] * grow(self.e[key], array(tvec)-self.start)
                yinterp = applylimits(par=self, y=yinterp, limits=self.limits, dt=dt)
            else:
                yinterp = zeros(len(tvec))
            if asarray: output[pop,:] = yinterp
            else:       output[key] = yinterp
        return output



class Yearpar(Par):
    ''' The definition of a single year parameter'''
    
    def __init__(self, t=None, **defaultargs):
        Par.__init__(self, **defaultargs)
        del self.m # These don't exist for this class
        del self.msample 
        self.prior = None
        self.t = t # y-value data, e.g. 0.3
    
    def keys(self):
        return None 
    
    def sample(self, randseed=None):
        ''' No sampling, so simply return the value '''
        return self.t
    
    def updateprior(self, verbose=2):
        '''No prior, so return nothing'''
        return None
    
    def interp(self, tvec=None, dt=None, smoothness=None, sample=None, randseed=None, asarray=True, popkeys=None): # Keyword arguments are for consistency but not actually used
        '''No interpolation, so simply return the value'''
        return self.t



class Dist(object):
    ''' Define a distribution object for drawing samples from, usually to create a prior '''
    def __init__(self, dist=None, pars=None):
        self.dist = dist if dist is not None else 'uniform'
        self.pars = promotetoarray(pars) if pars is not None else array([0.9, 1.1]) # This is arbitrary, of course
    
    def __repr__(self):
        ''' Print out useful information when called'''
        output = defaultrepr(self)
        return output
    
    def sample(self, n=1, randseed=None):
        ''' Draw random samples from the specified distribution '''
        if randseed is not None: seed(randseed) # Reset the random seed, if specified
        if self.dist=='uniform':
            samples = uniform(low=self.pars[0], high=self.pars[1], size=n)
            return samples
        if self.dist=='normal':
            return normal(loc=self.pars[0], scale=self.pars[1], size=n)
        else:
            errormsg = 'Distribution "%s" not defined; available choices are: uniform, normal' % self.dist
            raise OptimaException(errormsg)



#############################################################################################################################
### Functions for handling the parameters
#############################################################################################################################

def getoutkeys(par=None, popkeys=None):
    ''' Small method to decide whether to return 'tot', a subset of population keys, or all population keys '''
    if par.by in ['mpop','fpop'] and popkeys is not None:
        return popkeys # Expand male or female only keys to all
    else:
        return par.keys() # Or just return the default
            

def grow(exponent, tvec):
    ''' Return a time vector for a population growth '''
    return exp(tvec*exponent) # Simple exponential growth


def getvalidyears(years, validdata, defaultind=0):
    ''' Return the years that are valid based on the validity of the input data '''
    if sum(validdata): # There's at least one data point entered
        if len(years)==len(validdata): # They're the same length: use for logical indexing
            validyears = array(array(years)[validdata]) # Store each year
        elif len(validdata)==1: # They're different lengths and it has length 1: it's an assumption
            validyears = array([array(years)[defaultind]]) # Use the default index; usually either 0 (start) or -1 (end)
    else: validyears = array([0.0]) # No valid years, return 0 -- NOT an empty array, as you might expect!
    return validyears



def data2prev(data=None, keys=None, index=0, blh=0, **defaultargs): # WARNING, "blh" means "best low high", currently upper and lower limits are being thrown away, which is OK here...?
    """ Take an array of data return either the first or last (...or some other) non-NaN entry -- used for initial HIV prevalence only so far... """
    par = Metapar(y=odict([(key,None) for key in keys]), **defaultargs) # Create structure -- need key:None for prior
    for row,key in enumerate(keys):
        par.y[key] = sanitize(data['hivprev'][blh][row])[index] # Return the specified index -- usually either the first [0] or last [-1]
        par.prior[key].pars *= par.y[key] # Get prior in right range
    return par



def data2popsize(data=None, keys=None, blh=0, uniformgrowth=False, doplot=False, **defaultargs):
    ''' Convert population size data into population size parameters '''
    par = Popsizepar(m=1, **defaultargs)
    
    # Parse data into consistent form
    sanitizedy = odict() # Initialize to be empty
    sanitizedt = odict() # Initialize to be empty
    for row,key in enumerate(keys):
        sanitizedy[key] = sanitize(data['popsize'][blh][row]) # Store each extant value
        sanitizedt[key] = array(data['years'])[~isnan(data['popsize'][blh][row])] # Store each year
    
    # Store a list of population sizes that have at least 2 data points
    atleast2datapoints = [] 
    for key in keys:
        if len(sanitizedy[key])>=2:
            atleast2datapoints.append(key)
    if len(atleast2datapoints)==0:
        errormsg = 'Not more than one data point entered for any population size\n'
        errormsg += 'To estimate growth trends, at least one population must have at least 2 data points'
        raise OptimaException(errormsg)
        
    largestpopkey = atleast2datapoints[argmax([mean(sanitizedy[key]) for key in atleast2datapoints])] # Find largest population size (for at least 2 data points)
    
    # Perform 2-parameter exponential fit to data
    startyear = data['years'][0]
    par.start = data['years'][0]
    tdata = odict()
    ydata = odict()
    for key in atleast2datapoints:
        tdata[key] = sanitizedt[key]-startyear
        ydata[key] = log(sanitizedy[key])
        try:
            fitpars = polyfit(tdata[key], ydata[key], 1)
            par.i[key] = exp(fitpars[1]) # Intercept/initial value
            par.e[key] = fitpars[0] # Exponent
        except:
            errormsg = 'Fitting population size data for population "%s" failed' % key
            raise OptimaException(errormsg)
    
    # Handle populations that have only a single data point
    only1datapoint = list(set(keys)-set(atleast2datapoints))
    thisyear = odict()
    thispopsize = odict()
    for key in only1datapoint:
        largest_i = par.i[largestpopkey] # Get the parameters from the largest population
        largest_e = par.e[largestpopkey]
        if len(sanitizedt[key]) != 1:
            errormsg = 'Error interpreting population size for population "%s"\n' % key
            errormsg += 'Please ensure at least one time point is entered'
            raise OptimaException(errormsg)
        thisyear[key] = sanitizedt[key][0]
        thispopsize[key] = sanitizedy[key][0]
        largestthatyear = largest_i*grow(largest_e, thisyear[key]-startyear)
        par.i[key] = largest_i*thispopsize[key]/largestthatyear # Scale population size
        par.e[key] = largest_e # Copy exponent
    par.i.sort(keys) # Sort to regain the original key order -- WARNING, causes horrendous problems later if this isn't done!
    par.e.sort(keys)
    
    if uniformgrowth:
        for key in keys:
            par.e[key] = par.e[largestpopkey] # Reset exponent to match the largest population
            meanpopulationsize = mean(sanitizedy[key]) # Calculate the mean of all the data
            weightedyear = mean(sanitizedy[key][:]*sanitizedt[key][:])/meanpopulationsize # Calculate the "mean year"
            par.i[key] = meanpopulationsize*(1+par.e[key])**(startyear-weightedyear) # Project backwards to starting population size
    
    for key in keys:
        par.i[key] = round(par.i[key]) # Fractional people look weird
        
    if doplot:
        from pylab import figure, subplot, plot, scatter, arange, show, title
        nplots = len(par.keys())
        figure()
        tvec = arange(data['years'][0], data['years'][-1]+1)
        yvec = par.interp(tvec=tvec)
        for k,key in enumerate(par.keys()):
            subplot(nplots,1,k+1)
            if key in atleast2datapoints: scatter(tdata[key]+startyear, exp(ydata[key]))
            elif key in only1datapoint: scatter(thisyear[key], thispopsize[key])
            else: raise OptimaException('This population is nonexistent')
            plot(tvec, yvec[k])
            title('Pop size: ' + key)
            print([par.i[key], par.e[key]])
            show()
    
    return par



def data2timepar(data=None, years=None, keys=None, defaultind=0, verbose=2, **defaultargs):
    """ Take an array of data and turn it into default parameters -- here, just take the means """
    # Check that at minimum, name and short were specified, since can't proceed otherwise
    try: 
        name, short = defaultargs['name'], defaultargs['short']
    except: 
        errormsg = 'Cannot create a time parameter without keyword arguments "name" and "short"! \n\nArguments:\n %s' % defaultargs.items()
        raise OptimaException(errormsg)
        
    # Process data
    if isinstance(data,dict): # The entire structure has been passed
        thisdata = data[short]
        years = data['years']
    elif isinstance(data,list): # Just the relevant entry has been passed
        thisdata = data
        
    par = Timepar(m=1.0, y=odict(), t=odict(), **defaultargs) # Create structure
    for row,key in enumerate(keys):
        try:
            validdata = ~isnan(thisdata[row]) # WARNING, this could all be greatly simplified!!!! Shouldn't need to call this and sanitize()
            par.t[key] = getvaliddata(years, validdata, defaultind=defaultind) 
            if sum(validdata): 
                par.y[key] = sanitize(thisdata[row])
            else:
                printv('data2timepar(): no data for parameter "%s", key "%s"' % (name, key), 3, verbose) # Probably ok...
                par.y[key] = array([0.0]) # Blank, assume zero -- WARNING, is this ok?
                par.t[key] = array([0.0])
        except:
            errormsg = 'Error converting time parameter "%s", key "%s"' % (name, key)
            printv(errormsg, 1, verbose)
            raise

    return par


## Acts
def balance(act=None, which=None, data=None, popkeys=None, limits=None, popsizepar=None, eps=None):
    ''' 
    Combine the different estimates for the number of acts or condom use and return the "average" value.
    
    Set which='numacts' to compute for number of acts, which='condom' to compute for condom.
    '''
    if eps is None: eps = Settings().eps   # If not supplied (it won't be), get from default settings  
    
    if which not in ['numacts','condom']: raise OptimaException('Can only balance numacts or condom, not "%s"' % which)
    mixmatrix = array(data['part'+act]) # Get the partnerships matrix
    npops = len(popkeys) # Figure out the number of populations
    symmetricmatrix = zeros((npops,npops));
    for pop1 in range(npops):
        for pop2 in range(npops):
            if which=='numacts': symmetricmatrix[pop1,pop2] = symmetricmatrix[pop1,pop2] + (mixmatrix[pop1,pop2] + mixmatrix[pop2,pop1]) / float(eps+((mixmatrix[pop1,pop2]>0)+(mixmatrix[pop2,pop1]>0)))
            if which=='condom': symmetricmatrix[pop1,pop2] = bool(symmetricmatrix[pop1,pop2] + mixmatrix[pop1,pop2] + mixmatrix[pop2,pop1])
        
    # Decide which years to use -- use the earliest year, the latest year, and the most time points available
    yearstouse = []    
    for row in range(npops): yearstouse.append(getvaliddata(data['years'], data[which+act][row]))
    minyear = Inf
    maxyear = -Inf
    npts = 1 # Don't use fewer than 1 point
    for row in range(npops):
        minyear = minimum(minyear, min(yearstouse[row]))
        maxyear = maximum(maxyear, max(yearstouse[row]))
        npts = maximum(npts, len(yearstouse[row]))
    if minyear==Inf:  minyear = data['years'][0] # If not set, reset to beginning
    if maxyear==-Inf: maxyear = data['years'][-1] # If not set, reset to end
    ctrlpts = linspace(minyear, maxyear, npts).round() # Force to be integer...WARNING, guess it doesn't have to be?
    
    # Interpolate over population acts data for each year
    tmppar = data2timepar(name='tmp', short=which+act, limits=(0,'maxacts'), data=data[which+act], years=data['years'], keys=popkeys, by='pop', verbose=0) # Temporary parameter for storing acts
    tmpsim = tmppar.interp(tvec=ctrlpts)
    if which=='numacts': popsize = popsizepar.interp(tvec=ctrlpts)
    npts = len(ctrlpts)
    
    # Compute the balanced acts
    output = zeros((npops,npops,npts))
    for t in range(npts):
        if which=='numacts':
            smatrix = dcp(symmetricmatrix) # Initialize
            psize = popsize[:,t]
            popacts = tmpsim[:,t]
            for pop1 in range(npops): smatrix[pop1,:] = smatrix[pop1,:]*psize[pop1] # Yes, this needs to be separate! Don't try to put in the next for loop, the indices are opposite!
            for pop1 in range(npops): smatrix[:,pop1] = psize[pop1]*popacts[pop1]*smatrix[:,pop1] / float(eps+sum(smatrix[:,pop1])) # Divide by the sum of the column to normalize the probability, then multiply by the number of acts and population size to get total number of acts
        
        # Reconcile different estimates of number of acts, which must balance
        thispoint = zeros((npops,npops));
        for pop1 in range(npops):
            for pop2 in range(npops):
                if which=='numacts':
                    balanced = (smatrix[pop1,pop2] * psize[pop1] + smatrix[pop2,pop1] * psize[pop2])/(psize[pop1]+psize[pop2]) # here are two estimates for each interaction; reconcile them here
                    thispoint[pop2,pop1] = balanced/psize[pop2] # Divide by population size to get per-person estimate
                    thispoint[pop1,pop2] = balanced/psize[pop1] # ...and for the other population
                if which=='condom':
                    thispoint[pop1,pop2] = (tmpsim[pop1,t]+tmpsim[pop2,t])/2.0
                    thispoint[pop2,pop1] = thispoint[pop1,pop2]
    
        output[:,:,t] = thispoint
    
    return output, ctrlpts








def makepars(data=None, verbose=2, die=True, fixprops=None):
    """
    Translates the raw data (which were read from the spreadsheet) into
    parameters that can be used in the model. These data are then used to update 
    the corresponding model (project). This method should be called before a 
    simulation is run.
    
    Version: 2017jun03
    """
    
    printv('Converting data to parameters...', 1, verbose)
    
    ###############################################################################
    ## Loop over quantities
    ###############################################################################
    
    pars = odict()
    
    # Shorten information on which populations are male, which are female, which inject, which provide commercial sex
    pars['male'] = array(data['pops']['male']).astype(bool) # Male populations 
    pars['female'] = array(data['pops']['female']).astype(bool) # Female populations
    
    # Set up keys
    totkey = ['tot'] # Define a key for when not separated by population
    popkeys = data['pops']['short'] # Convert to a normal string and to lower case...maybe not necessary
    fpopkeys = [popkey for popno,popkey in enumerate(popkeys) if data['pops']['female'][popno]]
    mpopkeys = [popkey for popno,popkey in enumerate(popkeys) if data['pops']['male'][popno]]
    pars['popkeys'] = dcp(popkeys)
    pars['age'] = array(data['pops']['age'])
    
    
    # Read in parameters automatically
    try: 
        rawpars = loadpartable() # Read the parameters structure
    except OptimaException as E: 
        errormsg = 'Could not load parameter table: "%s"' % repr(E)
        raise OptimaException(errormsg)
        
    pars['fromto'], pars['transmatrix'] = loadtranstable(npops=len(popkeys)) # Read the transitions
        
    for rawpar in rawpars: # Iterate over all automatically read in parameters
        printv('Converting data parameter "%s"...' % rawpar['short'], 3, verbose)
        
        try: # Optionally keep going if some parameters fail
        
            # Shorten key variables
            partype = rawpar.pop('partype')
            parname = rawpar['short']
            by = rawpar['by']
            fromdata = rawpar['fromdata']
            rawpar['verbose'] = verbose # Easiest way to pass it in
            
            # Decide what the keys are
            if   by=='tot' : keys = totkey
            elif by=='pop' : keys = popkeys
            elif by=='fpop': keys = fpopkeys
            elif by=='mpop': keys = mpopkeys
            else:            keys = [] # They're not necessarily empty, e.g. by partnership, but too complicated to figure out here
            
            # Decide how to handle it based on parameter type
            if partype=='initprev': # Initialize prevalence only
                pars['initprev'] = data2prev(data=data, keys=keys, **rawpar) # Pull out first available HIV prevalence point
            
            elif partype=='popsize': # Population size only
                pars['popsize'] = data2popsize(data=data, keys=keys, **rawpar)
            
            elif partype=='timepar': # Otherwise it's a regular time par, made from data
                domake = False # By default, don't make the parameter
                if by!='pship' and fromdata: domake = True # If it's not a partnership parameter and it's made from data, then make it
                if domake:
                    pars[parname] = data2timepar(data=data, keys=keys, **rawpar) 
                else:
                    pars[parname] = Timepar(m=1.0, y=odict([(key,array([nan])) for key in keys]), t=odict([(key,array([0.0])) for key in keys]), **rawpar) # Create structure
            
            elif partype=='constant': # The constants, e.g. transmfi
                best = data[parname][0] if fromdata else nan
                low = data[parname][1] if fromdata else nan
                high = data[parname][2] if fromdata else nan
                thisprior = {'dist':'uniform', 'pars':(low, high)} if fromdata else None
                pars[parname] = Constant(y=best, prior=thisprior, **rawpar)
            
            elif partype=='meta': # Force-of-infection, inhomogeneity, relative HIV-related death rates, and transitions
                pars[parname] = Metapar(y=odict([(key,None) for key in keys]), **rawpar)
                
            elif partype=='yearpar': # Years to fix proportions of people at different cascade stages
                pars[parname] = Yearpar(t=nan, **rawpar)
            
        except Exception as E:
            errormsg = 'Failed to convert parameter %s:\n%s' % (parname, repr(E))
            if die: raise OptimaException(errormsg)
            else: printv(errormsg, 1, verbose)

    
    ###############################################################################
    ## Tidy up -- things that can't be converted automatically
    ###############################################################################
    
    # Birth transitions - these are stored as the proportion of transitions, which is constant, and is multiplied by time-varying birth rates in model.py
    npopkeys = len(popkeys)
    birthtransit = zeros((npopkeys,npopkeys))
    c = 0
    for pkno,popkey in enumerate(popkeys):
        if data['pops']['female'][pkno]: # WARNING, really ugly
            for colno,col in enumerate(data['birthtransit'][c]):
                if sum(data['birthtransit'][c]):
                    birthtransit[pkno,colno] = col/sum(data['birthtransit'][c])
            c += 1
    pars['birthtransit'] = birthtransit 

    # Aging transitions - these are time-constant
    agetransit = zeros((npopkeys,npopkeys))
    duration = array([age[1]-age[0]+1.0 for age in data['pops']['age']])
    for rowno,row in enumerate(data['agetransit']):
        if sum(row):
            for colno,colval in enumerate(row):
                if colval:
                    agetransit[rowno,colno] = sum(row)*duration[rowno]/colval
    pars['agetransit'] = agetransit

    # Risk transitions - these are time-constant
    pars['risktransit'] = array(data['risktransit'])
    
    # Circumcision
    for key in pars['numcirc'].keys():
        pars['numcirc'].y[key] = array([0.0]) # Set to 0 for all populations, since program parameter only
    
    # Fix treatment from final data year
    for key in ['fixproptx', 'fixpropsupp', 'fixpropdx', 'fixpropcare', 'fixproppmtct']:
        pars[key].t = 2100 # TODO: don't use these, so just set to (hopefully) well past the end of the analysis

    # Set the values of parameters that aren't from data
    pars['transnorm'].y = 0.43 # See analyses/misc/calculatecd4transnorm.py for calculation
    pars['transnorm'].prior.pars *= pars['transnorm'].y # Scale default range
    for key in popkeys: # Define values for each population
        pars['force'].y[key] = 1.0
        pars['hivdeath'].y[key] = 1.0
        pars['inhomo'].y[key] = 0.0
        pars['inhomo'].prior[key].pars = array([0.0, 0.3]) # Arbitrary
    
    # Impose limits on force and transnorm so their values don't get too extreme (note, force.m functions identically to transnorm.y, but use the latter)
    for foipar in ['force','transnorm']:
        pars[foipar].limits = (0.05, 50) # Arbitrary
    
    
    # Handle acts
    tmpacts = odict()
    tmpcond = odict()
    tmpactspts = odict()
    tmpcondpts = odict()
    for act in ['reg','cas','com', 'inj']: # Number of acts
        actsname = 'acts'+act
        tmpacts[act], tmpactspts[act] = balance(act=act, which='numacts', data=data, popkeys=popkeys, popsizepar=pars['popsize'])
    for act in ['reg','cas','com']: # Condom use
        condname = 'cond'+act
        tmpcond[act], tmpcondpts[act] = balance(act=act, which='condom', data=data, popkeys=popkeys)
        
    # Convert matrices to lists of of population-pair keys
    for act in ['reg', 'cas', 'com', 'inj']: # Will probably include birth matrices in here too...
        actsname = 'acts'+act
        condname = 'cond'+act
        for i,key1 in enumerate(popkeys):
            for j,key2 in enumerate(popkeys):
                if sum(array(tmpacts[act])[i,j,:])>0:
                    pars[actsname].y[(key1,key2)] = array(tmpacts[act])[i,j,:]
                    pars[actsname].t[(key1,key2)] = array(tmpactspts[act])
                    if act!='inj':
                        if key1 in mpopkeys or key1 not in fpopkeys: # For condom use, only store one of the pair -- and store male first -- WARNING, would this fail with multiple MSM populations?
                            pars[condname].y[(key1,key2)] = array(tmpcond[act])[i,j,:]
                            pars[condname].t[(key1,key2)] = array(tmpcondpts[act])
    
    # Store information about injecting populations -- needs to be here since relies on other calculations
    pars['injects'] = array([pop in [pop1 for (pop1,pop2) in pars['actsinj'].keys()] for pop in pars['popkeys']])
    
    return pars



def makesimpars(pars, name=None, keys=None, start=None, end=None, dt=None, tvec=None, settings=None, smoothness=None, asarray=True, sample=None, tosample=None, randseed=None, verbose=2):
    ''' 
    A function for taking a single set of parameters and returning the interpolated versions -- used
    very directly in Parameterset.
    
    Version: 2017mar01
    '''
    
    # Handle inputs and initialization
    simpars = odict() 
    simpars['parsetname'] = name
    if keys is None: keys = list(pars.keys()) # Just get all keys
    if type(keys)==str: keys = [keys] # Listify if string
    if tvec is not None: simpars['tvec'] = tvec
    elif settings is not None: simpars['tvec'] = settings.maketvec(start=start, end=end, dt=dt)
    else: simpars['tvec'] = inclusiverange(start=start, stop=end, step=dt) # Store time vector with the model parameters
    if len(simpars['tvec'])>1: dt = simpars['tvec'][1] - simpars['tvec'][0] # Recalculate dt since must match tvec
    simpars['dt'] = dt  # Store dt
    if smoothness is None: smoothness = int(defaultsmoothness/dt)
    tosample = promotetolist(tosample) # Convert to list
    popkeys = pars['popkeys'] # Used for interpolation
    
    # Copy default keys by default
    for key in generalkeys: simpars[key] = dcp(pars[key])
    for key in staticmatrixkeys: simpars[key] = dcp(array(pars[key]))

    # Loop over requested keys
    for key in keys: # Loop over all keys
        if isinstance(pars[key], Par): # Check that it is actually a parameter -- it could be the popkeys odict, for example
            thissample = sample # Make a copy of it to check it against the list of things we are sampling
            if tosample and tosample[0] is not None and key not in tosample: thissample = False # Don't sample from unselected parameters -- tosample[0] since it's been promoted to a list
            try:
                simpars[key] = pars[key].interp(tvec=simpars['tvec'], dt=dt, popkeys=popkeys, smoothness=smoothness, asarray=asarray, sample=thissample, randseed=randseed)
            except OptimaException as E: 
                errormsg = 'Could not figure out how to interpolate parameter "%s"' % key
                errormsg += 'Error: "%s"' % repr(E)
                raise OptimaException(errormsg)


    return simpars




def applylimits(y, par=None, limits=None, dt=None, warn=True, verbose=2):
    ''' 
    A function to intelligently apply limits (supplied as [low, high] list or tuple) to an output.
    
    Needs dt as input since that determines maxrate.
    
    Version: 2016jan30
    '''
    
    # If parameter object is supplied, use it directly
    parname = ''
    if par is not None:
        if limits is None: limits = par.limits
        parname = par.name
        
    # If no limits supplied, don't do anything
    if limits is None:
        printv('No limits supplied for parameter "%s"' % parname, 4, verbose)
        return y
    
    if dt is None:
        if warn: raise OptimaException('No timestep specified: required for convertlimits()')
        else: dt = 0.2 # WARNING, should probably not hard code this, although with the warning, and being conservative, probably OK
    
    # Convert any text in limits to a numerical value
    limits = convertlimits(limits=limits, dt=dt, verbose=verbose)
    
    # Apply limits, preserving original class -- WARNING, need to handle nans
    if isnumber(y):
        if ~isfinite(y): return y # Give up
        newy = median([limits[0], y, limits[1]])
        if warn and newy!=y: printv('Note, parameter value "%s" reset from %f to %f' % (parname, y, newy), 3, verbose)
    elif shape(y):
        newy = array(y) # Make sure it's an array and not a list
        infiniteinds = findinds(~isfinite(newy))
        infinitevals = newy[infiniteinds] # Store these for safe keeping
        if len(infiniteinds): newy[infiniteinds] = limits[0] # Temporarily reset -- value shouldn't matter
        newy[newy<limits[0]] = limits[0]
        newy[newy>limits[1]] = limits[1]
        newy[infiniteinds] = infinitevals # And stick them back in
        if warn and any(newy!=array(y)):
            printv('Note, parameter "%s" value reset from:\n%s\nto:\n%s' % (parname, y, newy), 3, verbose)
    else:
        if warn: raise OptimaException('Data type "%s" not understood for applying limits for parameter "%s"' % (type(y), parname))
        else: newy = array(y)
    
    if shape(newy)!=shape(y):
        errormsg = 'Something went wrong with applying limits for parameter "%s":\ninput and output do not have the same shape:\n%s vs. %s' % (parname, shape(y), shape(newy))
        raise OptimaException(errormsg)
    
    return newy





def comparepars(pars1=None, pars2=None):
    ''' 
    Function to compare two sets of pars. Example usage:
    comparepars(P.parsets[0], P.parsets[1])
    '''
    if type(pars1)==Parameterset: pars1 = pars1.pars # If parset is supplied instead of pars, use that instead
    if type(pars2)==Parameterset: pars2 = pars2.pars
    keys = list(pars1.keys())
    nkeys = 0
    count = 0
    for key in keys:
        if hasattr(pars1[key],'y'):
            nkeys += 1
            if str(pars1[key].y) != str(pars2[key].y): # Convert to string representation for testing equality
                count += 1
                msg = 'Parameter "%s" differs:\n' % key
                msg += '%s\n' % pars1[key].y
                msg += 'vs\n'
                msg += '%s\n' % pars2[key].y
                msg += '\n\n'
                print(msg)
    if count==0: print('All %i parameters match' % nkeys)
    else:        print('%i of %i parameters did not match' % (count, nkeys))
    return None



def comparesimpars(pars1=None, pars2=None, inds=Ellipsis, inds2=Ellipsis):
    ''' 
    Function to compare two sets of simpars, like what's stored in results.
    
    Example:
        import optima as op
        P = op.demo(0)
        P.copyparset(0,'new')
        P.pars('new')['numtx'].y[:] *= 1.5
        R1 = P.runsim('default', keepraw=True)
        R2 = P.runsim('new', keepraw=True)
        op.comparesimpars(R1.simpars, R2.simpars)
    '''
    if type(pars1)==list: pars1 = pars1[0] # If a list is supplied, pull out just the dict
    if type(pars2)==list: pars2 = pars2[0]
    keys = pars1.keys()
    nkeys = 0
    count = 0
    for key in keys:
        nkeys += 1
        thispar1 = pars1[key]
        thispar2 = pars2[key]
        if isinstance(thispar1,dict): keys2 = thispar1.keys()
        else: keys2 = [None]
        for key2 in keys2:
            if key2 is not None:
                this1 = array(thispar1[key2])
                this2 = array(thispar2[key2])
                key2str = '(%s)' % str(key2)
            else:
                this1 = array(thispar1)
                this2 = array(thispar2)
                key2str = ''
            if len(shape(this1))==2:
                pars1str = str(this1[inds2][inds])
                pars2str = str(this2[inds2][inds])
            elif len(shape(this1))==1:
                pars1str = str(this1[inds])
                pars2str = str(this2[inds])
            else:
                pars1str = str(this1)
                pars2str = str(this2)
            if pars1str != pars2str: # Convert to string representation for testing equality
                count += 1
                dividerlen = 70
                bigdivide    = '='*dividerlen+'\n'
                littledivide = '-'*int(dividerlen/2.0-4)
                msg  = '\n\n'+bigdivide
                msg += 'Parameter "%s" %s differs:\n\n' % (key, key2str)
                msg += '%s\n' % pars1str
                msg += littledivide + ' vs ' + littledivide + '\n'
                msg += '%s\n\n' % pars2str
                msg += bigdivide
                print(msg)
    if count==0: print('All %i parameters match' % nkeys)
    else:        print('%i of %i parameters did not match' % (count, nkeys))
    return None


def sanitycheck(simpars=None, showdiff=True, threshold=0.1, eps=1e-6):
    '''
    Compare the current simpars with the default simpars, flagging
    potential differences. If simpars is None, generate it from the
    current parset. If showdiff is True, only show parameters that differ
    by more than the threshold amount (default, 10%). eps is just to
    avoid divide-by-zero errors and can be ignored, probably.
    
    Usage:
        sanitycheck(P)
    or
        result = P.runsim(keepraw=True) # Need keepraw or else it doesn't store simpars
        sanitycheck(result.simpars)
    '''
    if isinstance(simpars, op.Project): # It's actually a project
        thisproj = simpars # Rename so it's clearer
        try:    simpars = thisproj.result().simpars # Try to extract the simpars
        except: simpars = thisproj.runsim(keepraw=True, die=False).simpars # If not, rerun
            
    tmpproj = op.demo(dorun=False, doplot=False) # Can't import this earlier since not actually declared before
    tmpproj.runsim(keepraw=True)
    gsp = op.dcp(tmpproj.results[-1].simpars[0]) # "Good simpars"
    sp = simpars[0] # "Simpars"
    
    if set(sp.keys())!=set(gsp.keys()):
        errormsg = 'Keys do not match:'
        errormsg += 'Too many: %s' % (set(sp.keys())-set(gsp.keys()))
        errormsg += 'Missing: %s' % (set(gsp.keys())-set(sp.keys()))
        raise op.OptimaException(errormsg)
    
    outstr = ''
    skipped = []
    for k,key in enumerate(gsp.keys()):
        if op.checktype(sp[key], 'number'):
            spval = sp[key]
            gspval = gsp[key]
            ratio = (eps+spval)/(eps+gspval)
            if not showdiff or not(abs(1-ratio)<threshold):
                outstr += '\n\n%s\n%s (%03i/%03i)\n' % ('='*70, key, k, len(sp.keys())-1)
                outstr += 'Yours: %10s  Best: %10s Ratio: %10s\n' % (op.sigfig(spval), op.sigfig(gspval), op.sigfig(ratio))
            else:
                skipped.append(key)
        elif op.checktype(sp[key], 'arraylike') or op.checktype(sp[key], op.odict):
            if op.checktype(sp[key], op.odict):
                sp[key] = sp[key][:] # Try converting to an odict...
                gsp[key] = gsp[key][:] # Try converting to an odict...
            spmin = sp[key].min()
            spmax = sp[key].max()
            gspmin = gsp[key].min()
            gspmax = gsp[key].max()
            minratio = (eps+spmin)/(eps+gspmin)
            maxratio = (eps+spmax)/(eps+gspmax)
            if not showdiff or not(abs(1-minratio)<threshold) or not(abs(1-maxratio)<threshold):
                outstr += '\n\n%s\n%s (%03i/%03i)\n' % ('='*70, key, k, len(sp.keys())-1)
                outstr += 'Min-Yours: %10s  Min-Best: %10s Min-Ratio: %10s\n' % (op.sigfig(spmin), op.sigfig(gspmin), op.sigfig(minratio))
                outstr += 'Max-Yours: %10s  Max-Best: %10s Max-Ratio: %10s\n' % (op.sigfig(spmax), op.sigfig(gspmax), op.sigfig(maxratio))
            else:
                skipped.append(key)
        else:
            if not showdiff:
                outstr += '\n\n%s\n%s (%03i/%03i)\n' % ('='*70, key, k, len(sp.keys())-1)
                outstr += str(type(sp[key]))
            else:
                skipped.append(key)
        
    print(outstr)
    return outstr


def sanitizeprior(pars, m):
    tmppars = array(pars) # Convert to array for numerical magic
    if (tmppars==0).all():
        tmppars = array([0.9, 1.1]) # Reset
    tmpparssum = max(Settings().eps, tmppars.mean()) # Ensure it's nonzero
    newpars = tuple(m*tmppars/tmpparssum) # Recenter the limits around the mean
    return newpars













            
                
from optima import odict, getdate, today, uuid, dcp, objrepr, printv, scaleratio, OptimaException, findinds # Import utilities
from optima import gitinfo, tic, toc # Import functions
from optima import __version__ # Get current version

from optima import defaultobjectives, asd, Project

from numpy import arange

#######################################################################################################
## Portfolio class -- this contains Projects and GA optimisations
#######################################################################################################

budgeteps = 1e-8        # Project optimisations will fail for budgets that are optimised by GA to be zero. This avoids zeros.
tol = 1.0 # Tolerance for checking that budgets match


class Portfolio(object):
    """
    PORTFOLIO

    The super Optima portfolio class.

    Version: 2016jan20 by davidkedz
    """
    
    #######################################################################################################
    ## Built-in methods -- initialization, and the thing to print if you call a portfolio
    #######################################################################################################

    def __init__(self, name='default', projects=None, gaoptims=None):
        ''' Initialize the portfolio '''

        ## Set name
        self.name = name

        ## Define the structure sets
        self.projects = odict()
        if projects is not None: self.addprojects(projects)
        self.gaoptims = gaoptims if gaoptims else odict()

        ## Define metadata
        self.uid = uuid()
        self.created = today()
        self.modified = today()
        self.version = __version__
        self.gitbranch, self.gitversion = gitinfo()

        return None


    def __repr__(self):
        ''' Print out useful information when called '''
        output = '============================================================\n'
        output += '            Portfolio name: %s\n' % self.name
        output += '\n'
        output += '        Number of projects: %i\n' % len(self.projects)
        output += 'Number of GA Optimizations: %i\n' % len(self.gaoptims)
        output += '\n'
        output += '            Optima version: %0.1f\n' % self.version
        output += '              Date created: %s\n'    % getdate(self.created)
        output += '             Date modified: %s\n'    % getdate(self.modified)
        output += '                Git branch: %s\n'    % self.gitbranch
        output += '               Git version: %s\n'    % self.gitversion
        output += '                       UID: %s\n'    % self.uid
        output += '============================================================\n'
        output += objrepr(self)
        return output


    #######################################################################################################
    ## Methods to handle common tasks
    #######################################################################################################

    def addprojects(self, projects, verbose=2):
        ''' Store a project within portfolio '''
        printv('Adding project to portfolio...', 2, verbose)
        if type(projects)==Project: projects = [projects]
        if type(projects)==list:
            for project in projects:
                print str(project.uid) # TEMPPPP
                self.projects[str(project.uid)] = project        
                printv('\nAdded project "%s" to portfolio "%s".' % (project.name, self.name), 2, verbose)
        
    def getdefaultbudgets(self, progsetnames=None, verbose=2):
        ''' Get the default allocation totals of each project, using the progset names or indices specified '''
        budgets = []
        printv('Getting budgets...', 2, verbose)
        
        # Validate inputs
        if progsetnames==None:
            printv('\nWARNING: no progsets specified. Using default budget from first saved progset for each project for portfolio "%s".' % (self.name), 4, verbose)
            progsetnames = [0]*len(self.projects)
        if not len(progsetnames)==len(self.projects):
            printv('WARNING: %i program set names/indices were provided, but portfolio "%s" contains %i projects. OVERWRITING INPUTS and using default budget from first saved progset for each project.' % (len(progsetnames), self.name, len(self.projects)), 4, verbose)
            progsetnames = [0]*len(self.projects)

        # Loop over projects & get defaul budget for each, if you can
        for pno, p in enumerate(self.projects.values()):

            # Crash if any project doesn't have progsets
            if not p.progsets: 
                errormsg = 'Project "%s" does not have a progset. Cannot get default budgets.'
                raise OptimaException(errormsg)

            # Check that the progsets that were specified are indeed valid. They could be a string or a list index, so must check both
            if isinstance(progsetnames[pno],str) and progsetnames[pno] not in [progset.name for progset in p.progsets]:
                printv('\nCannot find progset "%s" in project "%s". Using progset "%s" instead.' % (progsetnames[pno], p.name, p.progsets[progsetnames[0]].name), 3, verbose)
                pno=0
            elif isinstance(progsetnames[pno],int) and len(p.progsets)<=progsetnames[pno]:
                printv('\nCannot find progset number %i in project "%s", there are only %i progsets in that project. Using progset 0 instead.' % (progsetnames[pno], p.name, len(p.progsets)), 1, verbose)
                pno=0
            else: 
                printv('\nCannot understand what program set to use for project "%s". Using progset 0 instead.' % (p.name), 3, verbose)
                pno=0            
                
            printv('\nAdd default budget from progset "%s" for project "%s" and portfolio "%s".' % (p.progsets[progsetnames[pno]].name, p.name, self.name), 4, verbose)
            budgets.append(sum(p.progsets[progsetnames[pno]].getdefaultbudget().values()))
        
        return budgets
    
    
    #######################################################################################################
    ## Methods to perform major tasks
    #######################################################################################################
        
        
    def genBOCs(self, objectives=None, progsetnames=None, parsetnames=None, maxtime=None, forceregen=False, verbose=2):
        ''' Loop through stored projects and construct budget-outcome curves '''
        printv('Generating BOCs...', 1, verbose)
        
        # Validate inputs
        if objectives == None: 
            printv('genBOCs(): WARNING, you have called genBOCs on portfolio %s without specifying obejctives. Using default objectives... ' % (self.name), 2, verbose)
            objectives = defaultobjectives()
        if progsetnames==None:
            printv('\ngenBOCs(): WARNING: no progsets specified. Using first saved progset for each project for portfolio "%s".' % (self.name), 3, verbose)
            progsetnames = [0]*len(self.projects)
        if not len(progsetnames)==len(self.projects):
            printv('genBOCs(): WARNING: %i program set names/indices were provided, but portfolio "%s" contains %i projects. OVERWRITING INPUTS and using first saved progset for each project.' % (len(progsetnames), self.name, len(self.projects)), 1, verbose)
            progsetnames = [0]*len(self.projects)
        if parsetnames==None:
            printv('\ngenBOCs(): WARNING: no parsets specified. Using first saved parset for each project for portfolio "%s".' % (self.name), 3, verbose)
            parsetnames = [0]*len(self.projects)
        if not len(parsetnames)==len(self.projects):
            printv('genBOCs(): WARNING: %i parset names/indices were provided, but portfolio "%s" contains %i projects. OVERWRITING INPUTS and using first saved parset for each project.' % (len(parsetnames), self.name, len(self.projects)), 1, verbose)
            parsetnames = [0]*len(self.projects)

        for pno, p in enumerate(self.projects.values()):
            getBOCtest = (p.getBOC(objectives) == None)
            if getBOCtest or forceregen:
                printv('genBOCs(): Regenerating BOC because getBOC=%s, forceregen=%s ' % (getBOCtest, forceregen), 2, verbose)

                # Crash if any project doesn't have progsets
                if not p.progsets or not p.parsets: 
                    errormsg = 'genBOCs(): Project "%s" does not have a progset and/or a parset, can''t generate a BOC.'
                    raise OptimaException(errormsg)
    
                # Check that the progsets that were specified are indeed valid. They could be a string or a list index, so must check both
                if isinstance(progsetnames[pno],str) and progsetnames[pno] not in [progset.name for progset in p.progsets.values()]:
                    printv('\ngenBOCs(): Cannot find progset "%s" in project "%s". Using progset "%s" instead.' % (progsetnames[pno], p.name, p.progsets[0].name), 1, verbose)
                    pno=0
                elif isinstance(progsetnames[pno],int) and len(p.progsets)<=progsetnames[pno]:
                    printv('\ngenBOCs(): Cannot find progset number %i in project "%s", there are only %i progsets in that project. Using progset 0 instead.' % (progsetnames[pno], p.name, len(p.progsets)), 1, verbose)
                    pno=0
                else: 
                    printv('\ngenBOCs(): Cannot understand what program set to use for project "%s". Using progset 0 instead.' % (p.name), 3, verbose)
                    pno=0            

                # Check that the progsets that were specified are indeed valid. They could be a string or a list index, so must check both
                if isinstance(parsetnames[pno],str) and parsetnames[pno] not in [parset.name for parset in p.parsets.values()]:
                    printv('\ngenBOCs(): Cannot find parset "%s" in project "%s". Using pargset "%s" instead.' % (progsetnames[pno], p.name, p.parsets[0].name), 1, verbose)
                    pno=0
                elif isinstance(parsetnames[pno],int) and len(p.parsets)<=parsetnames[pno]:
                    printv('\ngenBOCs(): Cannot find parset number %i in project "%s", there are only %i parsets in that project. Using parset 0 instead.' % (parsetnames[pno], p.name, len(p.parsets)), 1, verbose)
                    pno=0
                else: 
                    printv('\ngenBOCs(): Cannot understand what parset to use for project "%s". Using parset 0 instead.' % (p.name), 3, verbose)
                    pno=0            

                # Actually generate te BOCs
                printv('genBOCs(): WARNING, project %s does not have BOC, or else it does but you want to regenerate it. Generating one using parset %s and progset %s... ' % (p.name, p.parsets[parsetnames[pno]].name, p.progsets[progsetnames[pno]].name), 1, verbose)
                p.delBOC(objectives)    # Delete BOCs in case forcing regeneration.
                p.genBOC(parsetname=p.parsets[parsetnames[pno]].name, progsetname=p.progsets[progsetnames[pno]].name, objectives=objectives, maxtime=maxtime)

            else:
                printv('genBOCs(): Project %s contains a BOC, no need to generate... ' % p.name, 2, verbose)
                
                
    def plotBOCs(self, objectives=None, initbudgets=None, optbudgets=None, deriv=False, verbose=2):
        ''' Loop through stored projects and plot budget-outcome curves '''
        printv('Plotting BOCs...', 2, verbose)
        
        if initbudgets == None: initbudgets = [None]*len(self.projects)
        if optbudgets == None: optbudgets = [None]*len(self.projects)
        if objectives == None: 
            printv('WARNING, you have called plotBOCs on portfolio %s without specifying objectives. Using default objectives... ' % (self.name), 2, verbose)
            objectives = defaultobjectives()
            
        if not len(self.projects) == len(initbudgets) or not len(self.projects) == len(optbudgets):
            errormsg = 'Error: Plotting BOCs for %i projects with %i initial budgets (%i required) and %i optimal budgets (%i required).' % (len(self.projects), len(initbudgets), len(self.projects), len(optbudgets), len(self.projects))
            raise OptimaException(errormsg)
        
        # Loop for BOCs and then BOC derivatives.
        for c,p in enumerate(self.projects.values()):
            p.plotBOC(objectives=objectives, deriv=deriv, initbudget=initbudgets[c], optbudget=optbudgets[c])
            
            
    def minBOCoutcomes(self, objectives, progsetnames=None, parsetnames=None, seedbudgets=None, maxtime=None, verbose=2):
        ''' Loop through project BOCs corresponding to objectives and minimise net outcome '''
        printv('Calculating minimum BOC outcomes...', 2, verbose)

        # Check inputs
        if objectives == None: 
            printv('WARNING, you have called minBOCoutcomes on portfolio %s without specifying obejctives. Using default objectives... ' % (self.name), 2, verbose)
            objectives = defaultobjectives()
        if progsetnames==None:
            printv('\nWARNING: no progsets specified. Using first saved progset for each project for portfolio "%s".' % (self.name), 3, verbose)
            progsetnames = [0]*len(self.projects)
        if not len(progsetnames)==len(self.projects):
            printv('WARNING: %i program set names/indices were provided, but portfolio "%s" contains %i projects. OVERWRITING INPUTS and using first saved progset for each project.' % (len(progsetnames), self.name, len(self.projects)), 1, verbose)
            progsetnames = [0]*len(self.projects)
        if parsetnames==None:
            printv('\nWARNING: no parsets specified. Using first saved parset for each project for portfolio "%s".' % (self.name), 3, verbose)
            parsetnames = [0]*len(self.projects)
        if not len(parsetnames)==len(self.projects):
            printv('WARNING: %i parset names/indices were provided, but portfolio "%s" contains %i projects. OVERWRITING INPUTS and using first saved parset for each project.' % (len(parsetnames), self.name, len(self.projects)), 1, verbose)
            parsetnames = [0]*len(self.projects)
        
        # Initialise internal parameters
        BOClist = []
        grandtotal = objectives['budget']
        
        # Scale seedbudgets just in case they don't add up to the required total.
        if not seedbudgets == None:
            seedbudgets = scaleratio(seedbudgets, objectives['budget'])
            
        for pno,p in enumerate(self.projects.values()):
            
            if p.getBOC(objectives) is None:

                # Crash if any project doesn't have progsets
                if not p.progsets or not p.parsets: 
                    errormsg = 'Project "%s" does not have a progset and/or a parset, can''t generate a BOC.'
                    raise OptimaException(errormsg)
    
                # Check that the progsets that were specified are indeed valid. They could be a string or a list index, so must check both
                if isinstance(progsetnames[pno],str) and progsetnames[pno] not in [progset.name for progset in p.progsets]:
                    printv('\nCannot find progset "%s" in project "%s". Using progset "%s" instead.' % (progsetnames[pno], p.name, p.progsets[0].name), 3, verbose)
                    pno=0
                elif isinstance(progsetnames[pno],int) and len(p.progsets)<=progsetnames[pno]:
                    printv('\nCannot find progset number %i in project "%s", there are only %i progsets in that project. Using progset 0 instead.' % (progsetnames[pno], p.name, len(p.progsets)), 1, verbose)
                    pno=0
                else: 
                    printv('\nCannot understand what program set to use for project "%s". Using progset 0 instead.' % (p.name), 3, verbose)
                    pno=0            
    
                # Check that the parsets that were specified are indeed valid. They could be a string or a list index, so must check both
                if isinstance(parsetnames[pno],str) and parsetnames[pno] not in [parset.name for parset in p.parsets]:
                    printv('\nCannot find parset "%s" in project "%s". Using pargset "%s" instead.' % (progsetnames[pno], p.name, p.parsets[0].name), 1, verbose)
                    pno=0
                elif isinstance(parsetnames[pno],int) and len(p.parsets)<=parsetnames[pno]:
                    printv('\nCannot find parset number %i in project "%s", there are only %i parsets in that project. Using parset 0 instead.' % (parsetnames[pno], p.name, len(p.parsets)), 1, verbose)
                    pno=0
                else: 
                    printv('\nCannot understand what parset to use for project "%s". Using parset 0 instead.' % (p.name), 3, verbose)
                    pno=0
                
                    printv('WARNING, project "%s", parset "%s" does not have BOC. Generating one using parset %s and progset %s... ' % (p.name, pno, p.parsets[0].name, p.progsets[0].name), 1, verbose)
                p.genBOC(parsetname=p.parsets[parsetnames[pno]].name, progsetname=p.progsets[progsetnames[pno]].name, objectives=objectives, maxtime=maxtime)

            BOClist.append(p.getBOC(objectives))
            
        optbudgets = minBOCoutcomes(BOClist, grandtotal, budgetvec=seedbudgets, maxtime=maxtime)
            
        return optbudgets
        
        
    def fullGA(self, objectives=None, budgetratio=None, maxtime=None, doplotBOCs=False, verbose=2):
        ''' Complete geospatial analysis process applied to portfolio for a set of objectives '''
        printv('Performing full geospatial analysis', 1, verbose)
        
        GAstart = tic()

		# Check inputs
        if objectives == None: 
            printv('WARNING, you have called fullGA on portfolio %s without specifying obejctives. Using default objectives... ' % (self.name), 2, verbose)
            objectives = defaultobjectives()
        objectives = dcp(objectives)    # NOTE: Yuck. Somebody will need to check all of Optima for necessary dcps.
        
        gaoptim = GAOptim(objectives = objectives)
        self.gaoptims[str(gaoptim.uid)] = gaoptim
        
        if budgetratio == None: budgetratio = self.getdefaultbudgets()
        initbudgets = scaleratio(budgetratio,objectives['budget'])
        
        optbudgets = self.minBOCoutcomes(objectives, seedbudgets = initbudgets, maxtime = maxtime)
        if doplotBOCs: self.plotBOCs(objectives, initbudgets = initbudgets, optbudgets = optbudgets)
        
        gaoptim.complete(self.projects, initbudgets,optbudgets, maxtime=maxtime)
        outputstring = gaoptim.printresults()
        
        self.outputstring = outputstring # Store the results as an output string
        
        toc(GAstart)
        
        
        
        
#%% Functions for geospatial analysis

def constrainbudgets(x, grandtotal, minbound):
    
    # First make sure all values are not below the respective minimum bounds.
    for i in xrange(len(x)):
        if x[i] < minbound[i]:
            x[i] = minbound[i]
    
    # Then scale all excesses over the minimum bounds so that the new sum is grandtotal.
    constrainedx = []
    for i in xrange(len(x)):
        xi = (x[i] - minbound[i])*(grandtotal - sum(minbound))/(sum(x) - sum(minbound)) + minbound[i]
        constrainedx.append(xi)
    
    return constrainedx

def objectivecalc(x, BOClist, grandtotal, minbound):
    ''' Objective function. Sums outcomes from all projects corresponding to budget list x. '''
    x = constrainbudgets(x, grandtotal, minbound)
    
    totalobj = 0
    for i in xrange(len(x)):
        totalobj += BOClist[i].getoutcome([x[i]])[-1]     # Outcomes are currently passed to and from pchip as lists.
    return totalobj
    
def minBOCoutcomes(BOClist, grandtotal, budgetvec=None, minbound=None, maxiters=1000, maxtime=None, verbose=2):
    ''' Actual runs geospatial optimisation across provided BOCs. '''
    printv('Calculating minimum outcomes for grand total budget of %f' % grandtotal, 2, verbose)
    
    if minbound == None: minbound = [0]*len(BOClist)
    if budgetvec == None: budgetvec = [grandtotal/len(BOClist)]*len(BOClist)
    if not len(budgetvec) == len(BOClist): 
        errormsg = 'Geospatial analysis is minimising %i BOCs with %i initial budgets' % (len(BOClist), len(budgetvec))
        raise OptimaException(errormsg)
        
    args = {'BOClist':BOClist, 'grandtotal':grandtotal, 'minbound':minbound}    
    
#    budgetvecnew, fval, exitflag, output = asd(objectivecalc, budgetvec, args=args, xmin=budgetlower, xmax=budgethigher, timelimit=maxtime, MaxIter=maxiters, verbose=verbose)
    X, FVAL, EXITFLAG, OUTPUT = asd(objectivecalc, budgetvec, args=args, timelimit=maxtime, MaxIter=maxiters, verbose=verbose)
    X = constrainbudgets(X, grandtotal, minbound)
#    assert sum(X)==grandtotal      # Commenting out assertion for the time being, as it doesn't handle floats.

    return X



#%% Geospatial analysis runs are stored in a GAOptim object.

class GAOptim(object):
    """
    GAOPTIM

    Short for geospatial analysis optimisation. This class stores results from an optimisation run.

    Version: 2016jan26 by davidkedz
    """
        #######################################################################################################
    ## Built-in methods -- initialization, and the thing to print if you call a portfolio
    #######################################################################################################

    def __init__(self, name='default', objectives = None):
        ''' Initialize the GA optim object '''

        ## Define the structure sets
        self.objectives = objectives
        self.resultpairs = odict()

        ## Define other quantities
        self.name = name

        ## Define metadata
        self.uid = uuid()
        self.created = today()
        self.modified = today()
        self.version = __version__
        self.gitbranch, self.gitversion = gitinfo()

        return None


    def __repr__(self):
        ''' Print out useful information when called '''
        output = '============================================================\n'
        output += '      GAOptim name: %s\n'    % self.name
        output += '\n'
        output += '    Optima version: %0.1f\n' % self.version
        output += '      Date created: %s\n'    % getdate(self.created)
        output += '     Date modified: %s\n'    % getdate(self.modified)
        output += '        Git branch: %s\n'    % self.gitbranch
        output += '       Git version: %s\n'    % self.gitversion
        output += '               UID: %s\n'    % self.uid
        output += '============================================================\n'
        output += objrepr(self)
        return output
    
    
    def complete(self, projects, initbudgets, optbudgets, parsetnames=None, progsetnames=None, maxtime=None, parprogind=0, verbose=2):
        ''' Runs final optimisations for initbudgets and optbudgets so as to summarise GA optimisation '''
        printv('Finalizing geospatial analysis...', 1, verbose)
        printv('Warning, using default programset/programset!', 2, verbose)
        
        
        # Validate inputs
        if not len(projects) == len(initbudgets) or not len(projects) == len(optbudgets):
            errormsg = 'Cannot complete optimisations for %i projects given %i initial budgets (%i required) and %i optimal budgets (%i required).' % (len(self.projects), len(initbudgets), len(self.projects), len(optbudgets), len(self.projects))
            raise OptimaException(errormsg)
        if progsetnames==None:
            printv('\nWARNING: no progsets specified. Using first saved progset for each project.', 3, verbose)
            progsetnames = [0]*len(projects)
        if not len(progsetnames)==len(projects):
            printv('WARNING: %i program set names/indices were provided, but %i projects. OVERWRITING INPUTS and using first saved progset for each project.' % (len(progsetnames), len(self.projects)), 1, verbose)
            progsetnames = [0]*len(projects)
        if parsetnames==None:
            printv('\nWARNING: no parsets specified. Using first saved parset for each project.', 3, verbose)
            parsetnames = [0]*len(projects)
        if not len(parsetnames)==len(projects):
            printv('WARNING: %i parset names/indices were provided, but %i projects. OVERWRITING INPUTS and using first saved parset for each project.' % (len(parsetnames), len(self.projects)), 1, verbose)
            parsetnames = [0]*len(projects)

        # Project optimisation processes (e.g. Optims and Multiresults) are not saved to Project, only GA Optim.
        # This avoids name conflicts for Optims/Multiresults from multiple GAOptims (via project add methods) that we really don't need.
        for pind,p in enumerate(projects.values()):
            self.resultpairs[str(p.uid)] = odict()

            # Crash if any project doesn't have progsets
            if not p.progsets or not p.parsets: 
                errormsg = 'Project "%s" does not have a progset and/or a parset, can''t generate a BOC.'
                raise OptimaException(errormsg)

            initobjectives = dcp(self.objectives)
            initobjectives['budget'] = initbudgets[pind] + budgeteps
            printv("Generating initial-budget optimization for project '%s'." % p.name, 2, verbose)
            self.resultpairs[str(p.uid)]['init'] = p.optimize(name=p.name+' GA initial', parsetname=p.parsets[parsetnames[parprogind]].name, progsetname=p.progsets[progsetnames[parprogind]].name, objectives=initobjectives, maxtime=0.0, saveprocess=False) # WARNING TEMP
            preibudget = initobjectives['budget']
            postibudget = self.resultpairs[str(p.uid)]['init'].budget[-1]
#            assert abs(preibudget-sum(postibudget[:]))<tol
            
            optobjectives = dcp(self.objectives)
            optobjectives['budget'] = optbudgets[pind] + budgeteps
            printv("Generating optimal-budget optimization for project '%s'." % p.name, 2, verbose)
            self.resultpairs[str(p.uid)]['opt'] = p.optimize(name=p.name+' GA optimal', parsetname=p.parsets[parsetnames[parprogind]].name, progsetname=p.progsets[progsetnames[parprogind]].name, objectives=optobjectives, maxtime=maxtime, saveprocess=False)
            preobudget = optobjectives['budget']
            postobudget = self.resultpairs[str(p.uid)]['opt'].budget[-1]
#            assert abs(preobudget-sum(postobudget[:]))<tol

    def printresults(self, verbose=2):
        ''' Just displays results related to the GA run '''
        printv('Printing results...', 2, verbose)
        
        overallbudgetinit = 0
        overallbudgetopt = 0
        overalloutcomeinit = 0
        overalloutcomeopt = 0
        
        overalloutcomesplit = odict()
        for key in self.objectives['keys']:
            overalloutcomesplit['num'+key] = odict()
            overalloutcomesplit['num'+key]['init'] = 0
            overalloutcomesplit['num'+key]['opt'] = 0
        
        nprojects = len(self.resultpairs.keys())
        projnames = []
        projbudgets = []
        projcov = []
        projoutcomes = []
        projoutcomesplit = []
        ind = -1 # WARNING, should be a single index so doesn't actually matter
        
        for prj,x in enumerate(self.resultpairs.keys()):          # WARNING: Nervous about all this slicing. Problems foreseeable if format changes.
            # Figure out which indices to use
            tvector = self.resultpairs[x]['init'].tvec          # WARNING: NOT USING DT NORMALISATIONS LATER, SO ASSUME DT = 1 YEAR.
            initial = findinds(tvector, self.objectives['start'])
            final = findinds(tvector, self.objectives['end'])
            indices = arange(initial, final)
            
            projectname = self.resultpairs[x]['init'].project.name
            initalloc = self.resultpairs[x]['init'].budget[0]
            gaoptalloc = self.resultpairs[x]['opt'].budget[-1]
            initoutcome = self.resultpairs[x]['init'].improvement[-1][0]
            gaoptoutcome = self.resultpairs[x]['opt'].improvement[-1][-1]
            suminitalloc = sum(initalloc.values())
            sumgaoptalloc = sum(gaoptalloc.values())
            
            overallbudgetinit += suminitalloc
            overallbudgetopt += sumgaoptalloc
            overalloutcomeinit += initoutcome
            overalloutcomeopt += gaoptoutcome
            
            projnames.append(projectname)
            projbudgets.append(odict())
            projoutcomes.append(odict())
            projbudgets[prj]['init']  = initalloc
            projbudgets[prj]['opt']   = gaoptalloc
            projoutcomes[prj]['init'] = initoutcome
            projoutcomes[prj]['opt']  = gaoptoutcome
            
            projoutcomesplit.append(odict())
            projoutcomesplit[prj]['init'] = odict()
            projoutcomesplit[prj]['opt'] = odict()
            
            initpars = self.resultpairs[x]['init'].parset[0]
            optpars = self.resultpairs[x]['opt'].parset[-1]
            initprog = self.resultpairs[x]['init'].progset[0]
            optprog = self.resultpairs[x]['opt'].progset[-1]
            initcov = initprog.getprogcoverage(initalloc,self.objectives['start'],parset=initpars)
            optcov = optprog.getprogcoverage(gaoptalloc,self.objectives['start'],parset=optpars)
            
            projcov.append(odict())
            projcov[prj]['init']  = initcov
            projcov[prj]['opt']   = optcov
            
            
            for key in self.objectives['keys']:
                projoutcomesplit[prj]['init']['num'+key] = self.resultpairs[x]['init'].main['num'+key].tot[0][indices].sum()
                projoutcomesplit[prj]['opt']['num'+key] = self.resultpairs[x]['opt'].main['num'+key].tot[0][indices].sum()
                overalloutcomesplit['num'+key]['init'] += projoutcomesplit[prj]['init']['num'+key]
                overalloutcomesplit['num'+key]['opt'] += projoutcomesplit[prj]['opt']['num'+key]
                
                 
        ## Actually create the output
        output = ''
        output += 'Geospatial analysis results: minimize oucomes from %i to %i' % (self.objectives['start'], self.objectives['end'])
        output += '\n\n'
        output += '\n\t\tInitial\tOptimal'
        output += '\nOverall summary'
        output += '\n\tPortfolio budget:\t%0.0f\t%0.0f' % (overallbudgetinit, overallbudgetopt)
        output += '\n\tOutcome:\t%0.0f\t%0.0f' % (overalloutcomeinit, overalloutcomeopt)
        for key in self.objectives['keys']:
            output += '\n\t' + self.objectives['keylabels'][key] + ':\t%0.0f\t%0.0f' % (overalloutcomesplit['num'+key]['init'], overalloutcomesplit['num'+key]['opt'])
        for prj in range(nprojects):
            output += '\n'
            output += '\n'
            output += '\n\t\tInitial\tOptimal'
            output += '\nProject: "%s"' % projnames[prj]
            output += '\n'
            output += '\n\tBudget:\t%0.0f\t%0.0f' % (sum(projbudgets[prj]['init'][:]), sum(projbudgets[prj]['opt'][:]))
            output += '\n\tOutcome:\t%0.0f\t%0.0f' % (projoutcomes[prj]['init'], projoutcomes[prj]['opt'])
            for key in self.objectives['keys']:
                output += '\n\t' + self.objectives['keylabels'][key] + ':\t%0.0f\t%0.0f' % (projoutcomesplit[prj]['init']['num'+key], projoutcomesplit[prj]['opt']['num'+key])
            output += '\n'
            output += '\n\tAllocation:'
            for prg in projbudgets[prj]['init'].keys():
                output += '\n\t%s\t%0.0f\t%0.0f' % (prg, projbudgets[prj]['init'][prg], projbudgets[prj]['opt'][prg])
            output += '\n'
            output += '\n\tCoverage (%i):' % (self.objectives['start'])
            for prg in projbudgets[prj]['init'].keys():
                initval = projcov[prj]['init'][prg]
                optval = projcov[prj]['opt'][prg]
                if initval is None: initval = 0
                if optval is None: optval = 0
                output += '\n\t%s\t%0.0f\t%0.0f' % (prg, initval, optval)
        
        print(output)
        
        return output
        
    def getinitbudgets(self):
        bl = []
        for proj in self.resultpairs:
            bl.append(sum(self.resultpairs[proj]['init'].budget[0].values()))
        return bl
        
    def getoptbudgets(self):
        bl = []
        for proj in self.resultpairs:
            bl.append(sum(self.resultpairs[proj]['opt'].budget[-1].values()))
        return bl
        
    
#%% 'EASY' STACK-PLOTTING CODE PULLED FROM v1.5 FOR CLIFF. DOES NOT WORK.   :)
    

#    def sortfixed(somelist):
#        # Prog array help: [VMMC, FSW, MSM, HTC, ART, PMTCT, OVC, Other Care, MGMT, HR, ENV, SP, M&E, Other, SBCC, CT]
#        sortingarray = [2, 4, 5, 7, 8, 6, 1, 9, 10, 11, 12, 13, 14, 15, 3, 0]
#        return [y for x,y in sorted(zip(sortingarray,somelist), key = lambda x:x[0])]
#    #    return [y for x,y in sorted(enumerate(somelist), key = lambda x: -len(progs[x[0]]['effects']))]
        
    
    def superplot(self):
        
        from matplotlib.pylab import gca, xlabel, tick_params, xlim, figure, subplot, plot, pie, bar, title, legend, xticks, ylabel, show
        from gridcolormap import gridcolormap
        from matplotlib import gridspec
        import numpy
        
        progs = p1.regionlist[0].metadata['programs']    
        
        figure(figsize=(22,15))
        
        nprograms = len(p1.gpalist[-1][0].region.data['origalloc'])
#        colors = sortfixed(gridcolormap(nprograms))
#        colors[0] = numpy.array([ 0.20833333,  0.20833333,  0.54166667])  #CT
#        colors[1] = numpy.array([ 0.45833333,  0.875     ,  0.79166667])  #OVC
#        colors[2] = numpy.array([0.125, 0.125, 0.125])  #VMMC
#        colors[3] = numpy.array([ 0.79166667,  0.45833333,  0.875     ])+numpy.array([ 0.125,  0.125,  0.125     ])  #SBCC
#        colors[4] = numpy.array([ 0.54166667,  0.20833333,  0.20833333])  #FSW
#        colors[5] = numpy.array([ 0.875     ,  0.45833333,  0.125     ])  #MSM
#        colors[6] = numpy.array([ 0.125     ,  0.875     ,  0.45833333])+numpy.array([ 0.0,  0.125,  0.0     ])  #PMTCT
#        colors[7] = numpy.array([ 0.54166667,  0.875     ,  0.125     ])   #HTC
#        colors[8] = numpy.array([ 0.20833333,  0.54166667,  0.20833333])  #ART
#        for i in xrange(7): colors[-(i+1)] = numpy.array([0.25+0.5*i/6.0, 0.25+0.5*i/6.0, 0.25+0.5*i/6.0])
    
        
        
        gpl = sorted(p1.gpalist[-1], key=lambda sbo: sbo.name)
        ind = [val for pair in ([x, 0.25+x] for x in xrange(len(gpl))) for val in pair]
        width = [0.25, 0.55]*(len(gpl))       # the width of the bars: can also be len(x) sequence
        
        bar(ind, [val*1e-6 for pair in zip([sortfixed(sb.simlist[1].alloc)[-1] for sb in gpl], [sortfixed(sb.simlist[2].alloc)[-1] for sb in gpl]) for val in pair], width, color=colors[-1])
        for p in xrange(2,nprograms+1):
            bar(ind, [val*1e-6 for pair in zip([sortfixed(sb.simlist[1].alloc)[-p] for sb in gpl], [sortfixed(sb.simlist[2].alloc)[-p] for sb in gpl]) for val in pair], width, color=colors[-p], bottom=[val*1e-6 for pair in zip([sum(sortfixed(sb.simlist[1].alloc)[1-p:]) for sb in gpl], [sum(sortfixed(sb.simlist[2].alloc)[1-p:]) for sb in gpl]) for val in pair])
        
        xticks([x+0.5 for x in xrange(len(gpl))], [sb.region.getregionname() for sb in gpl], rotation=-60)
        xlim([0,32])
        tick_params(axis='both', which='major', labelsize=15)
        tick_params(axis='both', which='minor', labelsize=15)
        ylabel('Budget Allocation (US$m)', fontsize=15)
        
    
    
        fig = figure(figsize=(22,15))    
        
        gs = gridspec.GridSpec(3, 11) #, width_ratios=[len(sb.simlist[1:]), 2])
        
        for x in xrange(len(gpl)):
            sb = gpl[x]
            r = sb.region
    
            ind = xrange(len(sb.simlist[1:]))
            width = 0.8       # the width of the bars: can also be len(x) sequence
            
            if x < 10: subplot(gs[x])
            else: subplot(gs[x+1])
            bar(ind, [sortfixed([x*1e-6 for x in sim.alloc])[-1] for sim in sb.simlist[1:]], width, color=colors[-1])
            for p in xrange(2,nprograms+1):
                bar(ind, [sortfixed([x*1e-6 for x in sim.alloc])[-p] for sim in sb.simlist[1:]], width, color=colors[-p], bottom=[sum(sortfixed([x*1e-6 for x in sim.alloc])[1-p:]) for sim in sb.simlist[1:]])
            #xticks([index+width/2.0 for index in ind], [sim.getname() for sim in sb.simlist[1:]])
            xlabel(r.getregionname(), fontsize=18)
            if x in [0,10,21]: ylabel('Budget Allocation (US$m)', fontsize=18)
            tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
        
    #        ax = gca()
    #        ax.ticklabel_format(style='sci', axis='y')
    #        ax.yaxis.major.formatter.set_powerlimits((0,0))
            
            tick_params(axis='both', which='major', labelsize=14)
            tick_params(axis='both', which='minor', labelsize=14)
        
        fig.tight_layout()
        
        
        
    #    bar(ind, [val for pair in zip([sb.simlist[1].alloc[-1] for sb in gpl], [sb.simlist[2].alloc[-1] for sb in gpl]) for val in pair], width, color=colors[-1])
    #    for p in xrange(2,nprograms+1):
    #        bar(ind, [val for pair in zip([sb.simlist[1].alloc[-p] for sb in gpl], [sb.simlist[2].alloc[-p] for sb in gpl]) for val in pair], width, color=colors[-p], bottom=[val for pair in zip([sum(sb.simlist[1].alloc[1-p:]) for sb in gpl], [sum(sb.simlist[2].alloc[1-p:]) for sb in gpl]) for val in pair])
    #    
    #    xticks([x+0.5 for x in xrange(len(gpl))], [sb.region.getregionname() for sb in gpl], rotation=-60)
    #    xlim([0,32])
    #    tick_params(axis='both', which='major', labelsize=15)
    #    tick_params(axis='both', which='minor', labelsize=15)
    #    ylabel('Budget Allocation ($)', fontsize=15)
    #    
    #
    #
    #    fig = figure(figsize=(22,15))    
    #    
    #    gs = gridspec.GridSpec(3, 11) #, width_ratios=[len(sb.simlist[1:]), 2])
    #    
    #    for x in xrange(len(gpl)):
    #        sb = gpl[x]
    #        r = sb.region
    #
    #        ind = xrange(len(sb.simlist[1:]))
    #        width = 0.8       # the width of the bars: can also be len(x) sequence
    #        
    #        if x < 10: subplot(gs[x])
    #        else: subplot(gs[x+1])
    #        bar(ind, [sim.alloc[-1] for sim in sb.simlist[1:]], width, color=colors[-1])
    #        for p in xrange(2,nprograms+1):
    #            bar(ind, [sim.alloc[-p] for sim in sb.simlist[1:]], width, color=colors[-p], bottom=[sum(sim.alloc[1-p:]) for sim in sb.simlist[1:]])
    #        #xticks([index+width/2.0 for index in ind], [sim.getname() for sim in sb.simlist[1:]])
    #        xlabel(r.getregionname(), fontsize=18)
    #        if x in [0,10,21]: ylabel('Budget Allocation ($)', fontsize=18)
    #        tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    #    
    #        ax = gca()
    #        ax.ticklabel_format(style='sci', axis='y')
    #        ax.yaxis.major.formatter.set_powerlimits((0,0))
    #        
    #        tick_params(axis='both', which='major', labelsize=13)
    #        tick_params(axis='both', which='minor', labelsize=13)
    #    
    #    fig.tight_layout()    
        
        
        
    #    for x in xrange(len(p1.gpalist[-1])):
        for x in xrange(1):
            sb = p1.gpalist[-1][x]
            r = sb.region
            
            nprograms = len(r.data['origalloc'])
    #        colors = sortfixed(gridcolormap(nprograms))
            
            figure(figsize=(len(sb.simlist[1:])*2+4,nprograms/2))
            gs = gridspec.GridSpec(1, 2, width_ratios=[len(sb.simlist[1:]), 2]) 
            ind = xrange(len(sb.simlist[1:]))
            width = 0.8       # the width of the bars: can also be len(x) sequence
            
            subplot(gs[0])
            bar(ind, [sortfixed(sim.alloc)[-1] for sim in sb.simlist[1:]], width, color=colors[-1])
            for p in xrange(2,nprograms+1):
                bar(ind, [sortfixed(sim.alloc)[-p] for sim in sb.simlist[1:]], width, color=colors[-p], bottom=[sum(sortfixed(sim.alloc)[1-p:]) for sim in sb.simlist[1:]])
            xticks([index+width/2.0 for index in ind], [sim.getname() for sim in sb.simlist[1:]])
            ylabel('Budget Allocation ($)')
            
            subplot(gs[1])
            for prog in xrange(nprograms): plot(0, 0, linewidth=3, color=colors[prog])
            legend(sortfixed(r.data['meta']['progs']['short']))
            
        show()
                
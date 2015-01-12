def loadworkbook(filename='example.xlsx', verbose=2):
    """
    Loads the workbook (i.e. reads its contents into the data structure).
    This data structure is used in the next step to update the corresponding model.
    The workbook is assumed to be in the format specified in example.xlsx.
    
    Version: 2015jan13
    """
    

    ###########################################################################
    ## Preliminaries
    ###########################################################################
    
    from printv import printv
    from numpy import nan, array # For reading in empty values
    from xlrd import open_workbook # For opening Excel workbooks
    from bunch import Bunch as struct # Replicate Matlab-like structure behavior
    from time import strftime # For determining when a spreadsheet was last uploaded
    printv('Loading data from %s...' % filename, 1, verbose)
    from programs import programs_for_input_key
#    from makeccocs import default_init_coparams

        
    ###########################################################################
    ## Define the workbook and parameter names
    ###########################################################################
    
    # Metadata -- population and program names -- array sizes are (# populations) and (# programs)
    # groupname   sheetname                 name    thispar
    metadata = [['Populations & programs', 'meta', ['pops', 'progs']]]
    
    # Key data -- array sizes are time x population x uncertainty
    keydata =  [['Demographics & HIV prevalence', 'key', ['popsize', 'hivprev']]]
    
    # Cost-coverage data -- array sizes are time x programs x cost/coverage
    cocodata = [['Cost & coverage',     'costcov', ['cov', 'cost']]]
    
    # Time data -- array sizes are time x population
    timedata = [
                 ['Other epidemiology',  'epi',     ['death', 'stiprevulc', 'tbprev']], \
                 ['Optional indicators', 'opt',     ['numtest', 'numdiag', 'numinfect', 'prev', 'death', 'newtreat']], \
                 ['Testing & treatment', 'txrx',    ['hivtest', 'aidstest', 'numfirstline', 'numsecondline', 'txelig', 'prep', 'numpmtct', 'birth', 'breast']], \
                 ['Sexual behavior',     'sex',     ['numactsreg', 'numactscas', 'numactscom', 'condomreg', 'condomcas', 'condomcom', 'circum', 'numcircum']], \
                 ['Injecting behavior',  'inj',     ['numinject', 'sharing', 'numost']], \
                 ['Macroeconomics & costs',      'macro',   ['gdp', 'revenue', 'govtexpend', 'totalhealth', 'domestichealth', 'domestichiv', 'globalfund', 'pepfar', 'otherint', 'private']]
                ]
                 
    # Matrix data -- array sizes are population x population
    matrices = [
                ['Partnerships', 'pships',  ['reg','cas','com','inj']], \
                ['Transitions',  'transit', ['asym','sym']]
               ]
    
    # Constants -- array sizes are scalars x uncertainty
    constants = [
                 ['Constants', 'const',              [['trans',    ['mfi', 'mfr', 'mmi', 'mmr', 'inj', 'mtctbreast', 'mtctnobreast']], \
                                                      ['cd4trans', ['acute', 'gt500', 'gt350', 'gt200', 'aids']], \
                                                      ['prog',     ['acute', 'gt500', 'gt350', 'gt200']],\
                                                      ['recov',    ['gt500', 'gt350', 'gt200', 'aids']],\
                                                      ['fail',     ['first', 'second']],\
                                                      ['death',    ['acute', 'gt500', 'gt350', 'gt200', 'aids', 'treat', 'tb']],\
                                                      ['eff',      ['condom', 'circ', 'dx', 'sti', 'meth', 'pmtct', 'tx', 'prep', 'pep']]]], \
                 ['Disutilities & costs', 'cost',    [['disutil',  ['acute', 'gt500', 'gt350', 'gt200', 'aids','tx']], \
                                                      ['health',   ['acute', 'gt500', 'gt350', 'gt200', 'aids']], \
                                                      ['social',   ['acute', 'gt500', 'gt350', 'gt200', 'aids']]]]
                ]
    
    
    ## Ugly, but allow the list of groups to be used as name and also as variables
    sheetstructure = struct()
    sheetstructure.metadata = metadata
    sheetstructure.cocodata = cocodata
    sheetstructure.keydata = keydata
    sheetstructure.timedata = timedata
    sheetstructure.matrices = matrices
    sheetstructure.constants = constants
    


    ###########################################################################
    ## Load data sheets
    ###########################################################################
    

    ## Basic setup
    data = struct() # Create structure for holding data
    data.__doc__ = 'Raw data as loaded from the workbook, including both epidemiological and behavioral data, plus economics and velociraptors.'
    data.__date__ = strftime("%Y-%m-%d %H:%M:%S")
    programs = struct() # Create structure for holding program data
    programs.__doc__ = 'Parameters that define the HIV programs -- cost-coverage and coverage-outcome curves.'
    workbook = open_workbook(filename) # Open workbook
    
    sheetstructure_keys = sheetstructure.keys()
    metadata_index = sheetstructure_keys.index('metadata')
    #ensure that metadata is parsed first
    sheetstructure_keys = ['metadata']+ sheetstructure_keys[:metadata_index]+sheetstructure_keys[metadata_index+1:]
    
    ## Loop over each group of sheets
    for groupname in sheetstructure_keys: # Loop over each type of data, but treat constants differently
        sheetgroup = sheetstructure[groupname]
        for sheet in sheetgroup: # Loop over each workbook for that data -- just one for constants
            lastdatacol = None
            sheetname = sheet[0] # Name of the workbook
            name = sheet[1] # Pull out the name of this field, e.g. 'epi'
            subparlist = sheet[2] # List of subparameters
            data[name] = struct() # Create structure for holding data, e.g. data.epi
            sheetdata = workbook.sheet_by_name(sheetname) # Load this workbook
            parcount = -1 # Initialize the parameter count
            printv('  Loading "%s"...' % sheetname, 2, verbose)
            
            
            ## Calculate columns for which data are entered, and store the year ranges
            if groupname in ['keydata', 'cocodata', 'timedata']  and name != 'macro': # Need to gather year ranges for epidemic etc. data
                data.epiyears = [] # Initialize epidemiology data years
                for col in range(sheetdata.ncols):
                    thiscell = sheetdata.cell_value(1,col) # 1 is the 2nd row which is where the year data should be
                    if thiscell=='' and len(data.epiyears)>0: #  We've gotten to the end
                        lastdatacol = col # Store this column number
                        break # Quit
                    elif thiscell != '': # Nope, more years, keep going
                        data.epiyears.append(float(thiscell)) # Add this year
            
            if name == 'macro': # Need to gather year ranges for economic data
                data.econyears = [] # Initialize epidemiology data years
                for col in range(sheetdata.ncols):
                    thiscell = sheetdata.cell_value(1,col) # 1 is the 2nd row which is where the year data should be
                    if thiscell=='' and len(data.econyears)>0: #  We've gotten to the end
                        lastdatacol = col # Store this column number
                        break # Quit
                    elif thiscell != '': # Nope, more years, keep going
                        data.econyears.append(float(thiscell)) # Add this year
            
            if lastdatacol:    
                assumptioncol = lastdatacol + 1 # The "OR" space is in between
                ncolsperprog = 5 # Number of columns necessary for defining a single program; name, zero-spend-min, zero-spend-max, full-spend-min, full-spend-max
                nprogblocks = 4 # Number of program blocks
                programcols = assumptioncol + 3 + array([array(range(ncolsperprog))+(1+ncolsperprog)*i for i in range(nprogblocks)]) # Calculate which columns the program data is stored in
            
            
            
            
            ##################################################################
            ## Now, actually load the data
            ##################################################################
            
            
            # Loop over each row in the workbook
            for row in range(sheetdata.nrows): 
                paramcategory = sheetdata.cell_value(row,0) # See what's in the first column for this row
                
                if paramcategory != '': # It's not blank: e.g. "HIV prevalence"
                    printv('Loading "%s"...' % paramcategory, 3, verbose)
                    parcount += 1 # Increment the parameter count
                    
                    # It's metadata: pull out each of the pieces
                    if groupname=='metadata': 
                        thispar = subparlist[parcount] # Get the name of this parameter, e.g. 'pop'
                        data[name][thispar] = struct() # Initialize to empty list
                        data[name][thispar].short = [] # Store short population/program names, e.g. "FSW"
                        data[name][thispar].long = [] # Store long population/program names, e.g. "Female sex workers"
                        if thispar=='pops':
                            data[name][thispar].male = [] # Store whether or not this population is male
                            data[name][thispar].female = [] # Store whether or not this population is female
                            data[name][thispar].injects = [] # Store whether or not this population injects drugs
                            data[name][thispar].sexmen = [] # Store whether or not this population has sex with men
                            data[name][thispar].sexwomen = [] # Store whether or not this population has sex with women
                            data[name][thispar].sexworker = [] # Store whether or not this population is a sex worker
                            data[name][thispar].client = [] # Store whether or not this population is a client of sex workers
                    
                    # It's cost-coverage data: store cost and coverage for each program
                    elif groupname=='cocodata': 
                        data[name][subparlist[0]] = [] # Initialize coverage to an empty list -- i.e. data.costcov.cov
                        data[name][subparlist[1]] = [] # Initialize cost to an empty list -- i.e. data.costcov.cost
                        
                    elif groupname in ['keydata', 'timedata', 'matrices']: # It's basic data or a matrix: create an empty list
                        thispar = subparlist[parcount] # Get the name of this parameter, e.g. 'popsize'                    
                        data[name][thispar] = [] # Initialize to empty list
    
                    elif groupname=='constants': # It's a constant or a cost: create a structure
                        thispar = subparlist[parcount][0] # Get the name of this parameter, e.g. 'trans'
                        data[name][thispar] = struct() # Need yet another structure if it's a constant!
                    
                    else:
                        raise Exception('Group name %s not recognized!' % groupname)
                        

                
                if paramcategory == '': # The first column is blank: it's time for the data
                    subparam = sheetdata.cell_value(row, 1) # Get the name of a subparameter, e.g. 'FSW', population size for a given population
                    
                    if subparam != '': # The subparameter name isn't blank, load something!
                        printv('Parameter: %s' % subparam, 4, verbose)
                        
                        # It's meta-data, split into pieces
                        if groupname=='metadata': 
                            thesedata = sheetdata.row_values(row, start_colx=2, end_colx=11) # Data starts in 3rd column, finishes in 11th column
                            data[name][thispar].short.append(thesedata[0])

                            data[name][thispar].long.append(thesedata[1])
                            if thispar=='pops':
                                data[name][thispar].male.append(thesedata[2])
                                data[name][thispar].female.append(thesedata[3])
                                data[name][thispar].injects.append(thesedata[4])
                                data[name][thispar].sexmen.append(thesedata[5])
                                data[name][thispar].sexwomen.append(thesedata[6])
                                data[name][thispar].sexworker.append(thesedata[7])
                                data[name][thispar].client.append(thesedata[8])
                            if thispar=='progs':
                                if not thesedata[0] in programs: programs[thesedata[0]] = []
                                
                        # It's cost-coverage data, save the cost and coverage values separately
                        if groupname=='cocodata':
                            thesedata = sheetdata.row_values(row, start_colx=3, end_colx=lastdatacol) # Data starts in 4th column
                            thesedata = map(lambda val: nan if val=='' else val, thesedata) # Replace blanks with nan
                            assumptiondata = sheetdata.cell_value(row, assumptioncol)
                            if assumptiondata != '': thesedata = [assumptiondata] # Replace the (presumably blank) data if a non-blank assumption has been entered
                            ccindices = {'Coverage':0, 'Cost':1} # Define best-low-high indices
                            cc = sheetdata.cell_value(row, 2) # Read in whether indicator is best, low, or high
                            data[name][subparlist[ccindices[cc]]].append(thesedata) # Actually append the data
                        
                        
                        # It's key data, save both the values and uncertainties
                        if groupname=='keydata':
                            if len(data[name][thispar])==0: 
                                data[name][thispar] = [[] for z in range(3)] # Create new variable for best, low, high
                            thesedata = sheetdata.row_values(row, start_colx=3, end_colx=lastdatacol) # Data starts in 4th column
                            thesedata = map(lambda val: nan if val=='' else val, thesedata) # Replace blanks with nan
                            assumptiondata = sheetdata.cell_value(row, assumptioncol)
                            if assumptiondata != '': thesedata = [assumptiondata] # Replace the (presumably blank) data if a non-blank assumption has been entered
                            blhindices = {'best':0, 'low':1, 'high':2} # Define best-low-high indices
                            blh = sheetdata.cell_value(row, 2) # Read in whether indicator is best, low, or high
                            data[name][thispar][blhindices[blh]].append(thesedata) # Actually append the data
                            
                        
                        # It's basic data, append the data and check for programs
                        if groupname=='timedata': 
                            thesedata = sheetdata.row_values(row, start_colx=2, end_colx=lastdatacol) # Data starts in 3rd column
                            thesedata = map(lambda val: nan if val=='' else val, thesedata) # Replace blanks with nan
                            assumptiondata = sheetdata.cell_value(row, assumptioncol)
                            if assumptiondata != '': thesedata = [assumptiondata] # Replace the (presumably blank) data if a non-blank assumption has been entered
                            data[name][thispar].append(thesedata) # Store data

                            for programname in programs_for_input_key(thispar):
                                if programname in programs:
#                                    programs[programname].append([[name, thispar], [subparam], default_init_coparams])                        
                                    programs[programname].append([[name, thispar], [subparam]])
                        
                        # It's a matrix, append the data                                     
                        elif groupname=='matrices':
                            thesedata = sheetdata.row_values(row, start_colx=2, end_colx=sheetdata.ncols) # Data starts in 3rd column
                            thesedata = map(lambda val: 0 if val=='' else val, thesedata) # Replace blanks with 0
                            data[name][thispar].append(thesedata) # Store data
                        
                        # It's a constant, create a new dictionary entry
                        elif name=='const' or name=='cost':
                            thesedata = sheetdata.row_values(row, start_colx=2, end_colx=5) # Data starts in 3rd column, finishes in 5th column
                            thesedata = map(lambda val: nan if val=='' else val, thesedata) # Replace blanks with nan
                            subpar = subparlist[parcount][1].pop(0) # Pop first entry of subparameter list, which is namelist[parcount][1]
                            data[name][thispar][subpar] = thesedata # Store data
    
    printv('...done loading data.', 2, verbose)
    return data, programs

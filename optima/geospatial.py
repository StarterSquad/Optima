"""
GEOSPATIAL

This file defines everything needed for the Python GUI for geospatial analysis.

Version: 2016jan23
"""

from optima import Project, Portfolio, loadobj, saveobj, odict, defaultobjectives
from PyQt4 import QtGui
from pylab import figure, close
global geoguiwindow
geoguiwindow = None

def geogui():
    '''
    Open the GUI for doing geospatial analysis.
    
    Version: 2016jan23
    '''
    global geoguiwindow, portfolio, projectslist, objectives
    portfolio = None
    projectslist = []
    objectives = defaultobjectives()
    
    ## Set parameters
    wid = 650.0
    hei = 550.0
    top = 20
    spacing = 40
    left = 20.
    projext = '.prj'
    portext = '.prt'
    
    ## Housekeeping
    fig = figure(); close(fig) # Open and close figure...dumb, no? Otherwise get "QWidget: Must construct a QApplication before a QPaintDevice"
    geoguiwindow = QtGui.QWidget() # Create panel widget
    geoguiwindow.setGeometry(100, 100, wid, hei)
    geoguiwindow.setWindowTitle('Optima geospatial analysis')
    projectslist = []
    
    ##############################################################################################################################
    ## Define functions
    ##############################################################################################################################

    def _checkproj(filepath):
        ''' Check that the filename provided is a project; if so, return it; else, return None '''
        project = None
        try: 
            project = loadobj(filepath, verbose=0)
        except: 
            print('Could not load file "%s"' % filepath)
            return None
        try: 
            assert type(project)==Project
            return project
        except: 
            print('File "%s" is not an Optima project file' % filepath)
            return None
    
    
    def _loadproj():
        ''' Little helper function to load a project, since used more than once '''
        filepath = QtGui.QFileDialog.getOpenFileNames(caption='Choose project file', filter='*'+projext)
        return _checkproj(filepath)
        
        
    def makesheet():
        ''' Create a geospatial spreadsheet template based on a project file '''
        
        ## 1. Load a project file
        project = _loadproj()
            
        ## 2. Get destination filename
        spreadsheetpath = QtGui.QFileDialog.getSaveFileName(caption='Save geospatial spreadsheet file', filter='*.xlsx')
        
        ## 3. Extract data needed from project (population names, program names...)
        # ...
        
        ## 4. Generate and save spreadsheet
        # ...

        return None
        
    
    def makeproj():
        ''' Create a series of project files based on a seed file and a geospatial spreadsheet '''
        
        ## 1. Load a project file -- WARNING, could be combined with the above!
        project = _loadproj()
        
        return None


    def create():
        ''' Create a portfolio by selecting a list of projects; silently skip files that fail '''
        global projectslist
        projectslist = []
        projectpaths = []
        filepaths = QtGui.QFileDialog.getOpenFileNames(caption='Choose project files', filter='*'+projext)
        for filepath in filepaths:
            tmpproj = None
            try: tmpproj = loadobj(filepath, verbose=0)
            except: print('Could not load file "%s"; moving on...' % filepath)
            if tmpproj is not None: 
                try: 
                    assert type(tmpproj)==Project
                    projectslist.append(tmpproj)
                    projectpaths.append(filepath)
                    print('Project file "%s" loaded' % filepath)
                except: print('File "%s" is not an Optima project file; moving on...' % filepath)
        projectsbox.setText('\n'.join(projectpaths))
        portfolio = Portfolio()
        for project in projectslist: portfolio.addproject(project)
        return None
    
    def loadport():
        ''' Load an existing portfolio '''
        global portfolio
        filepath = QtGui.QFileDialog.getOpenFileName(caption='Choose portfolio file', filter='*'+portext)
        tmpport = None
        try: tmpport = loadobj(filepath, verbose=0)
        except: print('Could not load file "%s"' % filepath)
        if tmpport is not None: 
            try: 
                assert type(tmpport)==Portfolio
                portfolio = tmpport
                print('Portfolio file "%s" loaded' % filepath)
            except: print('File "%s" is not an Optima portfolio file' % filepath)
        projectsbox.setText('\n'.join([proj.name for proj in portfolio.projects.values()]))
        portfolio = Portfolio()
        for project in projectslist: portfolio.addproject(project)
        return None
    
    
    def rungeo():
        ''' Actually run geospatial analysis!!! '''
        global portfolio, objectives
        portfolio.genBOCs(objectives)
        portfolio.plotBOCs(objectives)
        # ...
        return None
    
    
    def export():
        ''' Save the current portfolio to disk '''
        global portfolio
        if type(portfolio)!=Portfolio: print('Warning, must load portfolio first!')
        
        # 1. Extract data needed from portfolio
        # ...
        
        # 2. Generate spreadsheet according to David's template to store these data
        # ...
        
        # 3. Create a new file dialog to save this spreadsheet
        # ...
        return None
        

    def saveport():
        ''' Save the current portfolio '''
        global portfolio
        filepath = QtGui.QFileDialog.getSaveFileName(caption='Save portfolio file', filter='*'+portext)
        saveobj(filepath, portfolio)
        return None


    def closewindow(): 
        ''' Close the control panel '''
        global geoguiwindow
        geoguiwindow.close()
        return None
    
    
    ##############################################################################################################################
    ## Define buttons
    ##############################################################################################################################
    buttons = odict()
    buttons['makesheet'] = QtGui.QPushButton('Make geospatial template from project', parent=geoguiwindow)
    buttons['makeproj']  = QtGui.QPushButton('Auto-generate projects template', parent=geoguiwindow)
    buttons['create']    = QtGui.QPushButton('Create portfolio from projects', parent=geoguiwindow)
    buttons['loadport']  = QtGui.QPushButton('Load existing portfolio', parent=geoguiwindow)
    buttons['rungeo']    = QtGui.QPushButton('Run geospatial analysis', parent=geoguiwindow)
    buttons['export']    = QtGui.QPushButton('Export results', parent=geoguiwindow)
    buttons['saveport']  = QtGui.QPushButton('Save portfolio', parent=geoguiwindow)
    buttons['close']     = QtGui.QPushButton('Close', parent=geoguiwindow)
    
    actions = odict()
    actions['makesheet'] = makesheet
    actions['makeproj']  = makeproj
    actions['create']    = create
    actions['loadport']  = loadport
    actions['rungeo']    = rungeo
    actions['export']    = export
    actions['saveport']  = saveport
    actions['close']     = closewindow
    
    ## Set button locations
    for b,button in enumerate(buttons.values()):
        button.move(left, top+spacing*b)
    
    ## Define button functions
    for k in buttons.keys():
        buttons[k].clicked.connect(actions[k])
    
    
    
    ## Define other objects
    projectsbox = QtGui.QTextEdit(parent=geoguiwindow)
    projectsbox.move(300, 20)
    projectsbox.resize(wid-320, hei-40)
    

    geoguiwindow.show()
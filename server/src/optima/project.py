import json
from flask import Blueprint, url_for, helpers, request, jsonify, redirect, current_app, Response
from werkzeug.utils import secure_filename
import os
import traceback
from sim.dataio import upload_dir_user, DATADIR, TEMPLATEDIR, fullpath
from sim.updatedata import updatedata
from sim.makeproject import makeproject, makeworkbook
from utils import allowed_file, project_exists, delete_spreadsheet, load_project
from utils import check_project_name, load_model, save_model, report_exception, model_as_bunch, model_as_dict
from flask.ext.login import login_required, current_user
from dbconn import db
from dbmodels import ProjectDb, WorkingProjectDb, ProjectDataDb
from utils import BAD_REPLY
import time,datetime
import dateutil.tz
from datetime import datetime


""" route prefix: /api/project """
project = Blueprint('project',  __name__, static_folder = '../static')
project.config = {}

@project.record
def record_params(setup_state):
  app = setup_state.app
  project.config = dict([(key,value) for (key,value) in app.config.iteritems()])

@project.route('/params')
@login_required
def get_project_params():
    """
    Gives back project params
    """
    from sim.parameters import parameters
    project_params = [p for p in parameters() if p['modifiable']]
    return json.dumps({"params":project_params})

@project.route('/predefined')
@login_required
def get_predefined():
    """
    Gives back default populations and programs
    """
    from sim.programs import programs
    from sim.populations import populations
    from sim.program_categories import program_categories
    programs = programs()
    populations = populations()
    program_categories = program_categories()
    category_per_program = {}
    for category in program_categories:
        for p in category['programs']:
            category_per_program[p['short_name']] = category['category']
    for p in populations: p['active']= False
    for p in programs:
        p['active'] = False
        p['category'] = category_per_program[p['short_name']]
        new_params = [dict([('value', param),('active',True)]) for param in p['parameters']]
        for np in new_params:
            if len(np['value']['pops'][0])==0: np['value']['pops']=['ALL_POPULATIONS']
        if new_params: p['parameters'] = new_params
    return json.dumps({"programs":programs, "populations": populations, "categories":program_categories})

@project.route('/create/<project_name>', methods=['POST'])
@login_required
@report_exception()
def createProject(project_name):
    """
    Creates the project with the given name and provided parameters.
    Result: on the backend, new project is stored,
    spreadsheet with specified name and parameters given back to the user.
    expects json with the following arguments (see example):
    {"npops":6,"nprogs":8, "datastart":2000, "dataend":2015}
    """
    from sim.makeproject import default_datastart, default_dataend, default_econ_dataend, default_pops, default_progs

    current_app.logger.debug("createProject %s" % project_name)
    data = request.form

    # get current user
    user_id = current_user.id
    edit_params = None
    if data:
        edit_params = json.loads(data['edit_params'])
        data = json.loads(data['params'])

    # check if current request is edit request
    is_edit = edit_params and edit_params.get('isEdit')
    can_update = edit_params and edit_params.get('canUpdate')

    makeproject_args = {"projectname":project_name, "savetofile":False}
    makeproject_args['datastart'] = data.get('datastart', default_datastart)
    makeproject_args['dataend'] = data.get('dataend', default_dataend)
    makeproject_args['econ_dataend'] = data.get('econ_dataend', default_econ_dataend)
    makeproject_args['progs'] = data.get('programs', default_progs)
    makeproject_args['pops'] = data.get('populations', default_pops)
    current_app.logger.debug("createProject(%s)" % makeproject_args)

    # See if there is matching project
    project = load_project(project_name)

    # update existing
    if project is not None:
        # set new project name if not none

        project.datastart = makeproject_args['datastart']
        project.dataend = makeproject_args['dataend']
        project.econ_dataend = makeproject_args['econ_dataend']
        project.programs = makeproject_args['progs']
        project.populations = makeproject_args['pops']
        current_app.logger.debug('Updating existing project %s' % project.name)
    else:
        user_id = current_user.id
        # create new project
        project = ProjectDb(project_name, user_id, makeproject_args['datastart'], makeproject_args['dataend'], \
            makeproject_args['econ_dataend'], makeproject_args['progs'], makeproject_args['pops'])
        current_app.logger.debug('Creating new project: %s' % project.name)

    D = makeproject(**makeproject_args) # makeproject is supposed to return the name of the existing file...
    project.model = D.toDict()
    if is_edit:
        db.session.query(WorkingProjectDb).filter_by(id=project.id).delete()
        if can_update and project.project_data is not None and project.project_data.meta is not None:
            # try to reload the data
            loaddir =  upload_dir_user(DATADIR)
            if not loaddir:
                loaddir = DATADIR
            filename = project_name + '.xlsx'
            server_filename = os.path.join(loaddir, filename)
            filedata = open(server_filename, 'wb')
            filedata.write(project.project_data.meta)
            filedata.close()
            D = model_as_bunch(project.model)
            D = updatedata(D, savetofile = False)
            model = model_as_dict(D)
            project.model = model
        else:
            db.session.query(ProjectDataDb).filter_by(id=project.id).delete()

    # Save to db
    db.session.add(project)
    db.session.commit()
    new_project_template = D.G.workbookname

    current_app.logger.debug("new_project_template: %s" % new_project_template)
    (dirname, basename) = (upload_dir_user(TEMPLATEDIR), new_project_template)
    return helpers.send_from_directory(dirname, basename)

@project.route('/open/<project_name>')
@login_required
def openProject(project_name):
    """
    Opens the project with the given name.
    If the project exists, notifies the user about success.
    expects project name,
    todo: only if it can be found
    """
    proj_exists = False
    try: #first check DB
        proj_exists = project_exists(project_name)
        current_app.logger.debug("proj_exists: %s" % proj_exists)
    except:
        proj_exists = False
    if not proj_exists:
        return jsonify({'status':'NOK','reason':'No such project %s' % project_name})
    else:
        return jsonify({'status':'OK'})

@project.route('/workbook/<project_name>')
@login_required
@report_exception()
def giveWorkbook(project_name):
    """
    Generates workbook for the project with the given name.
    expects project name (project should already exist)
    if project exists, regenerates workbook for it
    if project does not exist, returns an error.
    """
    reply = BAD_REPLY
    proj_exists = False
    cu = current_user
    current_app.logger.debug("giveWorkbook(%s %s)" % (cu.id, project_name))
    project = load_project(project_name)
    if project is None:
        reply['reason']='Project %s does not exist.' % project_name
        return jsonify(reply)
    else:
        # See if there is matching project data
        projdata = ProjectDataDb.query.get(project.id)

        if projdata is not None and len(projdata.meta)>0:
            return Response(projdata.meta,
                mimetype= 'application/octet-stream',
                headers={'Content-Disposition':'attachment;filename='+ project_name+'.xlsx'})
        else:
        # if no project data found
            D = project.model
            wb_name = D['G']['workbookname']
            makeworkbook(wb_name, project.populations, project.programs, \
                project.datastart, project.dataend, project.econ_dataend)
            current_app.logger.debug("project %s template created: %s" % (project.name, wb_name))
            (dirname, basename) = (upload_dir_user(TEMPLATEDIR), wb_name)
            #deliberately don't save the template as uploaded data
            return helpers.send_from_directory(dirname, basename)

@project.route('/info')
@login_required
@check_project_name
def getProjectInformation():
    """
    Returns information of the requested project. (Including status of the model)

    Returns:
        A jsonified project dictionary accessible to the current user.
        In case of an anonymous user an object with status "NOK" is returned.
    """

    # default response
    response_data = { "status": "NOK" }

    # see if there is matching project
    project = load_project(request.project_name)

    # update response
    if project is not None:
        response_data = {
            'status': "OK",
            'name': project.name,
            'dataStart': project.datastart,
            'dataEnd': project.dataend,
            'projectionStartYear': project.datastart,
            'projectionEndYear': project.econ_dataend,
            'programs': project.programs,
            'populations': project.populations,
            'creation_time': project.creation_time,
            'data_upload_time': project.data_upload_time(),
            'has_data': project.has_data(),
            'can_calibrate': project.can_calibrate(),
            'can_scenarios': project.can_scenarios(),
        }
    return jsonify(response_data)

@project.route('/list')
@login_required
def getProjectList():
    """
    Returns the list of existing projects from db.

    Returns:
        A jsonified list of project dictionaries if the user is logged in.
        In case of an anonymous user an empty list will be returned.

    """
    projects_data = []
    # Get current user
    if current_user.is_anonymous() == False:

        # Get projects for current user
        projects = ProjectDb.query.filter_by(user_id=current_user.id)
        for project in projects:
            data_upload_time = project.creation_time
            if project.project_data: data_upload_time = project.project_data.upload_time
            project_data = {
                'status': "OK",
                'name': project.name,
                'dataStart': project.datastart,
                'dataEnd': project.dataend,
                'projectionStartYear': project.datastart,
                'projectionEndYear': project.econ_dataend,
                'programs': project.programs,
                'populations': project.populations,
                'creation_time': project.creation_time,
                'data_upload_time': data_upload_time
            }
            projects_data.append(project_data)

    return jsonify({"projects": projects_data})

@project.route('/delete/<project_name>', methods=['DELETE'])
@login_required
@report_exception()
def deleteProject(project_name):
    """
    Deletes the given project (and eventually, corresponding excel files)
    """
    current_app.logger.debug("deleteProject %s" % project_name)
    delete_spreadsheet(project_name)
    current_app.logger.debug("spreadsheets for %s deleted" % project_name)
    # Get project row for current user with project name
    project = load_project(project_name)

    if project is not None:
        id = project.id
        #delete all relevant entries explicitly
        db.session.query(ProjectDataDb).filter_by(id=id).delete()
        db.session.query(WorkingProjectDb).filter_by(id=id).delete()
        db.session.query(ProjectDb).filter_by(id=id).delete()

    db.session.commit()

    return jsonify({'status':'OK','reason':'Project %s deleted.' % project_name})

@project.route('/copy/<project_name>', methods=['POST'])
@login_required
@report_exception()
def copyProject(project_name):
    """
    Copies the given project to a different name
    usage: /api/project/copy/<project_name>?to=<new_project_name>
    """
    from sqlalchemy.orm.session import make_transient, make_transient_to_detached
    reply = BAD_REPLY
    new_project_name = request.args.get('to')
    if not new_project_name:
        reply['reason'] = 'New project name is not given'
        return reply
    # Get project row for current user with project name
    project = load_project(project_name, all_data = True)
    if project is None:
        reply['reason'] = 'Project %s does not exist.' % project_name
        return reply
    project_data_exists = project.project_data #force loading it
    db.session.expunge(project)
    make_transient(project)
    project.id = None
    project.name = new_project_name
    db.session.add(project)
    db.session.flush() #this updates the project ID to the new value
    if project_data_exists:
        db.session.expunge(project.project_data) # it should have worked without that black magic. but it didn't.
        make_transient(project.project_data)
        db.session.add(project.project_data)
    db.session.commit()
    # let's not copy working project, it should be either saved or discarded
    return jsonify({'status':'OK','project':project_name, 'copied_to':new_project_name})

@project.route('/export', methods=['POST'])
@login_required
@report_exception()
def exportGraph():
    """
    saves data as Excel file
    """
    from sim.makeworkbook import OptimaGraphTable
    data = json.loads(request.data)
    name = data['name']
    filename = name+'.xlsx'
    columns = data['columns']
    path = fullpath(filename)
    table = OptimaGraphTable(name, columns)
    table.create(path)
    (dirname, basename) = os.path.split(path)
    return helpers.send_file(path, mimetype='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')


@project.route('/update', methods=['POST'])
@login_required
@check_project_name
@report_exception()
def uploadExcel():
    """
    Uploads Excel file, uses it to update the corresponding model.
    Precondition: model should exist.
    """
    current_app.logger.debug("api/project/update")
    project_name = request.project_name
    user_id = current_user.id
    current_app.logger.debug("uploadExcel(project name: %s user:%s)" % (project_name, user_id))

    reply = {'status':'NOK'}
    file = request.files['file']

    # getting current user path
    loaddir =  upload_dir_user(DATADIR)
    if not loaddir:
        loaddir = DATADIR
    if not file:
        reply['reason'] = 'No file is submitted!'
        return json.dumps(reply)

    source_filename = secure_filename(file.filename)
    if not allowed_file(source_filename):
        reply['reason'] = 'File type of %s is not accepted!' % source_filename
        return json.dumps(reply)

    reply['file'] = source_filename

    filename = project_name + '.xlsx'
    server_filename = os.path.join(loaddir, filename)
    file.save(server_filename)

    # See if there is matching project
    project = load_project(project_name)
    if project is not None:
        # update and save model
        D = model_as_bunch(project.model)
        D = updatedata(D, savetofile = False)
        model = model_as_dict(D)
        project.model = model
        db.session.add(project)

        # save data upload timestamp
        data_upload_time = datetime.now(dateutil.tz.tzutc())
        # get file data
        filedata = open(server_filename, 'rb').read()
        # See if there is matching project data
        projdata = ProjectDataDb.query.get(project.id)

        # update existing
        if projdata is not None:
            projdata.meta = filedata
        else:
            # create new project data
            projdata = ProjectDataDb(project.id, filedata, data_upload_time)

        # Save to db
        db.session.add(projdata)
        db.session.commit()

    reply['status'] = 'OK'
    reply['result'] = 'Project %s is updated' % project_name
    return json.dumps(reply)

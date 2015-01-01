import os
from sim.dataio import DATADIR, TEMPLATEDIR, upload_dir_user
from flask import helpers, current_app
from flask.ext.login import current_user
from functools import wraps
from flask import request, jsonify, abort
from dbconn import db
from dbmodels import ProjectDb, WorkingProjectDb, UserDb
import traceback

ALLOWED_EXTENSIONS = {'txt', 'xlsx', 'xls'}

BAD_REPLY = {"status":"NOK"}

def check_project_name(api_call):
    @wraps(api_call)
    def _check_project_name(*args, **kwargs):
        reply = BAD_REPLY
        try:
            project_name = request.headers['project']
        except:
            project_name = ''

        if project_name == '':
            reply['reason'] = 'No project is open'
            return jsonify(reply)
        else:
            request.project_name = project_name
            return api_call(*args, **kwargs)
    return _check_project_name

def report_exception(reason = None):
    def _report_exception(api_call):
        @wraps(api_call)
        def __report_exception(*args, **kwargs):
            try:
                return api_call(*args, **kwargs)
            except Exception, err:
                var = traceback.format_exc()
                reply = BAD_REPLY
                reply['exception'] = var
                if reason:
                    reply['reason'] = reason
                return jsonify(reply)
        return __report_exception
    return _report_exception

#verification by secret (hashed pw)
def verify_request(api_call):
    @wraps(api_call)
    def _verify_request(*args, **kwargs):
        secret = request.args.get('secret','')
        u = UserDb.query.filter_by(password = secret, is_admin=True).first()
        if u is None:
            abort(401)
        else:
            current_app.logger.debug("admin_user: %s %s %s" % (u.name, u.password, u.email))
            return api_call(*args, **kwargs)
    return _verify_request


""" Finds out if this file is allowed to be uploaded """
def allowed_file(filename):
    return '.' in filename and \
    filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

def loaddir(app):
    loaddir = app.config['UPLOAD_FOLDER']
    if not loaddir:
        loaddir = DATADIR
    return loaddir

def project_exists(name):
    cu = current_user
    return ProjectDb.query.filter_by(user_id=cu.id, name=name).count()>0

def load_project(name):
    cu = current_user
    current_app.logger.debug("getting project %s for user %s" % (name, cu.id))
    project = ProjectDb.query.filter_by(user_id=cu.id, name=name).first()
    if project is None:
        current_app.logger.warning("no such project found: %s for user %s %s" % (name, cu.id, cu.name))
    return project

def save_data_spreadsheet(name, folder=DATADIR):
    spreadsheet_file = name
    user_dir = upload_dir_user(folder)
    if not spreadsheet_file.startswith(user_dir):
        spreadsheet_file = helpers.safe_join(user_dir, name+ '.xlsx')

def delete_spreadsheet(name):
    spreadsheet_file = name
    for parent_dir in [TEMPLATEDIR, DATADIR]:
        user_dir = upload_dir_user(TEMPLATEDIR)
        if not spreadsheet_file.startswith(user_dir):
            spreadsheet_file = helpers.safe_join(user_dir, name+ '.xlsx')
        if os.path.exists(spreadsheet_file):
            os.remove(spreadsheet_file)

def model_as_dict(model):
    from sim.bunch import Bunch
    if isinstance(model, Bunch):
        model = model.toDict()
    return model

def model_as_bunch(model):
    from sim.bunch import Bunch
    return Bunch.fromDict(model)

"""
  loads the project with the given name
  returns the model (D).
"""
def load_model(name, as_bunch = True, working_model = False):
    current_app.logger.debug("load_model:%s" % name)
    model = None
    project = load_project(name)
    if project is not None:
        if project.working_project is None or working_model == False:
            current_app.logger.debug("project %s does not have working model" % name)
            model = project.model
        else:
            current_app.logger.debug("project %s has working model" % name)
            model = project.working_project.model
        if model is None or len(model.keys())==0:
            current_app.logger.debug("model %s is None" % name)
        else:
            if as_bunch:
                model = model_as_bunch(model)
    return model

def save_model_db(name, model):
    current_app.logger.debug("save_model_db %s" % name)

    model = model_as_dict(model)
    project = load_project(name)
    project.model = model #we want it to fail if there is no project...
    db.session.add(project)
    db.session.commit()

def save_working_model(name, model):

    model = model_as_dict(model)
    project = load_project(name)

    # If we do not have an instance for working project, make it now
    if project.working_project is None:
        working_project = WorkingProjectDb(project.id, model=model, is_calibrating=True)
    else:
        project.working_project.model = model
        working_project = project.working_project

    db.session.add(working_project)
    db.session.commit()

def save_working_model_as_default(name):
    current_app.logger.debug("save_working_model_as_default %s" % name)

    project = load_project(name)
    model = project.model

    # Make sure there is a working project
    if project.working_project is not None:
        project.model = project.working_project.model
        model = project.model
        db.session.add(project)
        db.session.commit()

    return model

def revert_working_model_to_default(name):
    current_app.logger.debug("revert_working_model_to_default %s" % name)

    project = load_project(name)
    model = project.model

    # Make sure there is a working project
    if project.working_project is not None:
        project.working_project.is_calibrating = False
        project.working_project.model = model
        db.session.add(project.working_project)
        db.session.commit()

    return model

def save_model(name, model):
  try:
    save_model_db(name, model)
  except:
    save_model_file(name, model)

def pick_params(params, data, args = {}):
    for param in params:
        the_value = data.get(param)
        if the_value:
            args[param] = the_value
    return args

def for_fe(item): #only for json
    import numpy as np
    from sim.bunch import Bunch as struct

    if isinstance(item, list):
        return [for_fe(v) for v in item]
    if isinstance(item, np.ndarray):
        return [for_fe(v) for v in item.tolist()]
    elif isinstance(item, struct):
        return item.toDict()
    elif isinstance(item, dict):

        return dict( (k, for_fe(v)) for k,v in item.iteritems() )
    elif isinstance(item, float) and np.isnan(item):
        return None
    else:
        return item

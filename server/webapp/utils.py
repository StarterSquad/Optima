import os
from functools import wraps
import traceback

from dataio import TEMPLATEDIR, upload_dir_user, fromjson, tojson

from flask import helpers, current_app, abort
from flask import request, jsonify, Response, make_response
from werkzeug.datastructures import FileStorage

from flask.ext.login import current_user
from flask_restful.reqparse import RequestParser as OrigReqParser

from server.webapp.dbconn import db
from server.webapp.dbmodels import ProjectDb, UserDb, ResultsDb, ParsetsDb

import optima as op

# json should probably removed from here since we are now using prj for up/download
ALLOWED_EXTENSIONS = {'txt', 'xlsx', 'xls', 'json', 'prj', 'prg', 'par'}  # TODO this should be checked per upload type


def check_project_name(api_call):
    @wraps(api_call)
    def _check_project_name(*args, **kwargs):
        project_name = None
        project_id = None
        try:
            project_name = request.headers['project']
            project_id = request.headers['project-id']
            request.project_name = project_name
            request.project_id = project_id
            return api_call(*args, **kwargs)
        except Exception:
            exception = traceback.format_exc()
            current_app.logger.error("Exception during request %s: %s" % (request, exception))
            reply = {'reason': 'No project is open', 'exception': exception}
            return jsonify(reply), 400
    return _check_project_name


#this should be run after check_project_name
def check_project_exists(api_call):
    @wraps(api_call)
    def _check_project_exists(*args, **kwargs):
        project_id = request.headers['project-id']
        project_name = request.headers['project']
        if not project_exists(project_id):
            error_msg = 'Project %s(%s) does not exist' % (project_id, project_name)
            current_app.logger.error(error_msg)
            reply = {'reason': error_msg}
            return jsonify(reply), 404
        else:
            return api_call(*args, **kwargs)
    return _check_project_exists


def report_exception(api_call):
    @wraps(api_call)
    def _report_exception(*args, **kwargs):
        from werkzeug.exceptions import HTTPException
        try:
            return api_call(*args, **kwargs)
        except Exception, e:
            exception = traceback.format_exc()
            # limiting the exception information to 10000 characters maximum
            # (to prevent monstrous sqlalchemy outputs)
            current_app.logger.error("Exception during request %s: %.10000s" % (request, exception))
            if isinstance(e, HTTPException):
                raise
            code = 500
            reply = {'exception': exception}
            return make_response(jsonify(reply), code)
    return _report_exception


def verify_admin_request(api_call):
    """
    verification by secret (hashed pw) or by being a user with admin rights
    """
    @wraps(api_call)
    def _verify_admin_request(*args, **kwargs):
        u = None
        if (not current_user.is_anonymous()) and current_user.is_authenticated() and current_user.is_admin:
            u = current_user
        else:
            secret = request.args.get('secret', '')
            u = UserDb.query.filter_by(password=secret, is_admin=True).first()
        if u is None:
            abort(403)
        else:
            current_app.logger.debug("admin_user: %s %s %s" % (u.name, u.password, u.email))
            return api_call(*args, **kwargs)
    return _verify_admin_request


def allowed_file(filename):
    """
    Finds out if this file is allowed to be uploaded
    """
    return '.' in filename and filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS


def loaddir(app):
    the_loaddir = app.config['UPLOAD_FOLDER']
    return the_loaddir


def send_as_json_file(data):
    import json
    the_loaddir = upload_dir_user(TEMPLATEDIR)
    if not the_loaddir:
        the_loaddir = TEMPLATEDIR
    filename = 'data.json'
    server_filename = os.path.join(the_loaddir, filename)
    print "server_filename", server_filename
    with open(server_filename, 'wb') as filedata:
        json.dump(data, filedata)

    response = helpers.send_from_directory(loaddir, filename)
#    response.headers.add('content-length', str(os.path.getsize(server_filename)))
    return response


def project_exists(project_id, raise_exception=False):
    from server.webapp.exceptions import ProjectDoesNotExist
    cu = current_user
    query_args = {'id': project_id}
    if not current_user.is_admin:
        query_args['user_id'] = cu.id
    count = ProjectDb.query(**query_args)

    if raise_exception and count == 0:
        raise ProjectDoesNotExist(id=project_id)

    return count > 0


def load_project(project_id, all_data=False, raise_exception=False, db_session=None):
    from sqlalchemy.orm import undefer, defaultload
    from server.webapp.exceptions import ProjectDoesNotExist
    if not db_session:
        db_session = db.session
    cu = current_user
    current_app.logger.debug("getting project {} for user {} (admin:{})".format(
        project_id,
        cu.id if not cu.is_anonymous() else None,
        cu.is_admin if not cu.is_anonymous else False
    ))
    if cu.is_anonymous():
        if raise_exception:
            abort(401)
        else:
            return None
    if cu.is_admin:
        query = db_session.query(ProjectDb).filter_by(id=project_id)
    else:
        query = db_session.query(ProjectDb).filter_by(id=project_id, user_id=cu.id)
    if all_data:
        query = query.options(
            # undefer('model'),
            # defaultload(ProjectDb.working_project).undefer('model'),
            defaultload(ProjectDb.project_data).undefer('meta'))
    project = query.first()
    if project is None:
        current_app.logger.warning("no such project found: %s for user %s %s" % (project_id, cu.id, cu.name))
        if raise_exception:
            raise ProjectDoesNotExist(id=project_id)
    return project


def _load_project_child(project_id, record_id, record_class, exception_class, raise_exception=True):
    cu = current_user
    current_app.logger.debug("getting {} {} for user {}".format(record_class.__name__, record_id, cu.id))

    print "record_id", record_id, "record_class", record_class
    entry = db.session.query(record_class).get(record_id)
    if entry is None:
        if raise_exception:
            raise exception_class(id=record_id)
        return None

    if entry.project_id != project_id:
        if raise_exception:
            raise exception_class(id=record_id)
        return None

    if not cu.is_admin and entry.project.user_id != cu.id:
        if raise_exception:
            raise exception_class(id=record_id)
        return None

    return entry


def load_progset(project_id, progset_id, raise_exception=True):
    from server.webapp.dbmodels import ProgsetsDb
    from server.webapp.exceptions import ProgsetDoesNotExist

    return _load_project_child(project_id, progset_id, ProgsetsDb, ProgsetDoesNotExist, raise_exception)


def load_parset(project_id, parset_id, raise_exception=True):
    from server.webapp.dbmodels import ParsetsDb
    from server.webapp.exceptions import ParsetDoesNotExist

    return _load_project_child(project_id, parset_id, ParsetsDb, ParsetDoesNotExist, raise_exception)


def load_parset(project_id, parset_id, raise_exception=True):
    from server.webapp.dbmodels import ParsetsDb
    from server.webapp.exceptions import ParsetDoesNotExist

    cu = current_user
    current_app.logger.debug("getting parset {} for user {}".format(parset_id, cu.id))

    parset_entry = db.session.query(ParsetsDb).get(parset_id)
    if parset_entry is None:
        if raise_exception:
            raise ParsetDoesNotExist(id=parset_id)
        return None

    if parset_entry.project_id != project_id:
        if raise_exception:
            raise ParsetDoesNotExist(id=parset_id)
        return None

    return parset_entry


def load_program(project_id, progset_id, program_id, raise_exception=True):
    from server.webapp.dbmodels import ProgramsDb
    from server.webapp.exceptions import ProgramDoesNotExist

    cu = current_user
    current_app.logger.debug("getting project {} for user {}".format(progset_id, cu.id))

    progset_entry = load_progset(project_id, progset_id,
                                 raise_exception=raise_exception)

    program_entry = db.session.query(ProgramsDb).get(program_id)

    if program_entry.progset_id != progset_entry.id:
        if raise_exception:
            raise ProgramDoesNotExist(id=program_id)
        return None

    return program_entry


def load_scenario(project_id, scenario_id, raise_exception=True):
    from server.webapp.dbmodels import ScenariosDb
    from server.webapp.exceptions import ScenarioDoesNotExist

    cu = current_user
    current_app.logger.debug("getting scenario {} for user {}".format(scenario_id, cu.id))

    scenario_entry = db.session.query(ScenariosDb).get(scenario_id)

    if scenario_entry is None:
        if raise_exception:
            raise ScenarioDoesNotExist(id=scenario_id)
        return None

    if scenario_entry.project_id != project_id:
        if raise_exception:
            raise ScenarioDoesNotExist(id=scenario_id)
        return None

    return scenario_entry



def save_data_spreadsheet(name, folder=None):
    if folder is None:
        folder = current_app.config['UPLOAD_FOLDER']
    spreadsheet_file = name
    user_dir = upload_dir_user(folder)
    if not spreadsheet_file.startswith(user_dir):
        spreadsheet_file = helpers.safe_join(user_dir, name + '.xlsx')


def delete_spreadsheet(name, user_id=None):
    spreadsheet_file = name
    for parent_dir in [TEMPLATEDIR, current_app.config['UPLOAD_FOLDER']]:
        user_dir = upload_dir_user(parent_dir, user_id)
        if not spreadsheet_file.startswith(user_dir):
            spreadsheet_file = helpers.safe_join(user_dir, name + '.xlsx')
        if os.path.exists(spreadsheet_file):
            os.remove(spreadsheet_file)


# TODO get rid of this
def model_as_dict(model):
    return tojson(model)


# TODO get rid of this
def model_as_bunch(model):
    return fromjson(model)


def load_model(project_id, from_json=True, working_model=False):  # todo rename
    """
      loads the project with the given name
      returns the hydrated project instance (Can't think of another name than "model" yet...).
    """
    # TODO we won't have to do this for working_model, because this concept won't make sense in Optima 2.0
    current_app.logger.debug("load_model:%s" % project_id)
    model = None
    project = load_project(project_id)
    if project is not None:
        if not working_model or project.working_project is None:
            current_app.logger.debug("project %s loading main model" % project_id)
            model = project.hydrate()
        else:  # this branch won't be needed
            current_app.logger.debug("project %s loading working model" % project_id)
            model = project.working_project.hydrate()
        if model is None:
            current_app.logger.debug("model %s is None" % project_id)
#        else: todo remove from_json
#            if from_json: model = model_as_bunch(model)
    return model


def save_working_model_as_default(project_id):  # TODO will be about results, not about the model
    current_app.logger.debug("save_working_model_as_default %s" % project_id)

    project = load_project(project_id)
    model = project.model

    # Make sure there is a working project
    if project.working_project is not None:
        project.model = project.working_project.model
        model = project.model
        db.session.add(project)
        db.session.commit()

    return model


def revert_working_model_to_default(project_id):  # TODO will be about results, not about the model
    current_app.logger.debug("revert_working_model_to_default %s" % project_id)

    project = load_project(project_id, all_data=True)
    model = project.model

    # Make sure there is a working project
    if project.working_project is not None:
        project.working_project.is_calibrating = False
        project.working_project.model = model
        db.session.add(project.working_project)
        db.session.commit()

    return model


def save_model(project_id, model, to_json=False):
    # model is given as json by default, no need to convert
    current_app.logger.debug("save_model %s" % project_id)

    if to_json:
        model = model_as_dict(model)
    project = load_project(project_id)
    project.model = model  # we want it to fail if there is no project...
    db.session.add(project)
    db.session.commit()


def pick_params(params, data, args=None):
    if args is None:
        args = {}
    for param in params:
        the_value = data.get(param)
        if the_value:
            args[param] = the_value
    return args


def for_fe(item):  # only for json
    import numpy as np

    if isinstance(item, list):
        return [for_fe(v) for v in item]
    if isinstance(item, np.ndarray):
        return [for_fe(v) for v in item.tolist()]
    elif isinstance(item, dict):
        return dict((k, for_fe(v)) for k, v in item.iteritems())
    elif isinstance(item, float) and np.isnan(item):
        return None
    else:
        return item


def update_or_create_parset(project_id, name, parset):

    from datetime import datetime
    import dateutil
    from server.webapp.dbmodels import ParsetsDb
    from optima.utils import saves

    parset_record = ParsetsDb.query.filter_by(id=parset.uid, project_id=project_id).first()
    if parset_record is None:
        parset_record = ParsetsDb(
            id=parset.uid,
            project_id=project_id,
            name=name,
            created=parset.created or datetime.now(dateutil.tz.tzutc()),
            updated=parset.modified or datetime.now(dateutil.tz.tzutc()),
            pars=saves(parset.pars)
        )

        db.session.add(parset_record)
    else:
        parset_record.updated = datetime.now(dateutil.tz.tzutc())
        parset_record.name = name
        parset_record.pars = saves(parset.pars)
        db.session.add(parset_record)


def update_or_create_progset(project_id, name, progset):

    from datetime import datetime
    import dateutil
    from server.webapp.dbmodels import ProgsetsDb

    progset_record = ProgsetsDb.query \
        .filter_by(project_id=project_id, name=name) \
        .first()

    if progset_record is None:
        progset_record = ProgsetsDb(
            project_id=project_id,
            name=name,
            created=progset.created or datetime.now(dateutil.tz.tzutc()),
            updated=datetime.now(dateutil.tz.tzutc())
        )

        db.session.add(progset_record)
        db.session.flush()
    else:
        progset_record.updated = datetime.now(dateutil.tz.tzutc())
        db.session.add(progset_record)

    return progset_record


def update_or_create_program(project_id, progset_id, name, program, active=False):

    from datetime import datetime
    import dateutil
    from server.webapp.dbmodels import ProgramsDb
    from optima.utils import saves

    program_record = ProgramsDb.query \
        .filter_by(
            short=program.get('short', None),
            project_id=project_id,
            progset_id=progset_id
        ).first()

    if program_record is None:
        program_record = ProgramsDb(
            project_id=project_id,
            progset_id=progset_id,
            name=name,
            short=program.get('short', ''),
            category=program.get('category', ''),
            created=datetime.now(dateutil.tz.tzutc()),
            updated=datetime.now(dateutil.tz.tzutc()),
            pars=ProgramsDb.program_pars_to_pars(program.get('targetpars', [])),
            targetpops=program.get('targetpops', []),
            active=active,
            criteria=program.get('criteria', None),
            costcov=program.get('costcov', [])
        )

    else:
        program_record.updated = datetime.now(dateutil.tz.tzutc())
        program_record.pars = ProgramsDb.program_pars_to_pars(program.get('targetpars', []))
        program_record.targetpops = program.get('targetpops', [])
        program_record.short = program.get('short', '')
        program_record.category = program.get('category', '')
        program_record.active = active
        program_record.criteria = program.get('criteria', None)
        program_record.costcov = program.get('costcov', [])

    program_record.blob = saves(program_record.hydrate())
    db.session.add(program_record)
    return program_record


def modify_program(project_id, progset_id, program_id, args, program_modifier):
    # looks up a program, hydrates it, calls a modifier
    # (function defined somewhere) with given args, saves the result
    # TODO could such things be done as decorators?
    program_entry = load_program(project_id, progset_id, program_id)
    if program_entry is None:
        raise ProgramDoesNotExist(id=program_id, project_id=project_id)
    program_instance = program_entry.hydrate()
    program_modifier(program_instance, args)
    program_entry.restore(program_instance)
    result = {"params": program_entry.ccopars or {},
              "data": program_entry.data_db_to_api()}
    db.session.add(program_entry)
    db.session.commit()
    return result


def update_or_create_scenario(project_id, project, name):  # project might have a different ID than the one we want

    from datetime import datetime
    import dateutil
    from server.webapp.dbmodels import ScenariosDb, ParsetsDb, ProgsetsDb
    import json

    parset_id = None
    progset_id = None
    blob = {}
    scenario_type = None

    scenario = project.scens[name]
    if not scenario:
        raise Exception("scenario {} not present in project {}!".format(name, project_id))

    if isinstance(scenario, op.Parscen):
        scenario_type = 'parameter'
    elif isinstance(scenario, op.Budgetscen):
        scenario_type = 'budget'
    elif isinstance(scenario, op.Coveragescen):
        scenario_type = 'coverage'

    if scenario.t:
        blob['years'] = scenario.t
    for key in ['budget', 'coverage', 'pars']:
        if hasattr(scenario, key) and getattr(scenario, key):
            blob[key] = json.loads(json.dumps(getattr(scenario, key)))

    parset_name = scenario.parsetname
    if parset_name:
        parset_record = ParsetsDb.query \
            .filter_by(project_id=project_id, name=parset_name) \
            .first()
        if parset_record:
            parset_id = parset_record.id

    progset_id = None
    if hasattr(scenario, 'progsetname') and scenario.progsetname:
        progset_name = scenario.progsetname
        progset_record = ProgsetsDb.query \
            .filter_by(project_id=project_id, name=progset_name) \
            .first()
        if progset_record:
            progset_id = progset_record.id

    scenario_record = ScenariosDb.query \
        .filter_by(project_id=project_id, name=name) \
        .first()

    if scenario_record is None:
        scenario_record = ScenariosDb(
            project_id=project_id,
            parset_id=parset_id,
            progset_id=progset_id,
            name=name,
            scenario_type=scenario_type,
            active=scenario.active,
            blob=blob
        )
        db.session.add(scenario_record)
        db.session.flush()
    else:
        scenario_record.parset_id = parset_id
        scenario_record.scenario_type = scenario_type
        scenario_record.active=scenario.active
        scenario_record.blob = blob
        db.session.add(scenario_record)

    return scenario_record


def update_or_create_optimization(project_id, project, name):

    from datetime import datetime
    import dateutil
    from server.webapp.dbmodels import OptimizationsDb, ProgsetsDb, ParsetsDb
    from optima.utils import saves

    parset_id = None
    progset_id = None

    optim = project.optims[name]
    if not optim:
        raise Exception("optimization {} not present in project {}!".format(name, project_id))

    parset_name = optim.parsetname
    if parset_name:
        parset_record = ParsetsDb.query \
        .filter_by(project_id=project_id, name=parset_name) \
        .first()
        if parset_record:
            parset_id = parset_record.id

    progset_name = optim.progsetname
    if progset_name:
        progset_record = ProgsetsDb.query \
        .filter_by(project_id=project_id, name=progset_name) \
        .first()
        if progset_record:
            progset_id = progset_record.id

    optimization_record = OptimizationsDb.query.filter_by(name=name, project_id=project_id).first()
    if optimization_record is None:
        optimization_record = OptimizationsDb(
            project_id=project_id,
            name=name,
            which = optim.objectives.get('which', 'outcome') if optim.objectives else 'outcome',
            parset_id=parset_id,
            progset_id=progset_id,
            objectives=(optim.objectives or {}),
            constraints=(optim.constraints or {})
        )
        db.session.add(optimization_record)
        db.session.flush()
    else:
        optimization_record.which = optim.objectives.get('which', 'outcome') if optim.objectives else 'outcome'
        optimization_record.parset_id = parset_id
        optimization_record.progset_id = progset_id
        optimization_record.objectives = (optim.objectives or {})
        optimization_record.constraints = (optim.constraints or {})
        db.session.add(optimization_record)

    return optimization_record

def save_result(project_id, result, parset_name='default', calculation_type = ResultsDb.CALIBRATION_TYPE,
    db_session=None):
    if not db_session:
        db_session=db.session
    # find relevant parset for the result
    print("save_result(%s, %s, %s" % (project_id, parset_name, calculation_type))
    project_parsets = db_session.query(ParsetsDb).filter_by(project_id=project_id)
    default_parset = [item for item in project_parsets if item.name == parset_name]
    if default_parset:
        default_parset = default_parset[0]
    else:
        raise Exception("parset '{}' not generated for the project {}!".format(parset_name, project_id))
    result_parset_id = default_parset.id

    # update results (after runsim is invoked)
    project_results = db_session.query(ResultsDb).filter_by(project_id=project_id)

    result_record = [item for item in project_results if
                     item.parset_id == result_parset_id and
                     item.calculation_type == calculation_type]
    if result_record:
        if len(result_record) > 1:
            abort(500, "Found multiple records for result")
        result_record = result_record[0]
        result_record.blob = op.saves(result)
    if not result_record:
        result_record = ResultsDb(
            parset_id=result_parset_id,
            project_id=project_id,
            calculation_type=calculation_type,
            blob=op.saves(result)
        )
    return result_record


def remove_nans(obj):
    import json
    # a hack to get rid of NaNs, javascript JSON parser doesn't like them
    json_string = json.dumps(obj).replace('NaN', 'null')
    return json.loads(json_string)



def init_login_manager(login_manager):

    @login_manager.user_loader
    def load_user(userid):
        from server.webapp.dbmodels import UserDb
        try:
            user = UserDb.query.filter_by(id=userid).first()
        except Exception:
            user = None
        return user

    @login_manager.request_loader
    def load_user_from_request(request):  # pylint: disable=redefined-outer-name

        # try to login using the secret url arg
        secret = request.args.get('secret')
        if secret:
            from server.webapp.dbmodels import UserDb
            user = UserDb.query.filter_by(password=secret, is_admin=True).first()
            if user:
                return user

        # finally, return None if both methods did not login the user
        return None

    @login_manager.unauthorized_handler
    def unauthorized_handler():
        abort(401)


class RequestParser(OrigReqParser):

    def __init__(self, *args, **kwargs):
        super(RequestParser, self).__init__(*args, **kwargs)
        self.abort_on_error = True

    def get_swagger_type(self, arg):
        try:
            if issubclass(arg.type, FileStorage):
                return 'file'
        except TypeError:
            ## this arg.type was not a class
            pass

        if callable(arg.type):
            return arg.type.__name__
        return arg.type

    def get_swagger_location(self, arg):

        if isinstance(arg.location, tuple):
            loc = arg.location[0]
        else:
            loc = arg.location.split(',')[0]

        if loc == "args":
            return "query"
        return loc


    def swagger_parameters(self):
        return [
            {
                'name': arg.name,
                'dataType': self.get_swagger_type(arg),
                'required': arg.required,
                'description': arg.help,
                'paramType': self.get_swagger_location(arg),
            }
            for arg in self.args
        ]

    def add_arguments(self, arguments_dict):
        for argument_name, kwargs in arguments_dict.iteritems():
            self.add_argument(argument_name, **kwargs)

    def parse_args(self, req=None, strict=False):
        from werkzeug.exceptions import HTTPException

        try:
            return super(RequestParser, self).parse_args(req, strict)
        except HTTPException as e:
            if self.abort_on_error:
                raise e
            else:
                raise ValueError(e.data['message'])

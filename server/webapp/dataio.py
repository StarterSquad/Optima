__doc__ = """

dataio.py
=========

Contains all the functions that fetches and saves optima objects to/from database
and the file system. These functions abstracts out the data i/o for the web-server
api calls.

Function call pairs are load_*/save_* and refers to saving to database.

Database record variables should have suffix _record

Parsed data structures should have suffix _summary

All parameters and return types are either id's, json-summaries, or mpld3 graphs
"""

import traceback
from functools import wraps

from pprint import pprint
from server.webapp.dbconn import db
from server.webapp.dbmodels import UserDb
from server.webapp.exceptions import UserAlreadyExists, UserDoesNotExist, InvalidCredentials
from server.webapp.utils import nullable_email, hashed_password

import os
from zipfile import ZipFile
from uuid import uuid4, UUID
from datetime import datetime
import dateutil

from flask import helpers, current_app, abort, request, session, make_response, jsonify
from flask.ext.login import current_user, login_user, logout_user

from werkzeug.utils import secure_filename
from flask.ext.restful import marshal

from optima.dataio import loadobj as loaddbobj
import optima as op
import optima

from .dbmodels import ProjectDb, ResultsDb, ProjectDataDb, ProjectEconDb, PyObjectDb
from .exceptions import ProjectDoesNotExist, ParsetAlreadyExists
from .parse import get_default_program_summaries, \
    get_parameters_for_edit_program, get_parameters_for_outcomes, \
    get_parameters_from_parset, set_parameters_on_parset, \
    get_progset_from_project, get_populations_from_project, \
    set_populations_on_project, set_project_summary_on_project, \
    get_project_summary_from_project, get_parset_from_project, \
    get_parset_summaries, set_scenario_summaries_on_project, \
    get_scenario_summaries, get_parameters_for_scenarios, \
    get_optimization_summaries, get_default_optimization_summaries, \
    set_optimization_summaries_on_project, get_optimization_from_project, \
    get_program_from_progset, get_project_years, get_progset_summaries, \
    set_progset_summary_on_project, get_progset_summary, \
    get_outcome_summaries_from_progset, set_outcome_summaries_on_progset, \
    set_program_summary_on_progset, parse_portfolio_summary, \
    set_progset_summary_on_progset
from .plot import make_mpld3_graph_dict, convert_to_mpld3
from .utils import TEMPLATEDIR, templatepath, upload_dir_user, normalize_obj
from . import parse

# USERS

def authenticate_current_user(raise_exception=True):
    current_app.logger.debug("authenticating user {} (admin:{})".format(
        current_user.id if not current_user.is_anonymous() else None,
        current_user.is_admin if not current_user.is_anonymous else False
    ))
    if current_user.is_anonymous():
        if raise_exception:
            abort(401)
        else:
            return None


def marshal_user(query):
    return marshal(query, UserDb.resource_fields)


def get_users():
    return marshal_user(UserDb.query.all())


def parse_user_args(args):
    return {
        'email': nullable_email(args.get('email', None)),
        'name': args.get('displayName', ''),
        'username': args.get('username', ''),
        'password': hashed_password(args.get('password')),
    }


def create_user(args):
    n_user = UserDb.query.filter_by(username=args['username']).count()
    if n_user > 0:
        raise UserAlreadyExists(args['username'])

    user = UserDb(**args)
    db.session.add(user)
    db.session.commit()

    return marshal_user(user)


def update_user(user_id, args):
    user = UserDb.query.get(user_id)
    if user is None:
        raise UserDoesNotExist(user_id)

    try:
        userisanonymous = current_user.is_anonymous()  # CK: WARNING, SUPER HACKY way of dealing with different Flask versions
    except:
        userisanonymous = current_user.is_anonymous

    if userisanonymous or (str(user_id) != str(current_user.id) and not current_user.is_admin):
        secret = request.args.get('secret', '')
        u = UserDb.query.filter_by(password=secret, is_admin=True).first()
        if u is None:
            abort(403)

    for key, value in args.iteritems():
        if value is not None:
            setattr(user, key, value)

    db.session.commit()

    return marshal_user(user)


def do_login_user(args):
    try:
        userisanonymous = current_user.is_anonymous()  # CK: WARNING, SUPER HACKY way of dealing with different Flask versions
    except:
        userisanonymous = current_user.is_anonymous


    if userisanonymous:
        current_app.logger.debug("current user anonymous, proceed with logging in")

        print("login args", args)
        try:
            # Get user for this username
            user = UserDb.query.filter_by(username=args['username']).first()

            print("user", user)
            # Make sure user is valid and password matches
            if user is not None and user.password == args['password']:
                login_user(user)
                return marshal_user(user)

        except Exception:
            var = traceback.format_exc()
            print("Exception when logging user {}: \n{}".format(args['username'], var))

        raise InvalidCredentials

    else:
        return marshal_user(current_user)


def delete_user(user_id):
    user = UserDb.query.get(user_id)

    if user is None:
        raise UserDoesNotExist(user_id)

    user_email = user.email
    user_name = user.username
    from server.webapp.dbmodels import ProjectDb
    from sqlalchemy.orm import load_only

    # delete all corresponding projects and working projects as well
    # project and related records delete should be on a method on the project model
    projects = ProjectDb.query.filter_by(user_id=user_id).options(load_only("id")).all()
    for project in projects:
        project.recursive_delete()

    db.session.delete(user)
    db.session.commit()

    print("deleted user:{} {} {}".format(user_id, user_name, user_email))


def do_logout_current_user():
    logout_user()
    session.clear()


def report_exception_decorator(api_call):
    @wraps(api_call)
    def _report_exception(*args, **kwargs):
        from werkzeug.exceptions import HTTPException
        try:
            return api_call(*args, **kwargs)
        except Exception as e:
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


def verify_admin_request_decorator(api_call):
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

## PROJECT

def load_project_record(project_id, raise_exception=True, db_session=None, authenticate=False):
    if not db_session:
        db_session = db.session

    if authenticate:
        authenticate_current_user()

    if authenticate is False or current_user.is_admin:
        query = db_session.query(ProjectDb).filter_by(id=project_id)
    else:
        query = db_session.query(ProjectDb).filter_by(
            id=project_id, user_id=current_user.id)

    project_record = query.first()

    if project_record is None:
        if raise_exception:
            raise ProjectDoesNotExist(id=project_id)

    return project_record


def load_project(project_id, raise_exception=True, db_session=None, authenticate=True):
    if not db_session:
        db_session = db.session
    project_record = load_project_record(
        project_id,
        raise_exception=raise_exception,
        db_session=db_session,
        authenticate=authenticate)
    if project_record is None:
        if raise_exception:
            raise ProjectDoesNotExist(id=project_id)
        else:
            return None
    return project_record.load()


def load_project_name(project_id):
    return load_project(project_id).name


def update_project(project, db_session=None):
    if db_session is None:
        db_session = db.session
    project_record = load_project_record(project.uid)
    project_record.save_obj(project)
    project.modified = datetime.now(dateutil.tz.tzutc())
    project_record.updated = project.modified
    db_session.add(project_record)
    db_session.commit()


def update_project_with_fn(project_id, update_project_fn, db_session=None):
    if db_session is None:
        db_session = db.session
    project_record = load_project_record(project_id, db_session=db_session)
    project = project_record.load()
    update_project_fn(project)
    project.modified = datetime.now(dateutil.tz.tzutc())
    project_record.updated = project.modified
    project_record.save_obj(project)
    db_session.add(project_record)
    db_session.commit()


def load_project_summary_from_project_record(project_record):
    try:
        project = project_record.load()
    except:
        return {
            'id': project_record.id,
            'name': "Failed loading"
        }
    project_summary = get_project_summary_from_project(project)
    project_summary['userId'] = project_record.user_id
    return project_summary


def load_project_summary(project_id):
    project_entry = load_project_record(project_id)
    return load_project_summary_from_project_record(project_entry)


def load_project_summaries(user_id=None):
    if user_id is None:
        query = ProjectDb.query
    else:
        query = ProjectDb.query.filter_by(user_id=current_user.id)
    return map(load_project_summary_from_project_record, query.all())


def create_project_with_spreadsheet_download(user_id, project_summary):
    project_entry = ProjectDb(user_id=user_id)
    db.session.add(project_entry)
    db.session.flush()

    project = op.Project(name=project_summary["name"])
    project.created = datetime.now(dateutil.tz.tzutc())
    project.modified = datetime.now(dateutil.tz.tzutc())
    project.uid = project_entry.id
    set_populations_on_project(project, project_summary["populations"])
    project.data["years"] = (
        project_summary['dataStart'], project_summary['dataEnd'])
    project_entry.save_obj(project)
    db.session.commit()

    new_project_template = secure_filename(
        "{}.xlsx".format(project_summary['name']))
    path = templatepath(new_project_template)
    op.makespreadsheet(
        path,
        pops=project_summary['populations'],
        datastart=project_summary['dataStart'],
        dataend=project_summary['dataEnd'])

    print("> Project data template: %s" % new_project_template)

    return project.uid, upload_dir_user(TEMPLATEDIR), new_project_template


def delete_projects(project_ids):
    for project_id in project_ids:
        record = load_project_record(project_id, raise_exception=True)
        record.recursive_delete()
    db.session.commit()


def update_project_followed_by_template_data_spreadsheet(project_id, project_summary):
    project_entry = load_project_record(project_id)
    project = project_entry.load()
    set_project_summary_on_project(project, project_summary)
    project_entry.save_obj(project)
    db.session.add(project_entry)
    db.session.commit()

    secure_project_name = secure_filename(project.name)
    new_project_template = secure_project_name
    path = templatepath(secure_project_name)
    op.makespreadsheet(
        path,
        pops=project_summary['populations'],
        datastart=project_summary["dataStart"],
        dataend=project_summary["dataEnd"])

    (dirname, basename) = (
        upload_dir_user(TEMPLATEDIR), new_project_template)

    return dirname, basename


def save_project_as_new(project, user_id):
    project_record = ProjectDb(user_id)
    db.session.add(project_record)
    db.session.flush()

    project.uid = project_record.id

    for attr in ['parsets','progsets','scens','optims']:
        for obj in getattr(project,attr).values():
            obj.uid = op.uuid()

    project.created = datetime.now(dateutil.tz.tzutc())
    project.modified = datetime.now(dateutil.tz.tzutc())
    project_record.created = project.created
    project_record.updated = project.modified

    project_record.save_obj(project)
    db.session.flush()

    db.session.commit()


def copy_project(project_id, new_project_name):
    # Get project row for current user with project name
    project_record = load_project_record(
        project_id, raise_exception=True)
    user_id = project_record.user_id

    project = project_record.load()
    project.name = new_project_name
    save_project_as_new(project, user_id)

    parset_name_by_id = {parset.uid: name for name, parset in project.parsets.items()}
    copy_project_id = project.uid

    # copy each result
    result_records = project_record.results
    if result_records:
        for result_record in result_records:
            # reset the parset_id in results to new project
            result = result_record.load()
            parset_id = result_record.parset_id
            if parset_id not in parset_name_by_id:
                continue
            parset_name = parset_name_by_id[parset_id]
            new_parset = [r for r in project.parsets.values() if r.name == parset_name]
            if not new_parset:
                continue
            copy_parset_id = new_parset[0].uid

            copy_result_record = ResultsDb(
                copy_parset_id, copy_project_id, result_record.calculation_type)
            db.session.add(copy_result_record)
            db.session.flush()

            # serializes result with new
            result.uid = copy_result_record.id
            copy_result_record.save_obj(result)

    db.session.commit()

    return copy_project_id


def create_project_from_prj(prj_filename, project_name, user_id):
    project = loaddbobj(prj_filename)
    print('>> Migrating project from version %s' % project.version)
    project = op.migrate(project)
    print('>> ...to version %s' % project.version)
    project.name = project_name
    resolve_project(project)
    save_project_as_new(project, user_id)
    return project.uid


def create_project_from_spreadsheet(prj_filename, project_name, user_id):
    project = op.Project(spreadsheet=prj_filename)
    project = op.migrate(project)
    project.name = project_name
    resolve_project(project)
    save_project_as_new(project, user_id)
    return project.uid


def download_project(project_id):
    project_record = load_project_record(project_id, raise_exception=True)
    dirname = upload_dir_user(TEMPLATEDIR)
    if not dirname:
        dirname = TEMPLATEDIR
    filename = project_record.as_file(dirname)
    return dirname, filename


def update_project_from_prj(project_id, prj_filename):
    project = loaddbobj(prj_filename)
    print('>> Migrating project from version %s' % project.version)
    project = op.migrate(project)
    print('>> ...to version %s' % project.version)
    project_record = load_project_record(project_id)
    project_record.save_obj(project)
    db.session.add(project_record)
    db.session.commit()


def load_zip_of_prj_files(project_ids):
    dirname = upload_dir_user(TEMPLATEDIR)
    if not dirname:
        dirname = TEMPLATEDIR

    prjs = [load_project_record(id).as_file(dirname) for id in project_ids]

    zip_fname = '{}.zip'.format(uuid4())
    server_zip_fname = os.path.join(dirname, zip_fname)
    with ZipFile(server_zip_fname, 'w') as zipfile:
        for prj in prjs:
            zipfile.write(os.path.join(dirname, prj), 'portfolio/{}'.format(prj))

    return dirname, zip_fname


## PORTFOLIO



def create_portfolio(name, db_session=None):
    if db_session is None:
        db_session = db.session
    print("> Create portfolio %s" % name)
    portfolio = op.Portfolio()
    portfolio.name = name
    objectives = optima.portfolio.defaultobjectives(verbose=0)
    gaoptim = optima.portfolio.GAOptim(objectives=objectives)
    portfolio.gaoptims[str(gaoptim.uid)] = gaoptim
    record = PyObjectDb(
        user_id=current_user.id, name=name, id=portfolio.uid, type="portfolio")
    # something about dates
    record.save_obj(portfolio)
    db_session.add(record)
    db_session.commit()
    return parse_portfolio_summary(portfolio)


def delete_portfolio(portfolio_id, db_session=None):
    if db_session is None:
        db_session = db.session
    print("> Delete portfolio %s" % portfolio_id)
    record = db_session.query(PyObjectDb).get(portfolio_id)
    record.cleanup()
    db_session.delete(record)
    db_session.commit()
    return load_portfolio_summaries()


def load_portfolio(portfolio_id, db_session=None):
    if db_session is None:
        db_session = db.session
    kwargs = {'id': portfolio_id, 'type': "portfolio"}
    record = db_session.query(PyObjectDb).filter_by(**kwargs).first()
    if record:
        print("> load portfolio %s" % portfolio_id)
        return record.load()
    return optima.loadobj("server/example/malawi-decent-two-state.prt", verbose=0)


def load_portfolio_summaries(db_session=None):
    import optima.portfolio
    if db_session is None:
        db_session = db.session

    query = db_session.query(PyObjectDb).filter_by(user_id=current_user.id)
    if query is None:
        portfolio = optima.loadobj("server/example/malawi-decent-two-state.prt", verbose=0)
        record = PyObjectDb(
            user_id=current_user.id, type="portfolio", name=portfolio.name, id=portfolio.uid)
        record.save_obj(portfolio)
        db_session.add(record)
        db_session.commit()
        print("> Crreated default portfolio %s" % portfolio.name)
        portfolios = [portfolio]
    else:
        portfolios = []
        for record in query:
            print(">> Portfolio id %s" % record.id)
            portfolio = record.load()
            if len(portfolio.gaoptims) == 0:
                objectives = optima.portfolio.defaultobjectives(verbose=0)
                gaoptim = optima.portfolio.GAOptim(objectives=objectives)
                portfolio.gaoptims[str(gaoptim.uid)] = gaoptim
                record.save_obj(portfolio)
            portfolios.append(portfolio)

    summaries = map(parse_portfolio_summary, portfolios)
    print("> Loading portfolio summaries")
    pprint(summaries, indent=2)

    return summaries


def save_portfolio(portfolio, db_session=None):
    if db_session is None:
        db_session = db.session
    portfolio_id = portfolio.uid
    kwargs = {'id': portfolio_id, 'type': "portfolio"}
    query = db_session.query(PyObjectDb).filter_by(**kwargs)
    if query:
        record = query.first()
    else:
        record = PyObjectDb(user_id=current_user.id)
        record.id = UUID(portfolio_id)
        record.type = "portfolio"
        record.name = portfolio.name
    print ">> Saved portfolio %s" % (portfolio_id)
    record.save_obj(portfolio)
    db_session.add(record)
    db_session.commit()


def set_portfolio_summary_on_portfolio(portfolio, summary):
    gaoptim_summaries = summary['gaoptims']
    gaoptims = portfolio.gaoptims
    for gaoptim_summary in gaoptim_summaries:
        gaoptim_id = str(gaoptim_summary['id'])
        objectives = optima.odict(gaoptim_summary["objectives"])
        if gaoptim_id in gaoptims:
            gaoptim = gaoptims[gaoptim_id]
            gaoptim.objectives = objectives
        else:
            gaoptim = optima.portfolio.GAOptim(objectives=objectives)
            gaoptims[gaoptim_id] = gaoptim
    old_project_ids = portfolio.projects.keys()
    print("> old project ids %s" % old_project_ids)
    new_project_ids = [s["id"] for s in summary["projects"]]
    print("> new project ids %s" % new_project_ids)
    for old_project_id in old_project_ids:
        if old_project_id not in new_project_ids:
            portfolio.projects.pop(old_project_id)
    for new_project_id in new_project_ids:
        if new_project_id not in portfolio.projects:
            project = load_project(new_project_id)
            portfolio.projects[new_project_id] = project


def load_or_create_portfolio(portfolio_id, db_session=None):
    if db_session is None:
        db_session = db.session
    kwargs = {'id': portfolio_id, 'type': "portfolio"}
    record = db_session.query(PyObjectDb).filter_by(**kwargs).first()
    if record:
        print("> load portfolio %s" % portfolio_id)
        portfolio = record.load()
    else:
        print("> Create portfolio %s" % portfolio_id)
        portfolio = optima.Portfolio()
        portfolio.uid = UUID(portfolio_id)
    return portfolio


def save_portfolio_by_summary(portfolio_id, portfolio_summary, db_session=None):
    portfolio = load_or_create_portfolio(portfolio_id)
    set_portfolio_summary_on_portfolio(portfolio, portfolio_summary)
    save_portfolio(portfolio, db_session)
    return load_portfolio_summaries()


def delete_portfolio_project(portfolio_id, project_id):
    portfolio = load_portfolio(portfolio_id)
    portfolio.projects.pop(str(project_id))
    print ">> Deleted project %s from portfolio %s" % (project_id, portfolio_id)
    save_portfolio(portfolio)
    return load_portfolio_summaries()

## PARSET


def copy_parset(project_id, parset_id, new_parset_name):

    def update_project_fn(project):
        original_parset = get_parset_from_project(project, parset_id)
        original_parset_name = original_parset.name
        project.copyparset(orig=original_parset_name, new=new_parset_name)
        project.parsets[new_parset_name].uid = op.uuid()

    update_project_with_fn(project_id, update_project_fn)


def delete_parset(project_id, parset_id):

    def update_project_fn(project):
        parset = get_parset_from_project(project, parset_id)
        project.parsets.pop(parset.name)
        resolve_project(project)

    update_project_with_fn(project_id, update_project_fn)
    delete_result_by_parset_id(project_id, parset_id)


def rename_parset(project_id, parset_id, new_parset_name):

    def update_project_fn(project):
        parset = get_parset_from_project(project, parset_id)
        old_parset_name = parset.name
        parset.name = new_parset_name
        print(">> old parsets '%s'" % project.parsets.keys())
        del project.parsets[old_parset_name]
        project.parsets[new_parset_name] = parset
        print(">> new parsets '%s'" % project.parsets.keys())

    update_project_with_fn(project_id, update_project_fn)


def create_parset(project_id, new_parset_name):

    def update_project_fn(project):
        if new_parset_name in project.parsets:
            raise ParsetAlreadyExists(project_id, new_parset_name)
        project.makeparset(new_parset_name, overwrite=False)

    update_project_with_fn(project_id, update_project_fn)


def load_parset_summaries(project_id):
    project = load_project(project_id)
    return get_parset_summaries(project)


def load_project_parameters(project_id):
    return get_parameters_for_edit_program(load_project(project_id))


def load_parameters_from_progset_parset(project_id, progset_id, parset_id):
    project = load_project(project_id)
    return get_parameters_for_outcomes(project, progset_id, parset_id)


def load_parameters(project_id, parset_id):
    project = load_project(project_id)
    parset = get_parset_from_project(project, parset_id)
    return get_parameters_from_parset(parset)


def save_parameters(project_id, parset_id, parameters):

    def update_project_fn(project):
        parset = get_parset_from_project(project, parset_id)
        print ">> Updating parset '%s'" % parset.name
        set_parameters_on_parset(parameters, parset)

    update_project_with_fn(project_id, update_project_fn)

    delete_result_by_parset_id(project_id, parset_id)


def load_parset_graphs(
        project_id, parset_id, calculation_type, which=None,
        parameters=None, start=None, end=None): # Adding hooks for selecting year range

    project = load_project(project_id)
    parset = get_parset_from_project(project, parset_id)

    if parameters is not None:
        print ">> Updating parset '%s'" % parset.name
        set_parameters_on_parset(parameters, parset)
        delete_result_by_parset_id(project_id, parset_id)
        update_project(project)

    result = load_result(project_id, parset_id, calculation_type)
    if result is None:
        print ">> Runsim for for parset '%s'" % parset.name
        result = project.runsim(name=parset.name, start=start, end=end) # ACTUALLY RUN THE MODEL
        result_record = update_or_create_result_record_by_id(result, project_id, parset_id, calculation_type)
        db.session.add(result_record)
        db.session.commit()

    assert result is not None

    print ">> Generating graphs for parset '%s'" % parset.name
    graph_dict = make_mpld3_graph_dict(result, which)

    return {
        "parameters": get_parameters_from_parset(parset),
        "graphs": graph_dict["graphs"]
    }


# RESULT


def load_result(project_id, parset_id, calculation_type=ResultsDb.DEFAULT_CALCULATION_TYPE):
    result_record = db.session.query(ResultsDb).filter_by(
        project_id=project_id, parset_id=parset_id, calculation_type=calculation_type).first()
    if result_record is None:
        return None
    return result_record.load()


def load_result_by_id(result_id):
    result_record = db.session.query(ResultsDb).get(result_id)
    if result_record is None:
        raise Exception("Results '%s' does not exist" % result_id)
    return result_record.load()


def update_or_create_result_record_by_id(
        result,
        project_id,
        parset_id,
        calculation_type=ResultsDb.DEFAULT_CALCULATION_TYPE,
        db_session=None):

    if db_session is None:
        db_session = db.session

    result_record = db_session.query(ResultsDb).get(result.uid)
    if result_record is not None:
        print ">> Updating record for result '%s'" % (result.name)
    else:
        result_record = ResultsDb(
            parset_id=parset_id,
            project_id=project_id,
            calculation_type=calculation_type)
        print ">> Creating record for result '%s'" % (result.name)

    result_record.id = result.uid
    result_record.save_obj(result)
    db_session.add(result_record)

    return result_record


def delete_result_by_parset_id(
        project_id, parset_id, calculation_type=None, db_session=None):
    if db_session is None:
        db_session = db.session
    if calculation_type is None:
        records = db_session.query(ResultsDb).filter_by(
            project_id=project_id, parset_id=parset_id)
    else:
        records = db_session.query(ResultsDb).filter_by(
            project_id=project_id, parset_id=parset_id,
            calculation_type=calculation_type)
    for record in records:
        record.cleanup()
    records.delete()
    db_session.commit()


def delete_result_by_name(
        project_id, result_name, db_session=None):
    if db_session is None:
        db_session = db.session

    records = db_session.query(ResultsDb).filter_by(project_id=project_id)
    for record in records:
        result = record.load()
        if result.name == result_name:
            print ">> Deleting outdated result '%s'" % result_name
            record.cleanup()
            db_session.delete(record)
    db_session.commit()


def load_result_csv(result_id):
    dirname = upload_dir_user(TEMPLATEDIR)
    if not dirname:
        dirname = TEMPLATEDIR
    filestem = 'results'
    filename = filestem + '.csv'

    result = load_result_by_id(result_id)
    result.export(filestem=os.path.join(dirname, filestem))

    return dirname, filename


def load_result_by_optimization(project, optimization):

    result_name = "optim-" + optimization.name
    parset_id = project.parsets[optimization.parsetname].uid

    print(">> Loading result '%s'" % result_name)
    result_records = db.session.query(ResultsDb).filter_by(
        project_id=project.uid,
        parset_id=parset_id,
        calculation_type="optimization")

    for result_record in result_records:
        result = result_record.load()
        if result.name == result_name:
            return result

    print(">> Not found result '%s'" % (optimization.name))

    return None


def load_result_mpld3_graphs(result_id, which):
    result = load_result_by_id(result_id)
    return make_mpld3_graph_dict(result, which)


## SCENARIOS


def make_scenarios_graphs(project_id):
    result = load_result(project_id, None, "scenarios")
    if result is None:
        project = load_project(project_id)
        if len(project.scens) == 0:
            print(">> No scenarios in project")
            return {}
        print(">> Run scenarios for project '%s'" % project_id)
        project.runscenarios()
        result = project.results[-1]
        record = update_or_create_result_record_by_id(
            result, project.uid, None, 'scenarios')
        db.session.add(record)
        db.session.commit()
    return make_mpld3_graph_dict(result)


def save_scenario_summaries(project_id, scenario_summaries):
    delete_result_by_parset_id(project_id, None, "scenarios")
    project_record = load_project_record(project_id)
    project = project_record.load()
    set_scenario_summaries_on_project(project, scenario_summaries)
    project_record.save_obj(project)
    return {'scenarios': get_scenario_summaries(project)}


def load_and_resolve_project(project_id):
    project_record = load_project_record(project_id)
    project = project_record.load()
    if resolve_project(project):
        print(">> Resolved project updated")
        project_record.save_obj(project)
    return project


def load_scenario_summaries(project_id):
    project = load_and_resolve_project(project_id)
    return {
        'scenarios': get_scenario_summaries(project),
        'ykeysByParsetId': get_parameters_for_scenarios(project),
        'defaultBudgetsByProgsetId': parse.get_budgets_for_scenarios(project),
        'defaultCoveragesByParsetIdyProgsetId': parse.get_coverages_for_scenarios(project),
        'years': get_project_years(project)
    }


## OPTIMIZATION

def load_optimization_summaries(project_id):
    project = load_and_resolve_project(project_id)
    return {
        'optimizations': get_optimization_summaries(project),
        'defaultOptimizationsByProgsetId': get_default_optimization_summaries(project)
    }


def save_optimization_summaries(project_id, optimization_summaries):
    project_record = load_project_record(project_id)
    project = project_record.load()
    old_names = [o.name for o in project.optims.values()]
    set_optimization_summaries_on_project(project, optimization_summaries)
    new_names = [o.name for o in project.optims.values()]
    deleted_names = [name for name in old_names if name not in new_names]
    deleted_result_names = ['optim-' + name for name in deleted_names]
    for result_name in deleted_result_names:
        delete_result_by_name(project.uid, result_name)
    project_record.save_obj(project)
    return {'optimizations': get_optimization_summaries(project)}


def upload_optimization_summary(project_id, optimization_id, optimization_summary):
    project_record = load_project_record(project_id)
    project = project_record.load()
    old_optim = get_optimization_from_project(project, optimization_id)
    optimization_summary['id'] = optimization_id
    optimization_summary['name'] = old_optim.name
    set_optimization_summaries_on_project(project, [optimization_summary])
    project_record.save_obj(project)
    return {'optimizations': get_optimization_summaries(project)}


def load_optimization_graphs(project_id, optimization_id, which):
    project = load_project(project_id)
    optimization = get_optimization_from_project(project, optimization_id)
    result = load_result_by_optimization(project, optimization)
    if result is None:
        return {}
    else:
        print(">> Loading graphs for result '%s'" % result.name)
        return make_mpld3_graph_dict(result, which)


## SPREADSHEETS

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


def load_data_spreadsheet_binary(project_id):
    data_record = ProjectDataDb.query.get(project_id)
    if data_record is not None:
        binary = data_record.meta
        if len(binary.meta) > 0:
            project = load_project(project_id)
            server_fname = secure_filename('{}.xlsx'.format(project.name))
            return server_fname, binary
    return None, None


def load_template_data_spreadsheet(project_id):
    project = load_project(project_id)
    fname = secure_filename('{}.xlsx'.format(project.name))
    server_fname = templatepath(fname)
    op.makespreadsheet(
        server_fname,
        pops=get_populations_from_project(project),
        datastart=int(project.data["years"][0]),
        dataend=int(project.data["years"][-1]))
    return upload_dir_user(TEMPLATEDIR), fname


def load_econ_spreadsheet_binary(project_id):
    econ_record = ProjectEconDb.query.get(project_id)
    if econ_record is not None:
        binary = econ_record.meta
        if len(binary.meta) > 0:
            project = load_project(project_id)
            server_fname = secure_filename('{}_economics.xlsx'.format(project.name))
            return server_fname, binary
    return None, None


def load_template_econ_spreadsheet(project_id):
    project = load_project(project_id)
    fname = secure_filename('{}_economics.xlsx'.format(project.name))
    server_fname = templatepath(fname)
    op.makeeconspreadsheet(
        server_fname,
        datastart=int(project.data["years"][0]),
        dataend=int(project.data["years"][-1]))
    return upload_dir_user(TEMPLATEDIR), fname


def update_project_from_data_spreadsheet(project_id, full_filename):
    project_record = load_project_record(project_id, raise_exception=True)
    project = project_record.load()

    parset_name = "default"
    parset_names = project.parsets.keys()
    basename = os.path.basename(full_filename)
    if parset_name in parset_names:
        parset_name = "uploaded from " + basename
        i = 0
        while parset_name in parset_names:
            i += 1
            parset_name = "uploaded_from_%s (%d)" % (basename, i)

    project.loadspreadsheet(full_filename, parset_name, makedefaults=True)
    project_record.save_obj(project)

    # importing spreadsheet will also runsim to project.results[-1]
    result = project.results[-1]
    parset = project.parsets[parset_name]
    result_record = update_or_create_result_record_by_id(
        result, project_id, parset.uid, "calibration")

    # save the binary data of spreadsheet for later download
    with open(full_filename, 'rb') as f:
        try:
            data_record = ProjectDataDb.query.get(project_id)
            data_record.meta = f.read()
        except:
            data_record = ProjectDataDb(project_id, f.read())

    db.session.add(data_record)
    db.session.add(project_record)
    db.session.add(result_record)
    db.session.commit()


def update_project_from_econ_spreadsheet(project_id, econ_spreadsheet_fname):
    project_record = load_project_record(project_id, raise_exception=True)
    project = project_record.load()

    project.loadeconomics(econ_spreadsheet_fname)
    project_record.save_obj(project)
    db.session.add(project_record)

    with open(econ_spreadsheet_fname, 'rb') as f:
        binary = f.read()
        upload_time = datetime.now(dateutil.tz.tzutc())
        econ_record = ProjectEconDb.query.get(project.id)
        if econ_record is not None:
            econ_record.meta = binary
            econ_record.updated = upload_time
        else:
            econ_record = ProjectEconDb(
                project_id=project.id,
                meta=binary,
                updated=upload_time)
        db.session.add(econ_record)

    db.session.commit()

    return econ_spreadsheet_fname


def delete_econ(project_id):
    project_record = load_project_record(project_id)
    project = project_record.load()
    if 'econ' not in project.data:
        raise Exception("No economoics data found in project %s" % project_id)
    del project.data['econ']
    project_record.save_obj(project)
    db.session.add(project_record)

    econ_record = ProjectEconDb.query.get(project.id)
    if econ_record is None or len(econ_record.meta) == 0:
        db.session.delete(econ_record)
    else:
        raise Exception("No economics data has been uploaded")

    db.session.commit()


def get_odict_item(odict, key):
    if type(key) == int:
        if 0 <= key < len(odict):
            return odict[key]
    elif key in odict:
        return odict[key]
    return None


def resolve_project(project):
    """
    This function does some checks on the project to ensure that all the cross-reference
    fields are properly specified and that defaults are sensibly populated.

    Returns: boolean to whether any changes needed to be made to the project
    """
    print(">> Resolve project")
    is_change = False

    del_scenario_keys = []
    for scenario_key, scenario in project.scens.items():
        print(">> Checking scenario", scenario_key)
        if type(scenario.parsetname) is int:
            i = scenario.parsetname
            try:
                scenario.parsetname = project.parsets[i].name
                is_change = True
            except:
                del_scenario_keys.append(scenario_key)
        elif scenario.parsetname not in project.parsets:
            del_scenario_keys.append(scenario_key)
        if hasattr(scenario, "progsetname"):
            if type(scenario.progsetname) is int:
                i = scenario.progsetname
                try:
                    scenario.progsetname = project.progsets[i].name
                    is_change = True
                except:
                    del_scenario_keys.append(scenario_key)
            if scenario.progsetname not in project.progsets:
                del_scenario_keys.append(scenario_key)
    print(">> Delete deprecated scenarios", del_scenario_keys)
    for scenario_key in del_scenario_keys:
        del project.scens[scenario_key]

    is_change = is_change or len(del_scenario_keys) > 0

    # makes sure there is a parset called default as defaultprograms requires this
    if "default" not in project.parsets and len(project.parsets) > 0:
        parsetname = project.parsets[0].name
        project.copyparset(orig=parsetname, new="default")
        is_change = True

    # check optimizations are good
    del_optim_keys = []
    for optim_key, optim in project.scens.items():
        print(">> Checking optim", optim_key)
        if type(optim.parsetname) is int:
            i = optim.parsetname
            try:
                optim.parsetname = project.parsets[i].name
                is_change = True
            except:
                del_optim_keys.append(optim_key)
        elif optim.parsetname not in project.parsets:
            del_optim_keys.append(optim_key)
        if hasattr(optim, "progsetname"):
            if type(optim.progsetname) is int:
                i = optim.progsetname
                try:
                    optim.progsetname = project.progsets[i].name
                    is_change = True
                except:
                    del_optim_keys.append(optim_key)
            if optim.progsetname not in project.progsets:
                del_optim_keys.append(optim_key)
    print(">> Delete deprecated optims", del_optim_keys)
    for optim_key in del_optim_keys:
        del project.optimis[optim_key]

    is_change = is_change or len(del_optim_keys) > 0

    # ensure constraints set to None are given a default
    for optim in project.optims.values():
        progset_name = optim.progsetname
        progset = project.progsets[progset_name]
        if optim.constraints is None:
            print(">> Fill out default constraints for constraints = None")
            optim.constraints = op.defaultconstraints(project=project, progset=progset)
            is_change = True

    results = db.session.query(ResultsDb).filter_by(project_id=project.uid)
    parset_ids = [parset.uid for parset in project.parsets.values()]
    is_delete_result = False
    for result in results:
        if result.parset_id is not None and result.parset_id not in parset_ids:
            print(">> Delete deprecated result", result.parset_id)
            db.session.delete(result)
            is_delete_result = True
    db.session.commit()

    is_change = is_change or is_delete_result

    return is_change



## PROGRAMS


def load_target_popsizes(project_id, parset_id, progset_id, program_id):
    project = load_project(project_id)
    parset = get_parset_from_project(project, parset_id)
    progset = get_progset_from_project(project, progset_id)
    program = get_program_from_progset(progset, program_id)
    years = get_project_years(project)
    popsizes = program.gettargetpopsize(t=years, parset=parset)
    return normalize_obj(dict(zip(years, popsizes)))


def load_project_program_summaries(project_id):
    project = load_project(project_id, raise_exception=True)
    return get_default_program_summaries(project)


def load_progset_summary(project_id, progset_id):
    project = load_project(project_id)
    progset = get_progset_from_project(project, progset_id)
    return get_progset_summary(project, progset.name)


def load_progset_summaries(project_id):
    project = load_project(project_id)
    return get_progset_summaries(project)


def create_progset(project_id, progset_summary):
    project_record = load_project_record(project_id)
    project = project_record.load()
    set_progset_summary_on_project(project, progset_summary)
    project_record.save_obj(project)
    return get_progset_summary(project, progset_summary["name"])


def save_progset(project_id, progset_id, progset_summary):
    project_record = load_project_record(project_id)
    project = project_record.load()
    set_progset_summary_on_project(project, progset_summary, progset_id=progset_id)
    project_record.save_obj(project)
    return get_progset_summary(project, progset_summary["name"])


def upload_progset(project_id, progset_id, progset_summary):
    project_record = load_project_record(project_id)
    project = project_record.load()
    old_progset = get_progset_from_project(project, progset_id)
    print(">> Upload progset '%s' into '%s'" % (progset_summary['name'], old_progset.name))
    progset_summary['id'] = progset_id
    progset_summary['name'] = old_progset.name
    set_progset_summary_on_project(project, progset_summary, progset_id=progset_id)
    project_record.save_obj(project)
    return get_progset_summary(project, progset_summary["name"])


def copy_progset(project_id, progset_id, new_progset_name):

    def update_project_fn(project):
        original_progset = get_progset_from_project(project, progset_id)
        project.copyprogset(orig=original_progset.name, new=new_progset_name)
        project.progsets[new_progset_name].uid = op.uuid()

    update_project_with_fn(project_id, update_project_fn)
    return load_progset_summaries(project_id)


def delete_progset(project_id, progset_id):
    project_record = load_project_record(project_id)
    project = project_record.load()

    progset = get_progset_from_project(project, progset_id)

    progset_name = progset.name
    optims = [o for o in project.optims.values() if o.progsetname == progset_name]

    for optim in optims:
        result_name = 'optim-' + optim.name
        delete_result_by_name(project.uid, result_name)
        project.optims.pop(optim.name)

    project.progsets.pop(progset.name)

    project_record.save_obj(project)


def load_progset_outcome_summaries(project_id, progset_id):
    project = load_project(project_id)
    progset = get_progset_from_project(project, progset_id)
    outcomes = get_outcome_summaries_from_progset(progset)
    return outcomes


def save_outcome_summaries(project_id, progset_id, outcome_summaries):
    project_record = load_project_record(project_id)
    project = project_record.load()
    progset = get_progset_from_project(project, progset_id)
    set_outcome_summaries_on_progset(outcome_summaries, progset)
    project_record.save_obj(project)
    return get_outcome_summaries_from_progset(progset)


def save_program(project_id, progset_id, program_summary):
    project_record = load_project_record(project_id)
    project = project_record.load()

    progset = get_progset_from_project(project, progset_id)

    print("> Saving program " + program_summary['name'])
    set_program_summary_on_progset(progset, program_summary)

    progset.updateprogset()

    project_record.save_obj(project)


def load_costcov_graph(project_id, progset_id, program_id, parset_id, t):
    project_record = load_project_record(project_id)
    project = project_record.load()
    progset = get_progset_from_project(project, progset_id)

    program = get_program_from_progset(progset, program_id)
    plotoptions = None
    if hasattr(program, "attr"):
        plotoptions = program.attr

    parset = get_parset_from_project(project, parset_id)
    plot = program.plotcoverage(t=t, parset=parset, plotoptions=plotoptions)

    return convert_to_mpld3(plot)



import datetime
import dateutil.tz

from flask.ext.sqlalchemy import SQLAlchemy
from sqlalchemy.orm import sessionmaker, scoped_session

import optima as op

from server.api import app
from server.webapp.dbmodels import WorkLogDb, WorkingProjectDb
from server.webapp.exceptions import ProjectDoesNotExist
from server.webapp.dataio import save_result, load_project_record, load_parset_record

from celery import Celery


db = SQLAlchemy(app)

celery_instance = Celery(app.import_name, broker=app.config['CELERY_BROKER_URL'])
celery_instance.conf.update(app.config)

TaskBase = celery_instance.Task


class ContextTask(TaskBase):
    abstract = True

    def __call__(self, *args, **kwargs):
        with app.app_context():
            return TaskBase.__call__(self, *args, **kwargs)


celery_instance.Task = ContextTask


def init_db_session():
    """
    Create scoped_session, eventually bound to engine
    """
    return scoped_session(sessionmaker(db.engine))


def close_db_session(db_session):
    # this line might be redundant (not 100% sure - not clearly described)
    db_session.connection().close() # pylint: disable=E1101
    db_session.remove()
    # black magic to actually close the connection by forcing the engine to dispose of garbage (I assume)
    db_session.bind.dispose() # pylint: disable=E1101


def start_or_report_calculation(project_id, parset_id, work_type):
    print "start_or_report_calculation(%s %s %s)" % (project_id, parset_id, work_type)
    # work type is supposed to correspond to the method of Project class
    # db_session = init_db_session() this does not work
    db_session = init_db_session()
    can_start = False
    can_join = False
    wp_parset_id = parset_id
    project_entry = load_project_record(project_id, raise_exception=False, db_session=db_session)
    if not project_entry:
        close_db_session(db_session)
        raise ProjectDoesNotExist(project_id)
    project_instance = project_entry.hydrate()
    project_pickled = op.saves(project_instance)
    wp = db_session.query(WorkingProjectDb).filter_by(id=project_entry.id).first()
    if wp is not None and wp.is_working and (wp.parset_id != parset_id or wp.work_type != work_type):
        print("wp is already present for %s" % project_entry.id)
        work_type = wp.work_type
    else:
        work_log = WorkLogDb(project_id=project_entry.id, parset_id=parset_id, work_type = work_type)
        work_log.start_time = datetime.datetime.now(dateutil.tz.tzutc())
        db_session.add(work_log)
        db_session.flush()
        if wp is None:
            print("No working project was found - creating new one")
            db_session.add(WorkingProjectDb(project_entry.id, parset_id = parset_id,
                project = project_pickled,
                is_working = True, work_type = work_type, work_log_id = work_log.id))
            can_start = True
            can_join = True
        else:
            print("Found working project for %s: %s %s %s" % (wp.id, wp.work_type, wp.parset_id, wp.is_working))
            can_start = not wp.is_working
            can_join = (not wp.is_working) or (wp.parset_id == parset_id and wp.work_type == work_type)
            if can_start:
                wp.work_type = work_type
                wp.parset_id = parset_id
                wp.is_working = True
                wp.project = project_pickled
                wp.work_log_id = work_log.id
                db_session.add(wp)
            else:
                wp_parset_id = wp.parset_id
                work_type = wp.work_type
    db_session.commit()
    close_db_session(db_session)
    return {
        'can_start': can_start,
        'can_join': can_join,
        'parset_id': wp_parset_id,
        'work_type': work_type
    }


def check_calculation_status(project_id):
    db_session = init_db_session()
    wp = db_session.query(WorkingProjectDb).get(project_id)
    work_log = db_session.query(WorkLogDb).get(wp.work_log_id)
    close_db_session(db_session)
    if work_log is not None:
        result = {
            'status': work_log.status,
            'error_text': work_log.error,
            'start_time': work_log.start_time,
            'stop_time': work_log.stop_time,
            'result_id': work_log.result_id
        }
    else:
        result = {
            'status': 'unknown',
            'error_text': None,
            'start_time': None,
            'stop_time': None,
            'result_id': None
        }
    return result


@celery_instance.task()
def run_autofit(project_id, parset_id, maxtime=60):
    import traceback
    error_text = ""
    status = 'completed'
    db_session = init_db_session()
    wp = db_session.query(WorkingProjectDb).filter_by(id=project_id).first()
    project_instance = op.loads(wp.project)
    parset_name = None
    for key in project_instance.parsets:
        parset = project_instance.parsets[key]
        if str(parset.uid) == str(parset_id):
            parset_name = parset.name
            break
    app.logger.debug("started autofit: {} {}".format(project_id, parset_name))
    close_db_session(db_session)
    result = None
    try:
        project_instance.autofit(
            name=str(parset_name),
            orig=str(parset_name),
            maxtime=maxtime
        )
        result = project_instance.parsets[str(parset_name)].getresults()
        print "result", result
    except Exception:
        var = traceback.format_exc()
        print("ERROR for project_id: %s, args: %s calculation: %s\n %s" % (project_id, parset_name, 'autofit', var))
        error_text = var
        status='error'

    db_session = init_db_session()
    wp = db_session.query(WorkingProjectDb).filter_by(id=project_id).first()
    wp.project = op.saves(project_instance)
    work_log = db_session.query(WorkLogDb).get(wp.work_log_id)
    work_log.status = status
    work_log.error = error_text
    work_log.stop_time = datetime.datetime.now(dateutil.tz.tzutc())
    if result:
        result_entry = save_result(project_id, result, parset_name, 'autofit', db_session=db_session)
        db_session.add(result_entry)
        db_session.flush()
        work_log.result_id = result_entry.id
    db_session.add(work_log)
    wp.is_working = False
    wp.work_type = None
    db_session.add(wp)
    db_session.commit()
    close_db_session(db_session)
    app.logger.debug("stopped autofit")


@celery_instance.task()
def run_optimization(project_id, optimization_name, parset_name, progset_name, objectives, constraints):
    import traceback
    app.logger.debug('started optimization: {} {} {} {} {} {}'.format(
        project_id, optimization_name, parset_name, progset_name, objectives, constraints))
    error_text = ""
    status = 'completed'
    db_session = init_db_session()
    wp = db_session.query(WorkingProjectDb).filter_by(id=project_id).first()
    project_instance = op.loads(wp.project)
    close_db_session(db_session)
    result = None
    if not objectives['budget']:
        objectives['budget'] = 1000000
    try:
        result = project_instance.optimize(
            name=optimization_name,
            parsetname=parset_name,
            progsetname=progset_name,
            objectives=objectives,
            constraints=constraints
        )
        # result = op.defaults.defaultproject('generalized').optimize()  # TODO replace with the real thing late
        print "result", result
    except Exception:
        var = traceback.format_exc()
        print("ERROR for project_id: %s, args: %s calculation: %s\n %s" % (project_id, optimization_name, 'optimization', var))
        error_text = var
        status='error'

    db_session = init_db_session()
    wp = db_session.query(WorkingProjectDb).filter_by(id=project_id).first()
    wp.project = op.saves(project_instance)
    work_log = db_session.query(WorkLogDb).get(wp.work_log_id)
    work_log.status = status
    work_log.error = error_text
    work_log.stop_time = datetime.datetime.now(dateutil.tz.tzutc())
    if result:
        result_entry = save_result(project_id, result, parset_name, 'optimization', db_session=db_session)
        db_session.add(result_entry)
        db_session.flush()
        work_log.result_id = result_entry.id
    db_session.add(work_log)
    wp.is_working = False
    wp.work_type = None
    db_session.add(wp)
    db_session.commit()
    close_db_session(db_session)
    app.logger.debug("stopped optimization")

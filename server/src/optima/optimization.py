#!/bin/env python
# -*- coding: utf-8 -*-
"""
Optimization Module
~~~~~~~~~~~~~~~~~~~

1. Start an optimization request
2. Stop an optimization
3. Save current optimization model
4. Revert to the last saved model
"""
from flask import request, jsonify, Blueprint, current_app
from dbconn import db
from async_calculate import CalculatingThread, start_or_report_calculation, cancel_calculation, check_calculation
from async_calculate import check_calculation_status, good_exit_status
from utils import check_project_name, project_exists, load_model, save_model, \
revert_working_model_to_default, save_working_model_as_default, report_exception
from sim.optimize import optimize
from sim.bunch import bunchify, unbunchify
import json
import traceback
from flask.ext.login import login_required, current_user

# route prefix: /api/analysis/optimization
optimization = Blueprint('optimization',  __name__, static_folder = '../static')

def get_optimization_results(D_dict):
    return {'graph': D_dict.get('plot',{}).get('OM',{}), 'pie':D_dict.get('plot',{}).get('OA',{})}

@optimization.route('/list')
@login_required
@check_project_name
@report_exception()
def getOptimizationParameters():
    """ retrieve list of optimizations defined by the user, with parameters """
    from sim.optimize import defaultoptimizations
    current_app.logger.debug("/api/analysis/optimization/list")
    # get project name
    project_name = request.project_name
    if not project_exists(project_name):
        reply['reason'] = 'Project %s does not exist' % project_name
        return reply
    D = load_model(project_name)
    if not 'optimizations' in D:
        optimizations = defaultoptimizations(D)
    else:
        optimizations = D.optimizations
    optimizations = unbunchify(optimizations)
    return json.dumps({'optimizations':optimizations})


@optimization.route('/start', methods=['POST'])
@login_required
@check_project_name
def startOptimization():
    """ Start optimization """
    data = json.loads(request.data)
    current_app.logger.debug("optimize: %s" % data)
    # get project name
    project_name = request.project_name
    if not project_exists(project_name):
        from utils import BAD_REPLY
        reply = BAD_REPLY
        reply['reason'] = 'File for project %s does not exist' % project_name
        return jsonify(reply)
    try:
        can_start, can_join, current_calculation = start_or_report_calculation(current_user.id, project_name, optimize, db.session)
        if can_start:
            # Prepare arguments
            args = {'verbose':0}
            objectives = data.get('objectives')
            if objectives:
                args['objectives'] = bunchify( objectives )
            constraints = data.get('constraints')
            if constraints:
                args['constraints'] = bunchify( constraints )
            timelimit = int(data.get("timelimit")) # for the thread
            args["timelimit"] = 10 # for the autocalibrate function
            CalculatingThread(db.engine, current_user, project_name, timelimit, optimize, args).start()
            msg = "Starting optimization thread for user %s project %s" % (current_user.name, project_name)
            current_app.logger.debug(msg)
            return json.dumps({"status":"OK", "result": msg, "join":True})
        else:
            msg = "Thread for user %s project %s (%s) has already started" % (current_user.name, project_name, current_calculation)
            return json.dumps({"status":"OK", "result": msg, "join":can_join})
    except Exception, err:
        var = traceback.format_exc()
        return jsonify({"status":"NOK", "exception":var})

@optimization.route('/stop')
@login_required
@check_project_name
def stopCalibration():
    """ Stops calibration """
    prj_name = request.project_name
    cancel_calculation(current_user.id, prj_name, optimize, db.session)
    return json.dumps({"status":"OK", "result": "optimize calculation for user %s project %s requested to stop" % (current_user.name, prj_name)})

@optimization.route('/working')
@login_required
@check_project_name
@report_exception()
def getWorkingModel():
    """ Returns the working model for optimization. """
    D_dict = {}
    # Get optimization working data
    prj_name = request.project_name
    error_text = None
    if check_calculation(current_user.id, prj_name, optimize, db.session):
        D_dict = load_model(prj_name, working_model = True, as_bunch = False)
        status = 'Running'
    else:
        current_app.logger.debug("no longer optimizing")
        status, error_text = check_calculation_status(current_user.id, prj_name, optimize, db.session)
        if status in good_exit_status:
            status = 'Done'
        else:
            status = 'NOK'
    result = get_optimization_results(D_dict)
    result['status'] = status
    if error_text:
        result['exception'] = error_text
    return jsonify(result)

@optimization.route('/save', methods=['POST'])
@login_required
@check_project_name
def saveModel():
    from sim.optimize import saveoptimization
    """ Saves working model as the default model """
    reply = {'status':'NOK'}
    data = json.loads(request.data)
    name = data.get('name', 'default')
    objectives = data.get('objectives')
    if objectives:objectives = bunchify( objectives )
    constraints = data.get('constraints')
    if constraints:constraints = bunchify( constraints )

    # get project name
    project_name = request.project_name
    if not project_exists(project_name):
        reply['reason'] = 'File for project %s does not exist' % project_name
        return jsonify(reply)

    try:
        # read the previously saved optmizations
        D_dict = load_model(project_name, as_bunch = False)
        optimizations = D_dict.get('optimizations')

        # now, save the working model, read results and save for the optimization with the given name
        D_dict = save_working_model_as_default(project_name)
        result = get_optimization_results(D_dict)
        D = bunchify(D_dict)

        # get the previously stored optimizations
        if optimizations: D.optimizations = bunchify(optimizations)

        #save the optimization parameters
        D = saveoptimization(D, name, objectives, constraints, result)
        D_dict = D.toDict()
        save_model(project_name, D_dict)

        reply['status']='OK'
        reply['optimizations'] = D_dict['optimizations']
        return jsonify(reply)
    except Exception, err:
        reply['exception'] = traceback.format_exc()
        return jsonify(reply)

@optimization.route('/revert', methods=['POST'])
@login_required
@check_project_name
def revertCalibrationModel():
    """ Revert working model to the default model """
    reply = {'status':'NOK'}

    # get project name
    project_name = request.project_name
    if not project_exists(project_name):
        reply['reason'] = 'File for project %s does not exist' % project_name
        return jsonify(reply)
    try:
        revert_working_model_to_default(project_name)
        return jsonify({"status":"OK"})
    except Exception, err:
        reply['exception'] = traceback.format_exc()
        return jsonify(reply)

@optimization.route('/remove/<name>', methods=['POST'])
@login_required
@check_project_name
def removeOptimizationSet(name):
    """ Removes given optimization from the optimization set """
    from sim.optimize import removeoptimization
    reply = {'status':'NOK'}

    # get project name
    project_name = request.project_name
    if not project_exists(project_name):
        reply['reason'] = 'File for project %s does not exist' % project_name
        return jsonify(reply)
    else:
        D = load_model(project_name, as_bunch = True)
        D = removeoptimization(name)
        D_dict = D.toDict()
        save_model(project_name, D_dict)
        reply['status']='OK'
        reply['name'] = 'deleted'
        reply['optimizations'] = D_dict['optimizations']
        return jsonify(reply)
    




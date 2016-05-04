from datetime import datetime
import dateutil
import uuid
import json
import pprint

from flask import current_app, helpers, make_response, request

from flask import jsonify

from flask.ext.login import login_required
from flask_restful import Resource, marshal_with, abort
from flask_restful_swagger import swagger

from server.webapp.inputs import (SubParser, secure_filename_input, AllowedSafeFilenameStorage,
                                  Json as JsonInput)

from server.webapp.utils import (load_project_record, RequestParser, report_exception, TEMPLATEDIR,
                                 upload_dir_user, save_result)
from server.webapp.exceptions import ParsetDoesNotExist, ParsetAlreadyExists

from server.webapp.dbconn import db
from server.webapp.dbmodels import ParsetsDb, ResultsDb, WorkingProjectDb, ScenariosDb,OptimizationsDb
from server.webapp.fields import Json, Uuid

import math

import optima as op


copy_parser = RequestParser()
copy_parser.add_arguments({
    'name': {'required': True},
    'parset_id': {'type': uuid.UUID}
})

y_keys_fields = {
    'keys': Json
}

limits_fields = {
    'parsets': Json
}


class ParsetYkeys(Resource):

    @swagger.operation(
        summary='get parsets ykeys'
    )
    @marshal_with(y_keys_fields)
    def get(self, project_id):
        project_entry = load_project_record(project_id, raise_exception=True)

        reply = db.session.query(ParsetsDb).filter_by(project_id=project_entry.id).all()
        parsets = {str(item.id): item.hydrate() for item in reply}
        y_keys = {
           id: {par.short: [{
                    'val': k,
                    'label': ' - '.join(k) if isinstance(k, tuple) else k
                } for k in par.y.keys()] for par in parset.pars[0].values()
           if hasattr(par, 'y') and par.visible}
           for id, parset in parsets.iteritems()
        }
        return {'keys': y_keys}


class ParsetLimits(Resource):
    @swagger.operation(
        summary='get parameters limits'
    )
    @marshal_with(limits_fields)
    def get(self, project_id):
        project_entry = load_project_record(project_id, raise_exception=True)
        be_project = project_entry.hydrate()

        reply = db.session.query(ParsetsDb).filter_by(project_id=project_entry.id).all()
        parsets = {str(item.id): item.hydrate() for item in reply}
        limits = {
           id: {par.short: [
                    be_project.settings.convertlimits(limits=limit) if isinstance(limit, str) else limit
                    for limit in par.limits
                ] for par in parset.pars[0].values()
           if hasattr(par, 'y') and par.visible}
           for id, parset in parsets.iteritems()
        }
        return {'parsets': limits}


class Parsets(Resource):
    """
    Parsets for a given project.
    """
    method_decorators = [report_exception, login_required]

    @swagger.operation(
        description='Download parsets for the project with the given id.',
        notes="""
            if project exists, returns parsets for it
            if project does not exist, returns an error.
        """,
        responseClass=ParsetsDb.__name__
    )

    @marshal_with(ParsetsDb.resource_fields, envelope='parsets')
    def get(self, project_id):

        current_app.logger.debug("/api/project/%s/parsets" % str(project_id))
        project_entry = load_project_record(project_id, raise_exception=True)
        reply = db.session.query(ParsetsDb).filter_by(project_id=project_entry.id).all()
        result = [item.hydrate() for item in reply]

        return result

    @swagger.operation(
        description='Create new parset with default settings or copy existing parset',
        notes="""
            If parset_id argument is given, copy from the existing parset.
            Otherwise, create a parset with default settings
            """
    )
    @report_exception
    def post(self, project_id):
        current_app.logger.debug("POST /api/project/{}/parsets".format(project_id))
        args = copy_parser.parse_args()
        print "args", args
        name = args['name']
        parset_id = args.get('parset_id')

        project_entry = load_project_record(project_id, raise_exception=True)
        project_instance = project_entry.hydrate()
        if name in project_instance.parsets:
            raise ParsetAlreadyExists(project_id, name)
        if not parset_id:
            # create new parset with default settings
            project_instance.makeparset(name, overwrite=False)
            new_result = project_instance.runsim(name)
            project_entry.restore(project_instance)
            db.session.add(project_entry)

            result_record = save_result(project_entry.id, new_result, name)
            db.session.add(result_record)
        else:
            # dealing with uid's directly might be messy...
            original_parset = [item for item in project_entry.parsets if item.id == parset_id]
            if not original_parset:
                raise ParsetDoesNotExist(parset_id, project_id=project_id)
            original_parset = original_parset[0]
            parset_name = original_parset.name
            project_instance.copyparset(orig=parset_name, new=name)
            project_entry.restore(project_instance)
            db.session.add(project_entry)

            old_result_record = db.session.query(ResultsDb).filter_by(
                parset_id=str(parset_id), project_id=str(project_id),
                calculation_type=ResultsDb.CALIBRATION_TYPE).first()
            old_result = old_result_record.hydrate()
            new_result = op.dcp(old_result)
            new_result_record = save_result(project_entry.id, new_result, name)
            db.session.add(new_result_record)

        db.session.commit()

        rv = []
        for item in project_entry.parsets:
            rv_item = item.hydrate().__dict__
            rv_item['id'] = item.id
            rv.append(rv_item)

        return rv


rename_parser = RequestParser()
rename_parser.add_argument('name', required=True)


class ParsetsDetail(Resource):
    """
    Single Parset.
    """
    method_decorators = [report_exception, login_required]

    @swagger.operation(
        description='Delete parset with the given id.',
        notes="""
            if parset exists, delete it
            if parset does not exist, returns an error.
        """
    )
    @report_exception
    @marshal_with(ParsetsDb.resource_fields, envelope='parsets')
    def delete(self, project_id, parset_id):

        current_app.logger.debug("DELETE /api/project/{}/parsets/{}".format(project_id, parset_id))
        project_entry = load_project_record(project_id, raise_exception=True)

        parset = db.session.query(ParsetsDb).filter_by(project_id=project_entry.id, id=parset_id).first()
        if parset is None:
            raise ParsetDoesNotExist(id=parset_id, project_id=project_id)

        # # Is this how we should check for default parset?
        # if parset.name.lower() == 'default':  # TODO: it is lowercase
        #     abort(403)

        # TODO: also delete the corresponding calibration results
        db.session.query(ResultsDb).filter_by(project_id=project_id,
            id=parset_id, calculation_type=ResultsDb.CALIBRATION_TYPE).delete()
        db.session.query(ScenariosDb).filter_by(project_id=project_id,
            parset_id=parset_id).delete()
        db.session.query(OptimizationsDb).filter_by(project_id=project_id,
            parset_id=parset_id).delete()
        db.session.query(ParsetsDb).filter_by(project_id=project_id, id=parset_id).delete()
        db.session.commit()

        return '', 204

    @swagger.operation(
        description='Rename parset with the given id',
        notes="""
            if parset exists, rename it
            if parset does not exist, return an error.
            """
    )
    @report_exception
    @marshal_with(ParsetsDb.resource_fields, envelope='parsets')
    def put(self, project_id, parset_id):
        """
        For consistency, let's always return the updated parsets for operations on parsets
        (so that FE doesn't need to perform another GET call)
        """

        current_app.logger.debug("PUT /api/project/{}/parsets/{}".format(project_id, parset_id))
        args = rename_parser.parse_args()
        name = args['name']

        project_entry = load_project_record(project_id, raise_exception=True)
        target_parset = [item for item in project_entry.parsets if item.id == parset_id]
        if target_parset:
            target_parset = target_parset[0]
        if not target_parset:
            raise ParsetDoesNotExist(id=parset_id, project_id=project_id)
        target_parset.name = name
        db.session.add(target_parset)
        db.session.commit()
        return [item.hydrate() for item in project_entry.parsets]


calibration_fields = {
    "parset_id": Uuid,
    "parameters": Json,
    "graphs": Json,
    "selectors": Json,
    "result_id": Uuid,
}

calibration_parser = RequestParser()
calibration_parser.add_argument('which', location='args', default=None, action='append')
calibration_parser.add_argument('autofit', location='args', default=False, type=bool)


calibration_update_parser = RequestParser()
calibration_update_parser.add_arguments({
    'which': {'default': None, 'action': 'append'},
    'parameters': {'required': True, 'type': dict, 'action': 'append'},
    'doSave': {'default': False, 'type': bool, 'location': 'args'},
    'autofit': {'default': False, 'type': bool, 'location': 'args'}
})


parset_save_with_autofit_parser = RequestParser()
parset_save_with_autofit_parser.add_arguments({
    'parameters': {'type': JsonInput, 'required': True, 'action': 'append'},
    'result_id': {'type': uuid.UUID, 'required': True},
})


def get_parset_parameters(parset, ind=0):
    parameters = []
    for key, par in parset.pars[ind].items():
        if hasattr(par, 'fittable') and par.fittable != 'no':
            if par.fittable == 'meta':
                parameters.append({
                    "key": key,
                    "subkey": None,
                    "type": par.fittable,
                    "value": par.m,
                    "label": '%s -- meta' % par.name,
                })
            elif par.fittable == 'const':
                parameters.append({
                    "key": key,
                    "subkey": None,
                    "type": par.fittable,
                    "value": par.y,
                    "label": par.name,
                })
            elif par.fittable in ['pop', 'pship']:
                for subkey in par.y.keys():
                    parameters.append({
                        "key": key,
                        "subkey": subkey,
                        "type": par.fittable,
                        "value": par.y[subkey],
                        "label": '%s -- %s' % (par.name, str(subkey)),
                    })
            elif par.fittable == 'exp':
                for subkey in par.p.keys():
                    parameters.append({
                        "key": key,
                        "subkey": subkey,
                        "type": par.fittable,
                        "value": par.p[subkey][0],
                        "label": '%s -- %s' % (par.name, str(subkey)),
                    })
            else:
                print('Parameter type "%s" not implemented!' % par.fittable)
    return parameters


def put_parameters_in_parset(parameters, parset, ind=0):
    pars = parset.pars[ind]
    for p_dict in parameters:
        key = p_dict['key']
        value = p_dict['value']
        subkey = p_dict['subkey']
        par_type = p_dict['type']
        value = float(value)
        if par_type == 'meta':  # Metaparameters
            pars[key].m = value
        elif par_type in ['pop', 'pship']:  # Populations or partnerships
            pars[key].y[subkey] = value
        elif par_type == 'exp':  # Population growth
            pars[key].p[subkey][0] = value
        elif par_type == 'const':  # Metaparameters
            pars[key].y = value
        else:
            print('Parameter type "%s" not implemented!' % par_type)


class ParsetsCalibration(Resource):
    """
    Calibration info for the Parset.
    """

    method_decorators = [report_exception, login_required]

    def _result_to_jsons(self, result, which):
        import mpld3
        import json
        graphs = op.plotting.makeplots(result, figsize=(4, 3), toplot=[str(w) for w in which])  # TODO: store if that becomes an efficiency issue
        jsons = []
        for graph in graphs:
            # Add necessary plugins here
            mpld3.plugins.connect(graphs[graph], mpld3.plugins.MousePosition(fontsize=14, fmt='.4r'))
            # a hack to get rid of NaNs, javascript JSON parser doesn't like them
            json_string = json.dumps(mpld3.fig_to_dict(graphs[graph])).replace('NaN', 'null')
            jsons.append(json.loads(json_string))
        return jsons

    def _selectors_from_result(self, result, which):
        graph_selectors = op.getplotselections(result)
        keys = graph_selectors['keys']
        names = graph_selectors['names']
        if which is None:
            checks = graph_selectors['defaults']
        else:
            checks = [key in which for key in keys]
        selectors = [{'key': key, 'name': name, 'checked': checked}
                     for (key, name, checked) in zip(keys, names, checks)]
        return selectors

    def _which_from_selectors(self, graph_selectors):
        return [item['key'] for item in graph_selectors if item['checked']]

    @swagger.operation(
        description='Provides calibration information for the given parset',
        notes="""
        Returns data suitable for manual calibration and the set of corresponding graphs.
        """,
        parameters=calibration_parser.swagger_parameters()
    )
    @report_exception
    @marshal_with(calibration_fields, envelope="calibration")
    def get(self, project_id, parset_id):
        current_app.logger.debug("/api/project/{}/parsets/{}/calibration".format(project_id, parset_id))
        args = calibration_parser.parse_args()
        which = args.get('which')
        autofit = args.get('autofit', False)
        print "is autofit", autofit

        if not autofit:
            project_record = load_project_record(project_id, raise_exception=True)
            project = project_record.hydrate()
        else: # todo bail out if no working project
            wp = db.session.query(WorkingProjectDb).filter_by(id=project_id).first()
            project = op.loads(wp.project)
        
        parset = [project.parsets[item] for item in project.parsets
            if project.parsets[item].uid == parset_id]
        if not parset:
            raise ParsetDoesNotExist(project_id=project_id, id=parset_id)
        else:
            parset = parset[0]

        # store simulation results for plotting
        calculation_type = 'autofit' if autofit else ResultsDb.CALIBRATION_TYPE
        result_record = db.session.query(ResultsDb).filter_by(
            project_id=project_id, parset_id=parset_id, calculation_type=calculation_type).first()
        if result_record:
            result = result_record.hydrate()
        else:
            simparslist = parset.interp()
            result = project.runsim(simpars=simparslist)

        # generate graphs
        selectors = self._selectors_from_result(result, which)
        which = which or self._which_from_selectors(selectors)
        graphs = self._result_to_jsons(result, which)

        return {
            "parset_id": parset_id,
            "parameters": get_parset_parameters(parset),
            "graphs": graphs,
            "selectors": selectors,
            "result_id": result_record.id if result_record else None
        }

    @report_exception
    @marshal_with(calibration_fields, envelope="calibration")
    def put(self, project_id, parset_id):
        current_app.logger.debug("PUT /api/project/{}/parsets/{}/calibration".format(project_id, parset_id))
        args = calibration_update_parser.parse_args()
        parameters = args.get('parameters', [])
        which = args.get('which')
        doSave = args.get('doSave')
        autofit = args.get('autofit', False)

        parset_record = db.session.query(ParsetsDb).filter_by(id=parset_id).first()
        if parset_record is None or parset_record.project_id!=project_id:
            raise ParsetDoesNotExist(id=parset_id)

        # save parameters
        parset = parset_record.hydrate()
        put_parameters_in_parset(parameters, parset)

        # recalculate
        project_record = load_project_record(parset_record.project_id, raise_exception=True)
        project = project_record.hydrate()
        simparslist = parset.interp()
        result = project.runsim(simpars=simparslist)
        result_record = None

        if doSave:  # save the updated results
            parset_record.pars = op.saves(parset.pars)
            parset_record.updated = datetime.now(dateutil.tz.tzutc())
            db.session.add(parset_record)
            result_record = [item for item in project_record.results if
                            item.parset_id == parset_id and item.calculation_type == ResultsDb.CALIBRATION_TYPE]
            if result_record:
                result_record = result_record[-1]
                result_record.blob = op.saves(result)
            else:
                result_record = ResultsDb(
                    parset_id=parset_id,
                    project_id=project_record.id,
                    calculation_type=ResultsDb.CALIBRATION_TYPE,
                    blob=op.saves(result)
                )
            db.session.add(result_record)
            db.session.commit()

        # generate graphs
        selectors = self._selectors_from_result(result, which)
        which = which or self._which_from_selectors(selectors)
        graphs = self._result_to_jsons(result, which)

        return {
            "parset_id": parset_id,
            "parameters": get_parset_parameters(parset),
            "graphs": graphs,
            "selectors": selectors,
            "result_id": result_record.id if result_record is not None else None
        }

    @report_exception
    @marshal_with(calibration_fields, envelope="calibration")
    def post(self, project_id, parset_id):
        current_app.logger.debug("POST /api/project/{}/parsets/{}/calibration".format(project_id, parset_id))
        args = parset_save_with_autofit_parser.parse_args()
        parameters = args['parameters']

        parset_entry = db.session.query(ParsetsDb).filter_by(id=parset_id).first()
        if parset_entry is None or parset_entry.project_id != project_id:
            raise ParsetDoesNotExist(id=parset_id)

        # get manual parameters
        parset_instance = parset_entry.hydrate()
        mflists = {'keys': [], 'subkeys': [], 'types': [], 'values': [], 'labels': []}
        for param in parameters:
            mflists['keys'].append(param['key'])
            mflists['subkeys'].append(param['subkey'])
            mflists['types'].append(param['type'])
            mflists['labels'].append(param['label'])
            mflists['values'].append(param['value'])
        parset_instance.update(mflists)

        parset_entry.pars = op.saves(parset_instance.pars)
        parset_entry.updated = datetime.now(dateutil.tz.tzutc())
        db.session.add(parset_entry)
        ResultsDb.query.filter_by(parset_id=parset_id, 
            project_id=project_id, calculation_type=ResultsDb.CALIBRATION_TYPE).delete()
        result_entry = ResultsDb.query.filter_by(id=args['result_id']).first()
        result_entry.parset_id = parset_id
        result_entry.project_id = project_id
        result_entry.calculation_type = ResultsDb.CALIBRATION_TYPE
        db.session.add(result_entry)
        db.session.commit()

        return {
            "parset_id": parset_id,
            "parameters": args['parameters'],
            "result_id": result_entry.id
        }


manual_calibration_parser = RequestParser()
manual_calibration_parser.add_argument('maxtime', required=False, type=int, default=60)


class ParsetsAutomaticCalibration(Resource):

    @swagger.operation(
        summary='Launch auto calibration for the selected parset',
        parameters=manual_calibration_parser.swagger_parameters()
    )
    @report_exception
    def post(self, project_id, parset_id):
        from server.webapp.tasks import run_autofit, start_or_report_calculation
        from server.webapp.dbmodels import ParsetsDb

        args = manual_calibration_parser.parse_args()
        parset_entry = ParsetsDb.query.get(parset_id)
        parset_name = parset_entry.name

        can_start, can_join, wp_parset_id, work_type = start_or_report_calculation(project_id, parset_id, 'autofit')

        result = {'can_start': can_start, 'can_join': can_join, 'parset_id': wp_parset_id, 'work_type': work_type}
        if not can_start or not can_join:
            result['status'] = 'running'
            return result, 208
        else:
            run_autofit.delay(project_id, parset_name, args['maxtime'])
            result['status'] = 'started'
            result['maxtime'] = args['maxtime']
            return result, 201

    @report_exception
    def get(self, project_id, parset_id):
        from server.webapp.tasks import check_calculation_status
        from server.webapp.dbmodels import ParsetsDb

        parset_entry = ParsetsDb.query.get(parset_id)
        project_id = parset_entry.project_id

        status, error_text, start_time, stop_time, result_id = check_calculation_status(project_id, parset_id, 'autofit')
        return {'status': status, 'error_text': error_text, 'start_time': start_time, 'stop_time': stop_time, 'result_id': result_id}


file_upload_form_parser = RequestParser()
file_upload_form_parser.add_argument('file', type=AllowedSafeFilenameStorage, location='files', required=True)


class ParsetsData(Resource):
    """
    Export and import of the existing parset in / from pickled format.
    """
    method_decorators = [report_exception, login_required]

    @swagger.operation(
        produces='application/x-gzip',
        summary='Download data for the parset with the given id from project with the given id',
        notes="""
        if parset exists, returns data for it
        if parset does not exist, returns an error.
        """
    )
    @report_exception
    def get(self, project_id, parset_id):
        current_app.logger.debug("GET /api/project/{0}/parset/{1}/data".format(project_id, parset_id))
        parset_entry = db.session.query(ParsetsDb).filter_by(id=parset_id, project_id=project_id).first()
        if parset_entry is None:
            raise ParsetDoesNotExist(id=parset_id, project_id=project_id)

        # return result as a file
        loaddir = upload_dir_user(TEMPLATEDIR)
        if not loaddir:
            loaddir = TEMPLATEDIR

        filename = parset_entry.as_file(loaddir)

        response = helpers.send_from_directory(loaddir, filename)
        response.headers["Content-Disposition"] = "attachment; filename={}".format(filename)

        return response

    @swagger.operation(
        summary='Upload data for the parset with the given id in project with the given id',
        notes="""
        if parset exists, updates it with data from the file
        if parset does not exist, returns an error"""
    )
    @report_exception
    @marshal_with(ParsetsDb.resource_fields, envelope='parsets')
    def post(self, project_id, parset_id):
        # TODO replace this with app.config
        current_app.logger.debug("POST /api/project/{0}/parset/{1}/data".format(project_id, parset_id))

        print request.files, request.args
        args = file_upload_form_parser.parse_args()
        uploaded_file = args['file']

        project_entry = load_project_record(project_id, raise_exception=True)

        parset_entry = project_entry.find_parset(parset_id)
        parset_instance = op.loadobj(uploaded_file)

        parset_entry.restore(parset_instance)
        db.session.add(parset_entry)
        db.session.flush()

        # recalculate data (TODO: verify with Robyn if it's needed )
        project_instance = project_entry.hydrate()
        result = project_instance.runsim(parset_entry.name)
        current_app.logger.info("runsim result for project %s: %s" % (project_id, result))

        db.session.add(project_entry)  # todo: do we need to log that project was updated?
        db.session.flush()

        result_record = save_result(project_entry.id, result, parset_entry.name)
        db.session.add(result_record)

        db.session.commit()

        return [item.hydrate() for item in project_entry.parsets]

import uuid
import os
from datetime import datetime
import pprint
import json

import dateutil
from flask import current_app, helpers, request
from flask.ext.login import login_required
from flask_restful import Resource, marshal_with, fields
from flask_restful_swagger import swagger

import optima as op
from server.webapp.dataio import (
    load_project_record, TEMPLATEDIR, upload_dir_user, update_or_create_result_record,
    save_result, delete_result, save_parset, load_result, load_parset, load_project, load_result_record,
    load_parset_record, load_result_by_id, load_parset_list, get_parset_from_project)
from server.webapp.dbconn import db
from server.webapp.dbmodels import ParsetsDb, ResultsDb, ScenariosDb,OptimizationsDb
from server.webapp.exceptions import ParsetDoesNotExist, ParsetAlreadyExists
from server.webapp.parse import get_parameters_from_parset, put_parameters_in_parset
from server.webapp.resources.common import report_exception
from server.webapp.utils import AllowedSafeFilenameStorage, RequestParser, normalize_obj
from server.webapp.plot import make_mpld3_graph_dict



copy_parser = RequestParser()
copy_parser.add_arguments({
    'name': {'required': True},
    'parset_id': {'type': uuid.UUID}
})


class Parsets(Resource):
    """
    GET /api/project/<project_id>/parsets

    Returns all parsets of a project for display in dropdown menu in calibration

    POST /api/project/<project_id>/parsets

    Copy or make new parset that is passed in  body as parset_id
    """

    method_decorators = [report_exception, login_required]

    @swagger.operation(description='Download all parsets for project')
    @marshal_with(ParsetsDb.resource_fields, envelope='parsets')
    def get(self, project_id):
        current_app.logger.debug("/api/project/%s/parsets" % str(project_id))
        return load_parset_list(project_id)

    @swagger.operation(description='Create parset or copy existing parset')
    def post(self, project_id):
        current_app.logger.debug("POST /api/project/{}/parsets".format(project_id))
        args = copy_parser.parse_args()
        print "args", args
        name = args['name']
        parset_id = args.get('parset_id')

        project_record = load_project_record(project_id)
        project = project_record.hydrate()

        if name in project.parsets:
            raise ParsetAlreadyExists(project_id, name)

        if not parset_id:
            # CREATE parset with default settings
            project.makeparset(name, overwrite=False)
            result = project.runsim(name)
            project_record.restore(project)
            db.session.add(project_record)
            result_record = update_or_create_result_record(project_id, result, name)
            db.session.add(result_record)
        else:
            # COPY parset of parset_id to name
            original_parsets = [item for item in project_record.parsets if item.id == parset_id]
            if not original_parsets:
                raise ParsetDoesNotExist(parset_id, project_id=project_id)
            original_parset = original_parsets[0]
            parset_name = original_parset.name
            project.copyparset(orig=parset_name, new=name)
            project_record.restore(project)
            db.session.add(project_record)

        db.session.commit()

        rv = []
        for item in project_record.parsets:
            rv_item = item.hydrate().__dict__
            rv_item['id'] = item.id
            rv.append(rv_item)

        return rv


rename_parser = RequestParser()
rename_parser.add_argument('name', required=True)

class ParsetRenameDelete(Resource):
    """
    DELETE /api/project/<uuid:project_id>/parsets/<uuid:parset_id>

    Deletes given parset

    PUT /api/project/<uuid:project_id>/parsets/<uuid:parset_id>

    Renames the given parset
    """
    method_decorators = [report_exception, login_required]

    @swagger.operation(description='Delete parset with parset_id.')
    @report_exception
    def delete(self, project_id, parset_id):

        current_app.logger.debug("DELETE /api/project/{}/parsets/{}".format(project_id, parset_id))
        project_entry = load_project_record(project_id, raise_exception=True)

        parset = db.session.query(ParsetsDb).filter_by(project_id=project_entry.id, id=parset_id).first()
        if parset is None:
            raise ParsetDoesNotExist(id=parset_id, project_id=project_id)

        # TODO: also delete the corresponding calibration results
        db.session.query(ResultsDb).filter_by(
            project_id=project_id, id=parset_id, calculation_type=ResultsDb.DEFAULT_CALCULATION_TYPE).delete()
        db.session.query(ScenariosDb).filter_by(project_id=project_id,
            parset_id=parset_id).delete()
        db.session.query(OptimizationsDb).filter_by(project_id=project_id,
            parset_id=parset_id).delete()
        db.session.query(ParsetsDb).filter_by(project_id=project_id, id=parset_id).delete()
        db.session.commit()

        return '', 204

    @swagger.operation(description='Rename parset with parset_id')
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

        project_record = load_project_record(project_id, raise_exception=True)
        parset_records = [record for record in project_record.parsets if record.id == parset_id]
        if not parset_records:
            raise ParsetDoesNotExist(id=parset_id, project_id=project_id)
        parset_record = parset_records[0]
        parset_record.name = name
        db.session.add(parset_record)
        db.session.commit()
        return [record.hydrate() for record in project_record.parsets]



class ParsetCalibration(Resource):
    """
    /api/project/<uuid:project_id>/parsets/<uuid:parset_id>/calibration

    - GET: Returns parameter summaries and graphs for a project/parset, called on page init
    so doesn't really require a which
    - POST: Sends parameters to get graphs, with optional save, and optional autofit
    """

    method_decorators = [report_exception, login_required]

    @swagger.operation(description='Returns parameter summaries and graphs for a project/parset')
    def get(self, project_id, parset_id):
        """
        Returns the graphs for the parameter sets. If the results of the parameter
        sets exists, will use those results to generate the results, otherwise
        will the model, save the results, and then generate the graphs.

        Args:
            project_id: uid for project
            parset_id: uid for parset

        Url-query:
            autofit: boolean - true loads the results from the autofit parameters

        Returns:
             {
                 "calibration": {
                     "parset_id": uid_string,
                     "parameters": parameter_dictionary,
                     "graphs": mpl3_dict dictionary,
                     "resultId": uid_string
                 }
             }

        """
        current_app.logger.debug("/api/project/{}/parsets/{}/calibration".format(project_id, parset_id))
        autofit = request.args.get('autofit', False)
        calculation_type = 'autofit' if autofit else "calibration"

        print "> Calculation type: %s, autofit: %s" % (calculation_type, autofit)

        parset = load_parset(project_id, parset_id)
        parameters = get_parameters_from_parset(parset)

        result = load_result(project_id, parset_id, calculation_type)
        if result is not None:
            print "> Fetch result(%s) '%s' for parset '%s'" % (calculation_type, result.name, parset.name)
        else:
            print "> Runsim for new calibration results and store"
            project = load_project(project_id, autofit=autofit)
            result = project.runsim(simpars=parset.interp())
            save_result(project_id, result, parset.name, calculation_type)

        print "> Generating graphs"
        payload = {
            "calibration": {
                "parset_id": parset_id,
                "parameters": parameters,
                "resultId": result.uid,
            }
        }
        graph_dict = make_mpld3_graph_dict(result)
        payload["calibration"].update(graph_dict)
        return payload

    @swagger.operation(description='Updates parameters and returns graphs')
    def post(self, project_id, parset_id):
        """
        Generates selected graphs for an updated parameter set with optional save.
        Or returns graphs generated from an autofit pre-calculated result.

        Args:
            project_id: uid for project
            parset_id: uid for parset
        Post-body:
            parameters: list of model parameters and their values
            which: list of graphs to generate
            autofit: boolean indicates to fetch the autofit version of the results
        Returns:
             {
                 "calibration": {
                     "parset_id": uid_string,
                     "parameters": parameter_dictionary,
                     "graphs": mpl3_dict dictionary,
                     "resultId": uid_string
                 }
             }
        """
        args = normalize_obj(json.loads(request.data))
        autofit = args.get('autofit', False)
        which = args.get('which')
        if which is not None:
            which = map(str, which)

        if autofit:
            print "> Generating graphs from pre-autofitted results"
            # this is called after an autofit has run, which
            # will have deleted "calibration" results associated with the
            # optimized parset, written the optimized values to
            # the parset and saved the results to "autofit"
            result = load_result(project_id, parset_id, "autofit")
            if result is None:
                raise ValueError("Autofit results not found")
            parset = load_parset(project_id, parset_id)
            assert str(parset.uid) == str(parset_id)
            parameters = get_parameters_from_parset(parset)
            if which and 'improvement' not in which:
                which.insert(0, 'improvement')

        else:
            print '> Uploaded updated parameters for parset', parset_id
            project = load_project(project_id, raise_exception=True)
            parset = get_parset_from_project(project, parset_id)
            parameters = args.get('parameters')
            assert str(parset.uid) == str(parset_id)
            put_parameters_in_parset(parameters, parset)
            save_parset(project_id, parset)
            delete_result(project_id, parset.uid, "calibration")
            delete_result(project_id, parset.uid, "autofit")

            print "> Simulating model from uploaded parameters"
            result = project.runsim(simpars=parset.interp())
            save_result(project_id, result, parset.name, 'calibration')

        print "> Generate graphs"
        graphs = make_mpld3_graph_dict(result, which)

        payload = {
            'calibration': {
                "parset_id": parset_id,
                "parameters": parameters,
                "resultId": str(result.uid),
                "graphs": graphs["graphs"]
            }
        }
        return payload


class ParsetAutofit(Resource):
    """
    /api/project/<uuid:project_id>/parsets/<uuid:parset_id>/automatic_calibration

    - POST: Starts celery task to autofit parameters to historical data
    - GET: Returns the status for the current job:

    Returns:
        {
            'status': work_log.status,
            'error_text': work_log.error,
            'start_time': work_log.start_time,
            'stop_time': work_log.stop_time,
            'result_id': work_log.result_id
        }
    """

    method_decorators = [report_exception, login_required]

    @swagger.operation(summary='Launch auto calibration')
    def post(self, project_id, parset_id):
        """
        Launch auto calibration and returns the status code for the async job

        Args:
            project_id:
            parset_id:
        Post-query:
            maxtime: int - number of seconds to run

        """
        from server.webapp.tasks import run_autofit, start_or_report_calculation

        maxtime = json.loads(request.data).get('maxtime')
        calc_status = start_or_report_calculation(project_id, parset_id, 'autofit')
        if calc_status['status'] != "blocked":
            parset_name = load_parset_record(project_id, parset_id).name
            print "> Starting autofit for %s s" % maxtime
            run_autofit.delay(project_id, parset_name, maxtime)
            calc_status['status'] = 'started'
            calc_status['maxtime'] = maxtime
        return calc_status

    @swagger.operation(summary='Poll autofit status')
    def get(self, project_id, parset_id):
        """
        Returns status of auto calibration

        Args:
            project_id:
            parset_id:
        """
        from server.webapp.tasks import check_calculation_status
        print "> Checking calc state"
        calc_state = check_calculation_status(project_id, parset_id, 'autofit')
        pprint.pprint(calc_state, indent=2)
        if calc_state['status'] == 'error':
            raise Exception(calc_state['error_text'])
        return calc_state

file_upload_form_parser = RequestParser()
file_upload_form_parser.add_argument('file', type=AllowedSafeFilenameStorage, location='files', required=True)


class ParsetUploadDownload(Resource):
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

        result_record = update_or_create_result_record(project_entry.id, result, parset_entry.name)
        db.session.add(result_record)

        db.session.commit()

        return [item.hydrate() for item in project_entry.parsets]



class ResultsExportAsCsv(Resource):
    """
    /api/results/<results_id>
    - GET: returns a .csv file as blob
    """

    method_decorators = [report_exception, login_required]

    def get(self, result_id):
        """
        Export of data from an Optima Results object as a downloadable .csv file

        Args:
            result_id:
        """
        result = load_result_by_id(result_id)
        if result is None:
            raise Exception("Results '%s' does not exist" % result_id)

        load_dir = upload_dir_user(TEMPLATEDIR)
        if not load_dir:
            load_dir = TEMPLATEDIR
        filestem = 'results'
        filename = filestem + '.csv'

        result.export(filestem=os.path.join(load_dir, filestem))
        response = helpers.send_from_directory(load_dir, filename)
        response.headers["Content-Disposition"] = "attachment; filename={}".format(filename)
        return response


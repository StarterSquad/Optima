import pprint
import uuid
from pprint import pprint

import mpld3
from flask import current_app, request, helpers
from flask.ext.login import login_required
from flask_restful import Resource, marshal_with, marshal, fields
from flask_restful_swagger import swagger

from server.webapp.dataio import (
    load_project_record, load_progset_record, load_program, load_parset,
    update_or_create_program_record, get_target_popsizes, load_parameters_from_progset_parset)
from server.webapp.dbconn import db
from server.webapp.dbmodels import ProgsetsDb, ProgramsDb
from server.webapp.exceptions import (
    ProjectDoesNotExist, ProgsetDoesNotExist)
from server.webapp.resources.common import file_resource, file_upload_form_parser, report_exception
from server.webapp.utils import SubParser, Json, RequestParser, TEMPLATEDIR, upload_dir_user, normalize_obj


progset_parser = RequestParser()
progset_parser.add_arguments(
    {'name': {'required': True}, 'programs': {'type': Json, 'location': 'json'}})

class Progsets(Resource):
    """
    GET /api/project/<uuid:project_id>/progsets

    Download progsets for list in program-set manage page

    POST /api/project/<uuid:project_id>/progsets

    Save new project
    """
    method_decorators = [report_exception, login_required]

    @swagger.operation(description='Download progsets for the project with the given id.')
    @marshal_with(ProgsetsDb.resource_fields, envelope='progsets')
    def get(self, project_id):

        current_app.logger.debug("/api/project/%s/progsets" % project_id)
        project_record = load_project_record(project_id)
        if project_record is None:
            raise ProjectDoesNotExist(id=project_id)

        progsets_record = db.session.query(ProgsetsDb).filter_by(project_id=project_record.id).all()
        for progset_record in progsets_record:
            progset_record.get_extra_data()
            for program_record in progset_record.programs:
                program_record.get_optimizable()

        return progsets_record

    @swagger.operation(description='Create a progset for the project with the given id.')
    @marshal_with(ProgsetsDb.resource_fields)
    def post(self, project_id):
        current_app.logger.debug("/api/project/%s/progsets" % project_id)
        project_record = load_project_record(project_id)
        if project_record is None:
            raise ProjectDoesNotExist(id=project_id)
        args = progset_parser.parse_args()
        progset_record = ProgsetsDb(project_id, args['name'])
        progset_record.update_from_program_summaries(args['programs'], progset_record.id)
        progset_record.get_extra_data()
        db.session.add(progset_record)
        db.session.flush()
        db.session.commit()
        return progset_record, 201


class Progset(Resource):
    """
    GET /api/project/<uuid:project_id>/progsets/<uuid:progset_id>

    Download progset - is this ever used?

    PUT /api/project/<uuid:project_id>/progsets/<uuid:progset_id>

    Update existing project
    """
    method_decorators = [report_exception, login_required]

    @swagger.operation(description='Download progset with the given id.')
    @marshal_with(ProgsetsDb.resource_fields)
    def get(self, project_id, progset_id):
        current_app.logger.debug("/api/project/%s/progsets/%s" % (project_id, progset_id))
        progset_entry = load_progset_record(project_id, progset_id)
        progset_entry.get_extra_data()
        return progset_entry

    @swagger.operation(description='Update progset with the given id.')
    @marshal_with(ProgsetsDb.resource_fields)
    def put(self, project_id, progset_id):
        current_app.logger.debug("/api/project/%s/progsets/%s" % (project_id, progset_id))
        args = progset_parser.parse_args()

        progset_record = load_progset_record(project_id, progset_id)
        progset_record.name = args['name']

        program_summaries = normalize_obj(args.get('programs', []))
        progset_record.update_from_program_summaries(program_summaries, progset_id)
        progset_record.get_extra_data()

        db.session.commit()
        return progset_record

    @swagger.operation(description='Delete progset with the given id.')
    def delete(self, project_id, progset_id):
        current_app.logger.debug("/api/project/%s/progsets/%s" % (project_id, progset_id))
        progset_entry = db.session.query(ProgsetsDb).get(progset_id)
        if progset_entry is None:
            raise ProgsetDoesNotExist(id=progset_id)

        if progset_entry.project_id != project_id:
            raise ProgsetDoesNotExist(id=progset_id)

        db.session.query(ProgramsDb).filter_by(progset_id=progset_entry.id).delete()
        db.session.delete(progset_entry)
        db.session.commit()
        return '', 204


class ProgsetData(Resource):

    method_decorators = [report_exception, login_required]

    @swagger.operation(
        produces='application/x-gzip',
        description='Download progset with the given id as Binary.',
        notes="""
            if progset exists, returns it
            if progset does not exist, returns an error.
        """,

    )
    def get(self, project_id, progset_id):
        current_app.logger.debug("GET /api/project/{}/progsets/{}/data".format(project_id, progset_id))
        progset_entry = load_progset_record(project_id, progset_id)

        loaddir = upload_dir_user(TEMPLATEDIR)
        if not loaddir:
            loaddir = TEMPLATEDIR

        filename = progset_entry.as_file(loaddir)

        return helpers.send_from_directory(loaddir, filename)

    @swagger.operation(
        summary='Uploads data for already created progset',
        parameters=file_upload_form_parser.swagger_parameters()
    )
    @marshal_with(file_resource)
    def post(self, project_id, progset_id):
        """
        Uploads Data file, uses it to update the progrset and program models.
        Precondition: model should exist.
        """
        from server.webapp.parse import get_default_program_summaries

        current_app.logger.debug("POST /api/project/{}/progsets/{}/data".format(project_id, progset_id))

        args = file_upload_form_parser.parse_args()
        uploaded_file = args['file']

        source_filename = uploaded_file.source_filename

        progset_entry = load_progset_record(project_id, progset_id)

        project_entry = load_project_record(project_id)
        project = project_entry.hydrate()
        if project.data != {}:
            program_list = get_default_program_summaries(project)
        else:
            program_list = []

        from optima.utils import loadobj
        new_progset = loadobj(uploaded_file)
        progset_entry.restore(new_progset, program_list)
        db.session.add(progset_entry)

        db.session.commit()

        reply = {
            'file': source_filename,
            'result': 'Progset %s is updated' % progset_entry.name,
        }
        return reply


class ProgsetParameters(Resource):

    """
    GET /api/project/<uuid:project_id>/progsets/<uuid:progset_id>/parameters/<uuid:parset_id>

    Fetches parameters for a progset/parset combo to be used in outcome functions.
    """
    @swagger.operation(description='Get parameters sets for the selected progset')
    def get(self, project_id, progset_id, parset_id):
        return load_parameters_from_progset_parset(project_id, progset_id, parset_id)



class ProgsetEffects(Resource):
    """
    GET /api/project/<uuid:project_id>/progsets/<uuid:progset_id>/effects

    Fetch the effects of a given progset, given in the marshalled fields of
    a ProgsetsDB record, used in cost-coverage-ctrl.js

    PUT /api/project/<uuid:project_id>/progsets/<uuid:progset_id>/effects

    Saves the effects of a given progset, used in cost-coverage-ctrl.js
    """

    method_decorators = [report_exception, login_required]

    @swagger.operation(summary='Get List of existing Progset effects for the selected progset')
    def get(self, project_id, progset_id):
        from server.webapp.dataio import load_progset_record
        progset_record = load_progset_record(project_id, progset_id)
        return { 'effects': progset_record.effects }

    @swagger.operation(summary='Saves a list of outcomes')
    def put(self, project_id, progset_id):
        effects = request.get_json(force=True)
        from server.webapp.dataio import load_progset_record
        progset_record = load_progset_record(project_id, progset_id)
        progset_record.effects = normalize_obj(effects)
        db.session.add(progset_record)
        db.session.commit()
        return { 'effects': progset_record.effects }



query_program_parser = RequestParser()
query_program_parser.add_arguments({
    'program': {'required': True, 'type': Json, 'location': 'json'},
})

class Program(Resource):
    """
    POST /api/project/<project_id>/progsets/<progset_id>/program"

    Write program to web-server (for cost-coverage and outcome changes)
    The payload is JSON in the form:

    'program': {
      'active': True,
      'category': 'Care and treatment',
      'ccopars': { 'saturation': [[0.9, 0.9]],
                   't': [2016],
                   'unitcost': [[1.136849845773715, 1.136849845773715]]},
      'costcov': [ { 'cost': 16616289, 'coverage': 8173260, 'year': 2012},
                   { 'cost': 234234, 'coverage': 324234, 'year': 2013}],
      'created': 'Mon, 02 May 2016 05:27:48 -0000',
      'criteria': { 'hivstatus': 'allstates', 'pregnant': False},
      'id': '9b5db736-1026-11e6-8ffc-f36c0fc28d89',
      'name': 'HIV testing and counseling',
      'optimizable': True,
      'populations': [ 'FSW',
                       'Clients',
                       'Male Children 0-14',
                       'Female Children 0-14',
                       'Males 15-49',
                       'Females 15-49',
                       'Males 50+',
                       'Females 50+'],
      'progset_id': '9b55945c-1026-11e6-8ffc-130aba4858d2',
      'project_id': '9b118ef6-1026-11e6-8ffc-571b10a45a1c',
      'short': 'HTC',
      'targetpars': [ { 'active': True,
                        'param': 'hivtest',
                        'pops': [ 'FSW',
                                  'Clients',
                                  'Male Children 0-14',
                                  'Female Children 0-14',
                                  'Males 15-49',
                                  'Females 15-49',
                                  'Males 50+',
                                  'Females 50+']}],
      'updated': 'Mon, 02 May 2016 06:22:29 -0000'
    }
    """

    method_decorators = [report_exception, login_required]
    def post(self, project_id, progset_id):
        args = query_program_parser.parse_args()
        program_summary = normalize_obj(args['program'])
        program_record = update_or_create_program_record(
            project_id, progset_id, program_summary['short'],
            program_summary, program_summary['active'])
        current_app.logger.debug(
            "writing program = \n%s\n" % pprint.pformat(program_summary, indent=2))
        db.session.add(program_record)
        db.session.flush()
        db.session.commit()
        return 204


class ProgramPopSizes(Resource):
    """
    /api/project/{project_id}/progsets/{progset_id}/program/{program_id}/parset/{progset_id}/popsizes

    Return estimated popsize for a given program and parset. Used in
    cost-coverage function page to help estimate populations.
    """
    method_decorators = [report_exception, login_required]

    def get(self, project_id, progset_id, program_id, parset_id):
        payload = get_target_popsizes(project_id, parset_id, progset_id, program_id)
        return payload, 201



costcov_graph_parser = RequestParser()
costcov_graph_parser.add_arguments({
    't': {'required': True, 'type': str, 'location': 'args'},
    'parset_id': {'required': True, 'type': uuid.UUID, 'location': 'args'},
    'caption': {'type': str, 'location': 'args'},
    'xupperlim': {'type': long, 'location': 'args'},
    'perperson': {'type': bool, 'location': 'args'},
})

class ProgramCostcovGraph(Resource):
    """
    Costcoverage graph for a Program and a Parset (for population sizes).
    """

    method_decorators = [report_exception, login_required]

    def get(self, project_id, progset_id, program_id):
        """
        Args:
            t: comma-separated list of years (>= startyear in data)
            parset_id: parset ID of project (not related to program targetpars)
            caption: string to display in graph
            xupperlim: maximum dollar shown
            perperson: cost per person shown as data point

        Returns an mpld3 dict that can be displayed with the mpld3 plugin
        """
        args = costcov_graph_parser.parse_args()
        parset_id = args['parset_id']

        try:
            t = map(int, args['t'].split(','))
        except ValueError:
            raise ValueError("t must be a year or a comma-separated list of years.")

        plotoptions = {}
        for x in ['caption', 'xupperlim', 'perperson']:
            if args.get(x):
                plotoptions[x] = args[x]

        program = load_program(project_id, progset_id, program_id)
        parset = load_parset(project_id, parset_id)

        plot = program.plotcoverage(t=t, parset=parset, plotoptions=plotoptions)

        mpld3.plugins.connect(plot, mpld3.plugins.MousePosition(fontsize=14, fmt='.4r'))
        return normalize_obj(mpld3.fig_to_dict(plot))



import traceback

from flask import request, current_app, session, flash, redirect, url_for, abort

from flask_restful import Resource, marshal_with
from flask.ext.login import login_user, current_user, logout_user, login_required
from flask_restful_swagger import swagger

from server.webapp.dbconn import db
from server.webapp.dbmodels import UserDb

from server.webapp.utils import hashed_password, nullable_email, RequestParser
from server.webapp.exceptions import UserAlreadyExists, RecordDoesNotExist, InvalidCredentials
from server.webapp.resources.common import report_exception, verify_admin_request



user_parser = RequestParser()
user_parser.add_arguments({
    'email':       {'type': nullable_email, 'help': 'A valid e-mail address'},
    'displayName': {'dest': 'name'},
    'username':    {'required': True},
    'password':    {'type': hashed_password, 'required': True},
})


user_update_parser = RequestParser()
user_update_parser.add_arguments({
    'email':       {'type': nullable_email, 'help': 'A valid e-mail address'},
    'displayName': {'dest': 'name'},
    'username':    {'required': True},
    'password':    {'type': hashed_password},
})


class UserDoesNotExist(RecordDoesNotExist):
    _model = 'user'


class User(Resource):
    @swagger.operation(
        responseClass=UserDb.__name__,
        summary='List users'
    )
    @report_exception
    @marshal_with(UserDb.resource_fields, envelope='users')
    @verify_admin_request
    def get(self):
        current_app.logger.debug('/api/user/list {}'.format(request.args))
        return UserDb.query.all()

    @swagger.operation(
        responseClass=UserDb.__name__,
        summary='Create a user',
        parameters=user_parser.swagger_parameters()
    )
    @report_exception
    @marshal_with(UserDb.resource_fields)
    def post(self):
        current_app.logger.info("create request: {} {}".format(request, request.data))
        args = user_parser.parse_args()

        same_user_count = UserDb.query.filter_by(username=args.username).count()

        if same_user_count > 0:
            raise UserAlreadyExists(args.username)

        user = UserDb(**args)
        db.session.add(user)
        db.session.commit()

        return user, 201


class UserDetail(Resource):

    @swagger.operation(
        summary='Delete a user',
        notes='Requires admin privileges'
    )
    @report_exception
    @verify_admin_request
    def delete(self, user_id):
        current_app.logger.debug('/api/user/delete/{}'.format(user_id))
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

        current_app.logger.info("deleted user:{} {} {}".format(user_id, user_name, user_email))

        return '', 204

    @swagger.operation(
        responseClass=UserDb.__name__,
        summary='Update a user',
        notes='Requires admin privileges',
        parameters=user_update_parser.swagger_parameters(),
    )
    @report_exception
    @marshal_with(UserDb.resource_fields)
    def put(self, user_id):
        current_app.logger.debug('/api/user/{}'.format(user_id))

        user = UserDb.query.get(user_id)
        if user is None:
            raise UserDoesNotExist(user_id)

        try: userisanonymous = current_user.is_anonymous() # CK: WARNING, SUPER HACKY way of dealing with different Flask versions
    	except: userisanonymous = current_user.is_anonymous
        if userisanonymous or (str(user_id) != str(current_user.id) and not current_user.is_admin):
            secret = request.args.get('secret', '')
            u = UserDb.query.filter_by(password=secret, is_admin=True).first()
            if u is None:
                abort(403)

        args = user_update_parser.parse_args()
        for key, value in args.iteritems():
            if value is not None:
                setattr(user, key, value)

        db.session.commit()

        current_app.logger.info("modified user: {}".format(user_id))

        return user


# Authentication


user_login_parser = RequestParser()
user_login_parser.add_arguments({
    'username': {'required': True},
    'password': {'type': hashed_password, 'required': True},
})


class CurrentUser(Resource):
    method_decorators = [login_required]

    @swagger.operation(
        responseClass=UserDb.__name__,
        summary='Return the current user'
    )
    @report_exception
    @marshal_with(UserDb.resource_fields)
    def get(self):
        return current_user


class UserLogin(Resource):

    @swagger.operation(
        responseClass=UserDb.__name__,
        summary='Try to log a user in',
        parameters=user_login_parser.swagger_parameters()
    )
    @report_exception
    @marshal_with(UserDb.resource_fields)
    def post(self):
        current_app.logger.debug("/user/login {}".format(request.get_json(force=True)))

        try: userisanonymous = current_user.is_anonymous() # CK: WARNING, SUPER HACKY way of dealing with different Flask versions
    	except: userisanonymous = current_user.is_anonymous
        if userisanonymous:
            current_app.logger.debug("current user anonymous, proceed with logging in")

            args = user_login_parser.parse_args()
            try:
                # Get user for this username
                user = UserDb.query.filter_by(username=args['username']).first()

                # Make sure user is valid and password matches
                if user is not None and user.password == args['password']:
                    login_user(user)
                    return user

            except Exception:
                var = traceback.format_exc()
                print("Exception when logging user {}: \n{}".format(args['username'], var))

            raise InvalidCredentials

        else:
            return current_user


class UserLogout(Resource):
    method_decorators = [login_required]

    @swagger.operation(
        summary='Log the current user out'
    )
    @report_exception
    def get(self):
        current_app.logger.debug("logging out user {}".format(
            current_user.name
        ))
        logout_user()
        session.clear()
        flash(u'You have been signed out')

        return redirect(url_for("site"))

#!/bin/env python
# -*- coding: utf-8 -*-
"""
User Module
~~~~~~~~~~~~~~

1. Get current logged in user.
2. Login a user using username and password.
3. Logout.

"""
from flask import request, jsonify, g, session, flash, abort, Blueprint, url_for, current_app
from flask.ext.login import LoginManager, login_user, current_user, logout_user, redirect, login_required
from dbconn import db
from dbmodels import UserDb
from utils import verify_request
import logging

# route prefix: /api/user
user = Blueprint('user',  __name__, static_folder = '../static')

# System Imports
import hashlib

# Login Manager
login_manager = LoginManager()

@user.record
def record_params(setup_state):
    app = setup_state.app
    login_manager.init_app(app)

@user.before_request
def before_request():
    db.engine.dispose()
    g.user = None
    if 'user_id' in session:
        g.user = UserDb.query.filter_by(id=session['user_id']).first()

@user.route('/create', methods=['POST'])
def create_user():
    current_app.logger.info("create request: %s %s" % (request, request.data))
    current_app.logger.debug("/user/create %s" % request.get_json(force=True))
    # Check if the user already exists
    email = request.json['email']
    name = request.json['name']
#   password is now hashed on the client side using sha224
    password = request.json['password']

    if email is not None and name is not None and password is not None:
        # Get user for this username (if exists)
        try:
            no_of_users = UserDb.query.filter_by( email=email ).count()
        except:
            no_of_users = 0

        if no_of_users == 0:

            # Save to db
            u = UserDb(name, email, password)
            db.session.add( u )
            db.session.commit()

            # Login this user
            login_user(u)

            # Return user info
            return jsonify({'email': u.email, 'name': u.name })

    # We are here implies username is already taken
    return jsonify({'status': 'This email is already in use'})

@user.route('/login', methods=['POST'])
def login():
    current_app.logger.debug("/user/login %s" % request.get_json(force=True))
    # Make sure user is not logged in already.
    cu = current_user

    if cu.is_anonymous():
        current_app.logger.debug("current user anonymous")
        # Make sure user is valid.
        username = request.json['email']

        if username is not None:

            # Get hashsed password
            password = request.json['password']

            # Get user for this username
            try:
                u = UserDb.query.filter_by( email=username ).first()
            except:
                u = None

            # Make sure user is valid and password matches
            if u is not None and u.password == password:
                print("password:%s" % u.password)

                # Login the user
                login_user(u)

                # Return user info
                return jsonify({'email': u.email, 'name': u.name })

        # If we come here, login is not successful
        abort(401)
    else:
        current_app.logger.warning("User already logged in:%s %s" % (cu.name, cu.password))

    # User logged in, redirect
    return redirect(request.args.get("next") or url_for("site"))

@user.route('/current', methods=['GET'])
def current_user_api():
    cu = current_user
    if not cu.is_anonymous():
        return jsonify({ 'email': cu.email, 'name': cu.name })
    abort(401)

@user.route('/logout')
@login_required
def logout():
    cu = current_user
    username = cu.name
    current_app.logger.debug("logging out user %s" % username)
    logout_user()
    session.clear()
    flash(u'You have been signed out')
    current_app.logger.debug("User %s is signed out" % username)
    return redirect(url_for("site"))


#lists all the users. For internal reasons, this is implemented as console-only functionality
#with user hashed password as secret (can be changed later)
@user.route('/list')
@verify_request
def list():
    current_app.logger.debug('/api/user/list %s' % request.args)
    result = []
    users = UserDb.query.all()
    for u in users:
        result.append({'id':u.id, 'name':u.name, 'email':u.email})
    return jsonify({'users':result}) 

#deletes the given user by ID
@user.route('/delete/<user_id>', methods=['DELETE'])
@verify_request
def delete(user_id):
    current_app.logger.debug('/api/user/delete/%s' % user_id)
    user = UserDb.query.get(user_id)
    if not user:
        abort(404)
    else:
        user_email = user.email
        from dbmodels import ProjectDb, WorkingProjectDb
        from dbconn import db
        from sqlalchemy.orm import load_only
        #delete all corresponding projects and working projects as well
        projects = ProjectDb.query.filter_by(user_id=user_id).options(load_only("id")).all()
        project_ids = [project.id for project in projects]
        current_app.logger.debug("project_ids for user %s:%s" % (user_id, project_ids))
        WorkingProjectDb.query.filter(WorkingProjectDb.id.in_(project_ids)).delete(synchronize_session=False)
        ProjectDb.query.filter_by(user_id=user_id).delete()
        db.session.delete(user)
        db.session.commit()
        current_app.logger.info("deleted user:%s %s" % (user_id, user_email))
        return jsonify({'status':'OK','deleted':user_id})

#modify user by ID (can change email, name and/or password)
@user.route('/modify/<user_id>', methods=['PUT'])
@verify_request
def modify(user_id):
    current_app.logger.debug('/api/user/modify/%s' % user_id)
    user = UserDb.query.get(user_id)
    if not user:
        abort(404)
    else:
        from dbconn import db
        new_email = request.args.get('email')
        if new_email is not None: 
            user.email = new_email
        new_name = request.args.get('name')
        if new_name is not None:
            user.name = new_name
        new_password = request.args.get('password')
        #might change if we decide to hash PW twice
        if new_password is not None:
            user.password = new_password
        db.session.add(user)
        db.session.commit()
        current_app.logger.info("modified user:%s" % user_id)
        return jsonify({'status':'OK','modified':user_id})

#For Login Manager
@login_manager.user_loader
def load_user(userid):
    try:
        u = UserDb.query.filter_by(id=userid).first()
    except:
        u = None
    return u

@login_manager.unauthorized_handler
def unauthorized_handler():
    abort(401)

#!/bin/env python
# -*- coding: utf-8 -*-
from server.api import app, init_db
print "test package", __package__
from server.webapp.dbconn import db
import unittest
import hashlib
import json
from sqlalchemy.engine import reflection
from sqlalchemy.schema import (
        MetaData,
        Table,
        DropTable,
        ForeignKeyConstraint,
        DropConstraint,
    )



class OptimaTestCase(unittest.TestCase):
    """
    Baseclass for Optima API endpoint tests.

    This baseclass provides a setup that can be inherited to run API tests.

    Example:

        class SomeTestCase(OptimaTestCase):

            def setUp(self):
                super(SomeTestCase, self).setUp()

            def tearDown(self):
                super(SomeTestCase, self).tearDown()

            def test_a_response(self):
                response = self.test_client.get('/api/path')
                self.assertEqual(response.status_code, 200)


    """

    default_name = 'test'
    default_email = 'test@test.com'

    default_pops = [{"name": "Female sex workers", "short_name": "FSW", "sexworker": True, "injects": False, "sexmen": True, "client": False, "female": True, "male": False, "sexwomen": False}, \
        {"name": "Clients of sex workers", "short_name": "Clients", "sexworker": False, "injects": False, "sexmen": False, "client": True, "female": False, "male": True, "sexwomen": True}, \
        {"name": "Men who have sex with men", "short_name": "MSM", "sexworker": False, "injects": False, "sexmen": True, "client": False, "female": False, "male": True, "sexwomen": False}, \
        {"name": "Males who inject drugs", "short_name": "Male PWID", "sexworker": False, "injects": True, "sexmen": False, "client": False, "female": False, "male": True, "sexwomen": True}, \
        {"name": "Other males [enter age]", "short_name": "Other males", "sexworker": False, "injects": False, "sexmen": False, "client": False, "female": False, "male": True, "sexwomen": True}, \
        {"name": "Other females [enter age]", "short_name": "Other females", "sexworker": False, "injects": False, "sexmen": True, "client": False, "female": True, "male": False, "sexwomen": False}]

    def create_user(self, name = default_name, email = default_email):
        headers = {'Content-Type' : 'application/json'}
        create_data = '{"email":"%s","password":"%s","name":"%s"}' % (email, self.test_password, name)
        print ("create_user data: %s" % create_data)
        response = self.client.post('/api/user/create', data = create_data)
        return response

    def create_project(self, name):
        from server.webapp.dbmodels import ProjectDb
        """ Helper method to create project and save it to the database """
        project = ProjectDb(name, 1, '2000', '2010', OptimaTestCase.default_pops, '{}')
        db.session.add(project)
        db.session.flush()
        id = project.id
        db.session.commit()
        return id

    def api_create_project(self):
        project_data = json.dumps({
            'params': {
                'name': 'test',
                'datastart': 2000,
                'dataend': 2015,
                'populations': OptimaTestCase.default_pops
            },
        })
        headers = {'Content-Type': 'application/json'}
        response = self.client.post('/api/project/create/', data=project_data, headers=headers)
        return response

    def list_projects(self, user_id):
        from server.webapp.dbmodels import ProjectDb
        """ Helper method to list projects for the given user id"""
        projects = ProjectDb.query.filter_by(user_id=user_id).all()
        return [project for project in projects]

    def login(self, email=default_email, password=None):
        if not password: password = self.test_password
        headers = {'Content-Type' : 'application/json'}
        login_data = '{"email":"%s","password":"%s"}' % (email, password)
        self.client.post('/api/user/login', data=login_data, follow_redirects=True)

    def logout(self):
        self.client.get('/api/user/logout', follow_redirects=True)

    def setUp(self):
        self.test_password = hashlib.sha224("test").hexdigest()
        app.config['SQLALCHEMY_DATABASE_URI'] = 'postgresql+psycopg2://test:test@localhost:5432/optima_test'
        app.config['TESTING'] = True
        print "app created %s" % app
        init_db()
        print("db created. db: %s" % db)
        self.client = app.test_client()

    def _db_DropEverything(self, db):
        # From http://www.sqlalchemy.org/trac/wiki/UsageRecipes/DropEverything

        conn=db.engine.connect()

        # the transaction only applies if the DB supports
        # transactional DDL, i.e. Postgresql, MS SQL Server
        trans = conn.begin()

        inspector = reflection.Inspector.from_engine(db.engine)

        # gather all data first before dropping anything.
        # some DBs lock after things have been dropped in 
        # a transaction.
        metadata = MetaData()

        tbs = []
        all_fks = []

        for table_name in inspector.get_table_names():
            fks = []
            for fk in inspector.get_foreign_keys(table_name):
                if not fk['name']:
                    continue
                fks.append(
                    ForeignKeyConstraint((),(),name=fk['name'])
                    )
            t = Table(table_name,metadata,*fks)
            tbs.append(t)
            all_fks.extend(fks)

        for fkc in all_fks:
            conn.execute(DropConstraint(fkc))

        for table in tbs:
            conn.execute(DropTable(table))

        trans.commit()

    def tearDown(self):
        self.logout()
        db.session.remove()
        self._db_DropEverything(db)
        db.drop_all()
        db.get_engine(app).dispose()
        print "db dropped"

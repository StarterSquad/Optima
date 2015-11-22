from sqlalchemy.dialects.postgresql import JSON, UUID
from server.webapp.dbconn import db
from sqlalchemy import text
from sqlalchemy.orm import deferred

class UserDb(db.Model):
    __tablename__ = 'users'
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(60))
    email = db.Column(db.String(200))
    password = db.Column(db.String(200))
    is_admin = db.Column(db.Boolean, server_default=text('FALSE'))
    projects = db.relationship('ProjectDb', backref='users',
                                lazy='dynamic')

    def __init__(self, name, email, password, is_admin = False):
        self.name = name
        self.email = email
        self.password = password
        self.is_admin = is_admin

    def get_id(self):
        return self.id

    def is_active(self): # pylint: disable=R0201
        return True

    def is_anonymous(self): # pylint: disable=R0201
        return False

    def is_authenticated(self): # pylint: disable=R0201
        return True

class ProjectDb(db.Model):
    __tablename__ = 'projects'
    id = db.Column(UUID(True), server_default = text("uuid_generate_v1mc()"), primary_key = True)
    name = db.Column(db.String(60))
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'))
    datastart = db.Column(db.Integer)
    dataend = db.Column(db.Integer)
    populations = db.Column(JSON)
    created = db.Column(db.DateTime(timezone=True), server_default=text('now()'))
    updated = db.Column(db.DateTime(timezone=True), onupdate=db.func.now())
    version = db.Column(db.Text)
    settings = db.Column(db.LargeBinary)
    data = db.Column(db.LargeBinary)
    working_project = db.relationship('WorkingProjectDb', backref='projects',
                                uselist=False)
    project_data = db.relationship('ProjectDataDb', backref='projects',
                                uselist=False)
    parsets = db.relationship('ParsetsDb', backref = 'projects')
    results = db.relationship('ResultsDb', backref = 'results')

    def __init__(self, name, user_id, datastart, dataend, populations, version, 
        created = None, settings = None, data = None, parsets = None, results = None): # pylint: disable=R0913
        self.name = name
        self.user_id = user_id
        self.datastart = datastart
        self.dataend = dataend
        self.populations = populations
        self.created = created
        self.version = version
        self.settings = settings
        self.data = data
        self.parsets = parsets
        self.results = results

    def has_data(self):
        return self.data is not None

    def has_model_parameters(self):
        return self.parsets is not None

    def data_upload_time(self):
        return self.project_data.updated if self.project_data else None

class ParsetsDb(db.Model):
    __tablename__ = 'parsets'
    id = db.Column(UUID(True), server_default = text("uuid_generate_v1mc()"), primary_key = True)
    project_id = db.Column(UUID(True), db.ForeignKey('projects.id'))
    name = db.Column(db.Text)
    created = db.Column(db.DateTime(timezone=True), server_default=text('now()'))
    updated = db.Column(db.DateTime(timezone=True), onupdate=db.func.now())

    def __init__(self, project_id, name, created = None):
        self.project_id = project_id
        self.name = name
        self.created = created

class ResultsDb(db.Model):
    __tablename__ = 'results'
    id = db.Column(UUID(True), server_default = text("uuid_generate_v1mc()"), primary_key = True)
    parset_id = db.Column(UUID(True), db.ForeignKey('parsets.id'))
    project_id = db.Column(UUID(True), db.ForeignKey('projects.id'))
    calculation_type = db.Column(db.Text)
    pars = db.Column(db.LargeBinary)

    def __init__(self, parset_id, project_id, calculation_type, pars):
        self.parset_id = parset_id
        self.project_id = project_id
        self.calculation_type = calculation_type
        self.pars = pars


class WorkingProjectDb(db.Model): # pylint: disable=R0903
    __tablename__ = 'working_projects'
    id = db.Column(UUID(True),db.ForeignKey('projects.id'), primary_key=True )
    is_working = db.Column(db.Boolean, unique=False, default=False)
    work_type = db.Column(db.String(32), default=None)
    project = db.Column(db.LargeBinary)
    work_log_id = db.Column(UUID(True), default = None)

    def __init__(self, project_id, is_working=False, project = None, work_type = None, work_log_id = None): # pylint: disable=R0913
        self.id = project_id
        self.project = project
        self.is_working = is_working
        self.work_type = work_type
        self.work_log_id = work_log_id

class WorkLogDb(db.Model): # pylint: disable=R0903
    __tablename__ = "work_log"

    work_status = db.Enum('started', 'completed', 'cancelled', 'error' , name='work_status')

    id = db.Column(db.Integer, primary_key=True)
    work_type = db.Column(db.String(32), default = None)
    project_id = db.Column(db.Integer, db.ForeignKey('projects.id'), index = True)
    start_time = db.Column(db.DateTime(timezone=True), server_default=text('now()'))
    stop_time = db.Column(db.DateTime(timezone=True), default = None)
    status = db.Column(work_status, default='started')
    error = db.Column(db.Text, default = None)

    def __init__(self, project_id, work_type = None):
        self.project_id = project_id
        self.work_type = work_type

class ProjectDataDb(db.Model): # pylint: disable=R0903
    __tablename__ = 'project_data'
    id = db.Column(UUID(True),db.ForeignKey('projects.id'), primary_key=True )
    meta = deferred(db.Column(db.LargeBinary))
    upload_time = db.Column(db.DateTime(timezone=True), server_default=text('now()'))

    def __init__(self, project_id, meta, upload_time = None):
        self.id = project_id
        self.meta = meta
        self.upload_time = upload_time

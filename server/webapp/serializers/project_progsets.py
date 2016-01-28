from flask_restful import fields

from server.webapp.inputs import SubParser, Json as JsonInput
from server.webapp.fields import Json, Uuid
from server.webapp.utils import RequestParser


# Progset

costcov_parser = RequestParser()
costcov_parser.add_arguments({
    'year': {'required': True, 'location': 'json'},
    'spending': {'required': True, 'type': float, 'location': 'json', 'dest': 'cost'},
    'coverage': {'required': True, 'type': float, 'location': 'json', 'dest': 'cov'},
})

program_parser = RequestParser()
program_parser.add_arguments({
    'name': {'required': True, 'location': 'json'},
    'short': {'location': 'json'},
    'short_name': {'location': 'json'},
    'category': {'required': True, 'location': 'json'},
    'active': {'type': bool, 'default': False, 'location': 'json'},
    'parameters': {'type': list, 'dest': 'pars', 'location': 'json'},
    'populations': {'type': list, 'location': 'json', 'dest': 'targetpops'},
    'addData': {'type': SubParser(costcov_parser), 'dest': 'costcov', 'action': 'append', 'default': []},
})


progset_parser = RequestParser()
progset_parser.add_arguments({
    'name': {'required': True},
    'programs': {'required': True, 'type': SubParser(program_parser), 'action': 'append'}
})


param_fields = {
        'name': fields.String,
        'populations': Json,
        'coverage': fields.Boolean,
    }


program_effect_parser = RequestParser()
program_effect_parser.add_arguments({
    'name': {'required': False, 'location': 'json'},
    'intercept_lower': {'required': False, 'type': float, 'location': 'json'},
    'intercept_upper': {'required': False, 'type': float, 'location': 'json'},
})


param_year_effect_parser = RequestParser()
param_year_effect_parser.add_arguments({
    'year': {'required': False, 'location': 'json'},
    'intercept_lower': {'required': False, 'type': float, 'location': 'json'},
    'intercept_upper': {'required': False, 'type': float, 'location': 'json'},
    'interact': {'location': 'json', 'required': False, 'type': JsonInput},
    'programs': {
        'type': SubParser(program_effect_parser),
        # 'action': 'append',
        'default': [],
        'location': 'json',
        'required': False,
    },
})


param_effect_parser = RequestParser()
param_effect_parser.add_arguments({
    'name': {'required': False, 'location': 'json'},
    'pop': {'required': False, 'location': 'json', 'type': JsonInput},
    'years': {
        'type': SubParser(param_year_effect_parser),
        # 'action': 'append',
        'default': [],
        'location': 'json',
        'required': False,
    },
})


parset_effect_parser = RequestParser()
parset_effect_parser.add_arguments({
    'parset': {'required': False, 'location': 'json'},
    'parameters': {
        'type': SubParser(param_effect_parser),
        # 'action': 'append',
        'default': [],
        'location': 'json',
        'required': False
    }
})


effect_parser = RequestParser()
effect_parser.add_argument('effects', type=SubParser(parset_effect_parser), action='append')


program_effect_fields = {
    'name': fields.String,
    'intercept_lower': fields.Float,
    'intercept_upper': fields.Float
}


param_year_effect_fields = {
    'year': fields.Integer,
    'intercept_lower': fields.Float,
    'intercept_upper': fields.Float,
    'interact': fields.List(fields.String),
    'programs': fields.List(fields.Nested(program_effect_fields), default=[])
}

param_effect_fields = {
    'name': fields.String,
    'pop': fields.Raw,  # ToDo implement a field type that matches String or List of Strings
    'years': fields.List(fields.Nested(param_year_effect_fields), default=[])
}


parset_effect_fields = {
    'parset': Uuid,
    'parameters': fields.List(fields.Nested(param_effect_fields), default=[])
}


progset_effects_fields = {
    'effects': fields.List(fields.Nested(parset_effect_fields), default=[])
}

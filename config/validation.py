from cerberus import Validator, errors
import yaml

schema = {
    "region": {"type": "string", "allowed": ["Germany", "Europe"]},
    "forward" : {"type": "boolean"},
    "addtree": {"type": "boolean"},
    "backgrd" : {"type": "boolean"},
    "mcgrath" : {"type": "boolean"},
    "mcgback" : {"type": "boolean"},
    "irri" : {"type": "boolean"},
    "syear" : {"type": "integer", "min": 1900, "max": 2015},
    "eyear" : {"type": "integer", "min": 1900, "max": 2015},
    "npfts" : {"type": "integer"},
    "esayear" : {"type": "integer"},
    "mcgrath_eyear" : {"type": "integer"},
    "xsize" : {"type": "integer"},
    "ysize" : {"type": "integer"},
    "gradef" : {"type": "integer"},
    "crodef" : {"type": "integer"},
    "shrdef" : {"type": "integer"},
    "oname" : {"type": "string"},
    "lcd" : {"type": "string"},
    "mcg" : {"type": "string"},
    "vers" : {"type": "string"},
    "res" : {"type": "integer"},
    "remap" : {"type": "string"},
    "scenario" : {"type": "string"},
    "trans" : {"type": "boolean"},
    "state" : {"type": "boolean"},
}


class CustomValidator(Validator):
    def _validate_allowed(self, allowed, field, value):
        """{'type': 'list'}"""
        if value not in allowed:
            self._error(field, f"unallowed value {value} --> Select one of the following values: {allowed}")


def validate_config(config):
    v = CustomValidator(schema)
    v.validate(config)
    if v.errors:
        raise ValueError(v.errors)

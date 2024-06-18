from cerberus import Validator, errors
import yaml
import os
import xarray as xr

schema = {
    "region": {"type": "string", "allowed": ["Germany", "Europe"]},
    "forward" : {"type": "boolean"},
    "addtree": {"type": "boolean"},
    "backgrd" : {"type": "boolean"},
    "mcgrath" : {"type": "boolean"},
    "mcgback" : {"type": "boolean"},
    "irri" : {"type": "boolean"},
    "prepare_luh2_data": {"type": "boolean"},
    "prepare_mcgrath": {"type": "boolean"},
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
    "plot" : {"type": "boolean"},
    "plot_npft" : {"type": "integer", "nullable": True},
    "plot_year" : {"type": "integer", "nullable": True},
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

def validate_main_files(config):
    pass

def validate_prepared_files(namelist):
    for key, value in namelist.items():
        # checking if files exist
        if key.startswith("F_") and key not in ["F_IRRI_IN", "F_ADDTREE"]:
            if not os.path.isfile(value):
                raise ValueError(f"File {value} does not exist")
        # checking file sizes
        if key.startswith("F_") and key not in ["F_IRRI_IN", "F_ADDTREE"]:
            ds = xr.open_dataset(value, decode_times=False)
            if not (ds.dims.get("x") or ds.dims.get("lon")) and not (ds.dims.get("y") or ds.dims.get("lat")):
                raise ValueError(f"File {value} has wrong dimensions. Expected 2D but got {ds.dims}")
            else:
                x_dim = "x" if ds.dims.get("x") else "lon"
                y_dim = "y" if ds.dims.get("y") else "lat"
            if ds.dims.get(x_dim) != namelist["XSIZE"] or ds.dims.get(y_dim) != namelist["YSIZE"]:
                raise ValueError(f"File {value} has wrong dimensions. Expected {namelist['XSIZE']}x{namelist['YSIZE']} but got {ds.dims.get(x_dim)}x{ds.dims.get(y_dim)}")
        # checking irrigation file
        if (namelist.get("IRRI") and key=="F_IRRI_IN") or (namelist.get("ADDTREE") and key=="F_ADDTREE"):
            if not os.path.isfile(value):
                raise ValueError(f"File {value} does not exist")

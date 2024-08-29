from cerberus import Validator, errors
import yaml
import os
import xarray as xr
from src.config import datadir, scenario_dict

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
    "syear" : {"type": "integer"},
    "eyear" : {"type": "integer"},
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
    "remap" : {"type": "string", "allowed": ["bilinear", "con2"]},
    "scenario" : {"type": "string", "allowed": ["historical", "historical_high", "historical_low", "rcp19", "rcp26", "rcp34", "rcp45", "rcp60", "rcp70", "rcp85"]},
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
    if config.syear >= config.eyear:
        raise ValueError("Starting year (syear) must be smaller than ending year")

def validate_main_files(config):
    # IT MIGHT BE SIMPLIFIED AND INCLUDED JUST IN THE INIT FUNCTION OF LUT
    if config.scenario in ["historical", "historical_high", "historical_low"]:
        sfile="states"
        tfile="transitions"
        mfile="management"
    elif config.scenario in scenario_dict.keys():
        afile=f"added_tree_cover_input4MIPs_landState_ScenarioMIP_UofMD-{scenario_dict[self.scenario]}-2-1-f_gn_2015-2100"
        sfile=f"multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-{scenario_dict[self.scenario]}-2-1-f_gn_2015-2100"
        tfile=f"multiple-transitions_input4MIPs_landState_ScenarioMIP_UofMD-{scenario_dict[self.scenario]}-2-1-f_gn_2015-2100"
        mfile=f"multiple-management_input4MIPs_landState_ScenarioMIP_UofMD-{scenario_dict[self.scenario]}-2-1-f_gn_2015-2100"
    for value in [sfile, tfile, mfile]:
        if not os.path.isfile(os.path.join(datadir, value)+".nc"):
            raise ValueError(f"File {value} does not exist")   

def validate_prepared_files(namelist):
    for key, value in namelist.items():
        # checking if files exist
        if key == "F_GRID":
            if not os.path.isfile(value):
                raise ValueError(f"File {value} does not exist")
        elif key.startswith("F_") and key not in ["F_IRRI_IN", "F_ADDTREE", "F_MCGRATH", "F_LC_OUT"]:
            if not os.path.isfile(value):
                raise ValueError(f"File {value} does not exist")
        # checking file sizes
        if key.startswith("F_") and key not in ["F_IRRI_IN", "F_ADDTREE", "F_GRID", "F_LC_OUT"]:
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

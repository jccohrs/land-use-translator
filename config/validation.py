from cerberus import Validator, errors
import yaml
import os
import xarray as xr
from src.config import datadir, scenario_dict, mcg

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
    "mcgrath_eyear" : {"type": "integer", "nullable": True},
    "xsize" : {"type": "integer"},
    "ysize" : {"type": "integer"},
    "gradef" : {"type": "integer"},
    "crodef" : {"type": "integer"},
    "shrdef" : {"type": "integer"},
    "remap" : {"type": "string", "allowed": ["bilinear", "con2"]},
    "scenario" : {"type": "string", "allowed": ["historical", "historical_high", "historical_low", "rcp19", "rcp26", "rcp34", "rcp45", "rcp60", "rcp70", "rcp85"]},
    "grid": {"type": "float"},
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
        raise ValueError("Starting year (syear) must be smaller than ending year (eyear)")
    if config.mcgrath_eyear:
        if config.mcgrath_eyear > config.eyear:
            raise ValueError("Mcgrath year (mcgrath_eyear) must be equal or smaller than ending year (eyear)")
        if config.mcgrath_eyear < config.syear:
            raise ValueError("Mcgrath year (mcgrath_eyear) must be equal or bigger than starting year (syear)")
    if config.plot_year:
        if config.plot_year > (config.eyear-config.syear):
            raise ValueError("Plot year (plot_year) must be smaller than the number of years in the simulation")

def validation_outline(datadir, file):
    if not os.path.isfile(os.path.join(datadir, file)):
        path = os.path.join(datadir, file)
        raise ValueError(f"File {path} does not exist")

def validate_main_files(config):
    # IT MIGHT BE SIMPLIFIED AND INCLUDED JUST IN THE INIT FUNCTION OF LUT
    if config.scenario in ["historical", "historical_high", "historical_low"]:
        sfile="states"
        tfile="transitions"
        mfile="management"
    elif config.scenario in scenario_dict.keys():
        afile=f"added_tree_cover_input4MIPs_landState_ScenarioMIP_UofMD-{scenario_dict[config.scenario]}-2-1-f_gn_2015-2100"
        sfile=f"multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-{scenario_dict[config.scenario]}-2-1-f_gn_2015-2100"
        tfile=f"multiple-transitions_input4MIPs_landState_ScenarioMIP_UofMD-{scenario_dict[config.scenario]}-2-1-f_gn_2015-2100"
        mfile=f"multiple-management_input4MIPs_landState_ScenarioMIP_UofMD-{scenario_dict[config.scenario]}-2-1-f_gn_2015-2100"
    if config.mcgrath or config.prepare_mcgrath:
        ifile=f"{mcg}_{config.syear}_{config.mcgrath_eyear}.nc"
        validation_outline(datadir, ifile)
    if config.trans:
        validation_outline(datadir, tfile+".nc")
    if config.state:
        validation_outline(datadir, sfile+".nc")
    if config.irri:
        if config.scenario not in ["historical", "historical_low", "historical_high"]:
            validation_outline(datadir, sfile+".nc")
        validation_outline(datadir, mfile+".nc")
    if config.addtree:
        validation_outline(datadir, afile+".nc")

def validate_prepared_files(namelist, config):
    for key, value in namelist.items():
        # checking if files exist
        if key == "F_GRID":
            if not os.path.isfile(value):
                raise ValueError(f"File {value} does not exist")
        elif key.startswith("F_") and key not in ["F_IRRI_IN", "F_ADDTREE", "F_MCGRATH", "F_LC_OUT"]:
            if not os.path.isfile(value):
                raise ValueError(f"File {value} does not exist")
        # checking file sizes
        
        if key.startswith("F_") and key not in ["F_IRRI_IN", "F_ADDTREE", "F_MCGRATH", "F_GRID", "F_LC_OUT"]:
            ds = xr.open_dataset(value, decode_times=False)
            if not (ds.dims.get("x") or ds.dims.get("lon")) and not (ds.dims.get("y") or ds.dims.get("lat")):
                raise ValueError(f"File {value} has wrong dimensions. Expected 2D but got {ds.dims}")
            else:
                x_dim = "x" if ds.dims.get("x") else "lon"
                y_dim = "y" if ds.dims.get("y") else "lat"
            if ds.dims.get(x_dim) != config.xsize or ds.dims.get(y_dim) != config.ysize:
                raise ValueError(f"File {value} has wrong dimensions. Expected {config.xisze}x{config.ysize} but got {ds.dims.get(x_dim)}x{ds.dims.get(y_dim)}")
        # checking irrigation file
        if (namelist.get("IRRI") and key=="F_IRRI_IN") or (namelist.get("ADDTREE") and key=="F_ADDTREE") or (namelist.get("MCGRATH") and key=="F_MCGRATH"):
            if not os.path.isfile(value):
                raise ValueError(f"File {value} does not exist")


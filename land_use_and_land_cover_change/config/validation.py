from cerberus import Validator
import os
import xarray as xr
from cdo import Cdo
from lut_config import datadir, scriptsdir, scenario_dict, mcg, tf_file_syear, th_file_syear

cdo = Cdo()

schema = {
    "region": {"type": "string", "allowed": ["Germany", "Europe", "WestAfrica", "NorthAmerica", "Australasia"]},
    "forward" : {"type": "boolean"},
    "addtree": {"type": "boolean"},
    "backgrd" : {"type": "boolean"},
    "mcgrath" : {"type": "boolean"},
    "irri" : {"type": "boolean"},
    "prepare_luh2_data": {"type": "boolean"},
    "prepare_mcgrath": {"type": "boolean"},
    "syear" : {"type": "integer"},
    "eyear" : {"type": "integer"},
    "npfts" : {"type": "integer"},
    "mcgrath_eyear" : {"type": "integer", "nullable": True},
    "xsize" : {"type": "integer"},
    "ysize" : {"type": "integer"},
    "gradef" : {"type": "integer"},
    "crodef" : {"type": "integer"},
    "shrdef" : {"type": "integer"},
    "remap" : {"type": "string", "allowed": ["bilinear", "con2"]},
    "scenario" : {"type": "string", "allowed": ["historical", "historical_high", "historical_low", "rcp19", "rcp26", "rcp34", "rcp45", "rcp60", "rcp70", "rcp85"]},
    "grid": {"type": "float"},
    "rcm_lsm_var": {"type": "string"},
    "pfts_file_var": {"type": "string"},
    "path_file_states" : {"type": "string", "nullable": True},
    "path_file_trans" : {"type": "string", "nullable": True},
    "path_file_manag" : {"type": "string", "nullable": True},
    "path_file_addtree" : {"type": "string", "nullable": True},
    "path_file_lc_in" : {"type": "string", "nullable": True},
    "path_file_backgra_global" : {"type": "string", "nullable": True},
    "path_file_backshr_global" : {"type": "string", "nullable": True},
    "path_file_backfor_global" : {"type": "string", "nullable": True},
    "path_file_backurb_global" : {"type": "string", "nullable": True},
    "path_file_backcro_global" : {"type": "string", "nullable": True},
    "path_file_backgra" : {"type": "string", "nullable": True},
    "path_file_backshr" : {"type": "string", "nullable": True},
    "path_file_backfor" : {"type": "string", "nullable": True},
    "path_file_backurb" : {"type": "string", "nullable": True},
    "path_file_backcro" : {"type": "string", "nullable": True},
    "path_file_lsm" : {"type": "string", "nullable": True},
    "coords": {"type": "string", "nullable": True},
}


class CustomValidator(Validator):
    def _validate_allowed(self, allowed, field, value):
        """{'type': 'list'}"""
        if value not in allowed:
            self._error(field, f"unallowed value {value} --> Select one of the following values: {allowed}")

def validate_config(config) -> None:
    v = CustomValidator(schema)
    v.validate(config)
    if v.errors:
        raise ValueError(v.errors)
    if config.syear >= config.eyear:
        raise ValueError("Starting year (syear) must be smaller than ending year (eyear)")
    if config.mcgrath_eyear and config.mcgrath:
        if config.mcgrath_eyear > config.eyear:
            raise ValueError("Mcgrath year (mcgrath_eyear) must be equal or smaller than ending year (eyear)")
        if config.mcgrath_eyear < config.syear:
            raise ValueError("Mcgrath year (mcgrath_eyear) must be equal or bigger than starting year (syear)")
    if config.coords:
        if len(config.coords.split(",")) != 4:
            raise ValueError("Coordinates must given as 4 values (lonmin,lonmax,latmin,latmax) separated by commas")
        try:
            [float(config.coords.split(",")[i]) for i in range(4)]
        except ValueError:
            raise ValueError("Coordinates must be given as float values")

def validate_pfts_file(namelist, config):
    ds = xr.open_dataset(f"{namelist['F_LC_IN_REG'].replace('.nc','_tmp.nc')}")
    x_dim = "x" if ds.sizes.get("x") else "lon" if ds.sizes.get("lon") else "rlon"
    y_dim = "y" if ds.sizes.get("y") else "lat" if ds.sizes.get("lat") else "rlat"
    try: 
        year = config.syear if config.forward else config.eyear
        ds.sel(time=str(year))
    except KeyError:
        raise KeyError(f"Year {year} not found in file. Check if the file has the correct time dimensions or if syear/eyear are correctly set.")
    if config.xsize != ds.sizes.get(x_dim) or config.ysize != ds.sizes.get(y_dim):
        raise ValueError(f"Wrong sizes given. Expected {config.xsize}x{config.ysize} but got {ds.sizes.get(x_dim)}x{ds.sizes.get(y_dim)}. Check if the sizes coincide with the correspondig grid_reg file located in {scriptsdir}/.")

def validate_path(file, datadir=None, add_message=""):
    if datadir:
        if not os.path.isfile(os.path.join(datadir, file)):
            path = os.path.join(datadir, file)
            raise ValueError(f"File {path} does not exist. "+add_message)
    else:
        if not os.path.isfile(file):
            raise ValueError(f"File {file} does not exist. "+add_message)

def validate_main_files(namelist, config):
    # IT MIGHT BE SIMPLIFIED AND INCLUDED JUST IN THE INIT FUNCTION OF LUT
    for key, value in namelist.items():
        # checking if files exist
        if key == "F_GRID":
            validate_path(value)
        if key in ("F_LC_IN"):
            validate_path(value)
            validate_variable(value, config.pfts_file_var)
            validate_dimensions(value, config, type=2)
        if key.startswith("F_GLOBAL_BACK") and config.backgrd:
            validate_path(value)
            validate_dimensions(value, config, type=2)
        if config.path_file_lsm:
            validate_path(config.path_file_lsm)
            validate_dimensions(config.path_file_lsm, config)
            ds = xr.open_dataset(config.path_file_lsm, decode_times=False)
            try:
                ds[config.rcm_lsm_var]
            except KeyError:
                raise ValueError(f"Variable {config.rcm_lsm_var} not found in file {config.path_file_lsm}")
    if config.scenario in ["historical", "historical_high", "historical_low"]:
        sfile="states.nc" if not config.path_file_states else config.path_file_states
        tfile="transitions.nc" if not config.path_file_trans else config.path_file_trans
        mfile="management.nc" if not config.path_file_manag else config.path_file_manag
    elif config.scenario in scenario_dict.keys():
        afile=f"added_tree_cover_input4MIPs_landState_ScenarioMIP_UofMD-{scenario_dict[config.scenario]}-2-1-f_gn_2015-2100.nc" if not config.path_file_addtree else config.path_file_addtree
        sfile=f"multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-{scenario_dict[config.scenario]}-2-1-f_gn_2015-2100.nc" if not config.path_file_states else config.path_file_states
        tfile=f"multiple-transitions_input4MIPs_landState_ScenarioMIP_UofMD-{scenario_dict[config.scenario]}-2-1-f_gn_2015-2100.nc" if not config.path_file_trans else config.path_file_trans
        mfile=f"multiple-management_input4MIPs_landState_ScenarioMIP_UofMD-{scenario_dict[config.scenario]}-2-1-f_gn_2015-2100.nc" if not config.path_file_manag else config.path_file_manag
    if (config.prepare_mcgrath):
        ifile=f"{mcg}_{config.syear}_{config.mcgrath_eyear}.nc"
        validate_path(ifile, datadir, add_message="Consider adding the missing file or unable option 'prepare_mcgrath'.")
    if config.prepare_luh2_data:
        validate_path(tfile, datadir)
        input_file = f"{datadir}/{tfile}" if not config.path_file_trans else config.path_file_trans
        validate_timerange(input_file, config)
        validate_path(sfile, datadir)
    if config.irri:
        if config.scenario not in ["historical", "historical_low", "historical_high"]:
            validate_path(sfile, datadir)
        validate_path(mfile, datadir)
    if config.addtree:
        validate_path(afile, datadir)

def validate_mcgrath_prepared_files(namelist, config):
    for key, value in namelist.items():
        if (key=="F_MCGRATH"):
            validate_path(value, add_message="Consider enabling the option 'prepare_mcgrath' to prepare the missing file or unable 'mcgrath'.")
            validate_dimensions(value, config)

def validate_prepared_files(namelist, config):
    for key, value in namelist.items():
        # checking if files exist
        if key == "F_GRID":
            validate_path(value)
        elif key not in ["F_IRRI_IN", "F_ADDTREE", "F_MCGRATH", "F_LC_OUT", "F_LC_IN"] and not key.startswith("F_BACK") and not key.startswith("F_GLOBAL"):
            validate_path(value, add_message="Consider enabling the option 'prepare_luh2_data' to prepare the missing file.")
        # checking file sizes
        
        if key not in ["F_IRRI_IN", "F_ADDTREE", "F_MCGRATH", "F_GRID", "F_LC_OUT", "F_LC_IN"] and not key.startswith("F_BACK") and not key.startswith("F_GLOBAL"):
            validate_dimensions(value, config)

        if config.backgrd and key.startswith("F_BACK"):
            validate_dimensions(value, config)

        if (config.addtree and key=="F_ADDTREE"):
            validate_path(value)
            validate_dimensions(value, config)

        if (config.irri and key=="F_IRRI_IN"):
            validate_path(value, add_message="Consider enabling the option 'prepare_luh2_data' to prepare the missing file.")
            validate_dimensions(value, config)

def validate_dimensions(value, config, type=1):
    ds = xr.open_dataset(value, decode_times=False)
    if not (ds.sizes.get("x") or ds.sizes.get("lon") or ds.sizes.get("rlon")) and not (ds.sizes.get("y") or ds.sizes.get("lat") or ds.sizes.get("rlat")):
        raise ValueError(f"File {value} has wrong dimensions. Expected 2D but got {ds.sizes}")
    else:
        x_dim = "x" if ds.sizes.get("x") else "lon" if ds.sizes.get("lon") else "rlon"
        y_dim = "y" if ds.sizes.get("y") else "lat" if ds.sizes.get("lat") else "rlat"
    if type == 2:
        if ds.sizes.get(x_dim) < config.xsize or ds.sizes.get(y_dim) < config.ysize:
            raise ValueError(f"File {value} has wrong dimensions. Expected {config.xsize}x{config.ysize} but got {ds.sizes.get(x_dim)}x{ds.sizes.get(y_dim)}")
    else:
        if ds.sizes.get(x_dim) != config.xsize or ds.sizes.get(y_dim) != config.ysize:
            raise ValueError(f"File {value} has wrong dimensions. Expected {config.xsize}x{config.ysize} but got {ds.sizes.get(x_dim)}x{ds.sizes.get(y_dim)}")

def validate_variable(value, var):
    ds = xr.open_dataset(value, decode_times=False)
    if "var801" not in ds.variables and var not in ds.variables:
        raise ValueError(f"Neither Variable {var} or var8* found in file {value}")

def validate_timerange(input_file, config):
    init_year = tf_file_syear if config.forward else th_file_syear
    ds = xr.open_dataset(input_file, decode_times=False).sel(time=slice(config.syear-init_year, config.eyear-init_year))
    if not ds.time.units.split()[2].split("-")[0] == str(init_year):
        raise ValueError(f"File {input_file} has wrong time dimensions. Did not found time range {config.syear}/{config.eyear}.Check if the transitions file belong to another scenario.")
        
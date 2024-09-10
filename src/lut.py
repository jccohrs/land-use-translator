from pathlib import Path
import shutil
from cdo import *
from src.config import *
import os
import xarray as xr
import numpy as np
from src.utils import print_section_heading

import matplotlib.pyplot as plt
from scipy.stats import linregress
from matplotlib.patches import Patch
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap 

cdo = Cdo()

class LUT:

    def __init__(self, config):
        for key, value in config.items():
            setattr(self, key, value)
        self.grid = f"reg01_{self.region}"
        self.glc = f"{config.lcd}-{config.esayear}-{config.vers}"
        self.glc_lsm = f"{config.lcd}-2015-{config.vers}"
        self.namelist = self.generate_namelist()
        self.nr_crops = nr_crops
        self.nr_forest = nr_forest
        self.nr_grass = nr_grass
        self.nr_shrubs = nr_shrubs
        self.nr_urban = nr_urban
        self.years = abs(config.eyear - config.syear)

        # select period and self.region
        if self.region == 'Europe':
            self.reg = "-56,84,16,79"
        elif self.region == 'Australasia':
            self.reg = "102,218,-53,4"
        elif self.region == 'NorthAmerica':
            self.reg = "170,360,0,85"
        elif self.region == 'Germany':
            self.reg = "6,15.5,46.4,55.5"

        self.pfts_grass = GRAPFTS[0:self.nr_grass]
        self.pfts_crops = CROPFTS[0:self.nr_crops]
        self.pfts_shrubs = SHRPFTS[0:self.nr_shrubs]
        self.pfts_forest = FORPFTS[0:self.nr_forest]
        self.pfts_urban = URBPFTS[0:self.nr_urban]
        self.pft_grass_default = GRADEF
        self.pft_crops_default = CRODEF
        self.pft_shrubs_default = SHRDEF
        self.pft_forest_default = FORDEF
        self.pft_urban_default = URBDEF
        self.xsize = self.namelist["XSIZE"]
        self.ysize = self.namelist["YSIZE"]
        self.pft_frac_ts = np.zeros((self.xsize, self.ysize, self.npfts, self.years+1), dtype="float32")
        self.pfts_shrubs_grass = self.pfts_grass + self.pfts_shrubs

    def lucas_lut_forward(self):
        if self.addtree:
            rcm_lsm, pft_help, grass_backgr_help, crops_backgr_help, forest_backgr_help, shrubs_backgr_help, \
            shrubs_grass_backgr, nfv2cro, cro2nfv, for2cro, cro2for, ran2cro, cro2ran, pas2cro, cro2pas, cro2urb, nfv2urb, \
            for2urb, ran2urb, pas2urb, for2pas, pas2for, nfv2pas, ran2pas, pas2nfv, for2ran, ran2for, for2nfv, nfv2for, mcgrath_frac_help, \
            urb2cro, urb2nfv, urb2for, urb2ran, urb2pas, nfv2ran, nat2for, urban_backgr_help = self.lucas_lut_input()
        else:
            rcm_lsm, pft_help, grass_backgr_help, crops_backgr_help, forest_backgr_help, shrubs_backgr_help, \
            shrubs_grass_backgr, nfv2cro, cro2nfv, for2cro, cro2for, ran2cro, cro2ran, pas2cro, cro2pas, cro2urb, nfv2urb, \
            for2urb, ran2urb, pas2urb, for2pas, pas2for, nfv2pas, ran2pas, pas2nfv, for2ran, ran2for, for2nfv, nfv2for, mcgrath_frac_help, \
            urb2cro, urb2nfv, urb2for, urb2ran, urb2pas, nfv2ran, urban_backgr_help = self.lucas_lut_input()
        for z in range(self.years):
            pft_help = self.lucas_lut_transrules(for2cro[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_crops, self.pfts_forest, self.pfts_shrubs, self.pfts_grass, self.nr_crops, self.nr_forest, self.nr_shrubs, self.nr_grass, self.pft_crops_default, crops_backgr_help, False, 3, False)
            pft_help = self.lucas_lut_transrules(nfv2cro[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_crops, self.pfts_shrubs, self.pfts_grass, 0, self.nr_crops, self.nr_shrubs, self.nr_grass, 1, self.pft_crops_default, crops_backgr_help, False, 2, False)
            pft_help = self.lucas_lut_transrules(ran2cro[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_crops, self.pfts_shrubs, self.pfts_grass, 0, self.nr_crops, self.nr_shrubs, self.nr_grass, 1, self.pft_crops_default, crops_backgr_help, False, 2, False)
            pft_help = self.lucas_lut_transrules(pas2cro[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_crops, self.pfts_grass, 0, 0, self.nr_crops, self.nr_grass, 1, 1, self.pft_crops_default, crops_backgr_help, False, 1, False)
            pft_help = self.lucas_lut_transrules(cro2for[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_forest, self.pfts_crops, 0, 0, self.nr_forest, self.nr_crops, 1, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(cro2nfv[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_shrubs_grass, self.pfts_crops, 0, 0, self.nr_shrubs + self.nr_grass, self.nr_crops, 1, 1, self.pft_shrubs_default, shrubs_grass_backgr, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(cro2ran[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_grass, self.pfts_crops, 0, 0, self.nr_grass, self.nr_crops, 1, 1, self.pft_grass_default, grass_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(cro2pas[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_grass, self.pfts_crops, 0, 0, self.nr_grass, self.nr_crops, 1, 1, self.pft_grass_default, grass_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(cro2urb[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_urban, self.pfts_crops, 0, 0, self.nr_urban, self.nr_crops, 1, 1, self.pft_urban_default, urban_backgr_help, False, 1, False)
            pft_help = self.lucas_lut_transrules(for2urb[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_urban, self.pfts_forest, self.pfts_shrubs, self.pfts_grass, self.nr_urban, self.nr_forest, self.nr_shrubs, self.nr_grass, self.pft_urban_default, urban_backgr_help, False, 3, False)
            pft_help = self.lucas_lut_transrules(nfv2urb[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_urban, self.pfts_shrubs, self.pfts_grass, 0, self.nr_urban, self.nr_shrubs, self.nr_grass, 1, self.pft_urban_default, self.backgrd, 2, False)
            pft_help = self.lucas_lut_transrules(ran2urb[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_urban, self.pfts_shrubs_grass, 0, 0, self.nr_urban, self.nr_shrubs + self.nr_grass, 1, 1, self.pft_urban_default, urban_backgr_help, False, 1, False)
            pft_help = self.lucas_lut_transrules(pas2urb[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_urban, self.pfts_grass, 0, 0, self.nr_urban, self.nr_grass, 1, 1, self.pft_urban_default, urban_backgr_help, False, 1, False)
            pft_help = self.lucas_lut_transrules(urb2cro[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_crops, self.pfts_urban, 0, 0, self.nr_crops, self.nr_urban, 1, 1, self.pft_crops_default, crops_backgr_help, False, 1, False)
            pft_help = self.lucas_lut_transrules(urb2nfv[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_shrubs_grass, self.pfts_urban, 0, 0, self.nr_shrubs + self.nr_grass, self.nr_urban, 1, 1, self.pft_shrubs_default, shrubs_grass_backgr, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(urb2for[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_forest, self.pfts_urban, 0, 0, self.nr_forest, self.nr_urban, 1, 1, self.pft_forest_default, forest_backgr, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(urb2ran[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_shrubs_grass, self.pfts_urban, 0, 0, self.nr_shrubs + self.nr_grass, self.nr_urban, 1, 1, self.pft_shrubs_default, shrubs_grass_backgr, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(urb2pas[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_grass, self.pfts_urban, 0, 0, self.nr_grass, self.nr_urban, 1, 1, pft_grass_default, grass_backgr[i, j, :], self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(for2pas[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_grass, self.pfts_forest, 0, 0, self.nr_grass, self.nr_forest, 1, 1, self.pft_grass_default, grass_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(nfv2pas[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_grass, self.pfts_shrubs, 0, 0, self.nr_grass, self.nr_shrubs, 1, 1, self.pft_grass_default, grass_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(ran2pas[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_grass, self.pfts_shrubs, 0, 0, self.nr_grass, self.nr_shrubs, 1, 1, self.pft_grass_default, grass_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(pas2for[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_forest, self.pfts_grass, 0, 0, self.nr_forest, self.nr_grass, 1, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(pas2nfv[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_shrubs, self.pfts_grass, 0, 0, self.nr_shrubs, self.nr_grass, 1, 1, self.pft_shrubs_default, shrubs_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(for2ran[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_shrubs_grass, self.pfts_forest, 0, 0, self.nr_shrubs + self.nr_grass, self.nr_forest, 1, 1, self.pft_shrubs_default, shrubs_grass_backgr, self.backgrd, 1, False)
            #pft_help = self.lucas_lut_transrules(nfv2ran[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_grass, self.pfts_shrubs, 0, 0, self.nr_grass, self.nr_shrubs, 1, 1, self.pft_grass_default, grass_backgr, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(ran2for[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_forest, self.pfts_shrubs, self.pfts_grass, 0, self.nr_forest, self.nr_shrubs, self.nr_grass, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 2, False)
            #pft_help = self.lucas_lut_transrules(for2nfv[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_shrubs_grass, self.pfts_forest, 0, 0, self.nr_shrubs + self.nr_grass, self.nr_forest, 1, 1, self.pft_shrubs_default, shrubs_grass_backgr, self.backgrd, 1, False)
            #pft_help = self.lucas_lut_transrules(nfv2for[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_forest, self.pfts_shrubs, self.pfts_grass, 0, self.nr_forest, self.nr_shrubs, self.nr_grass, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 2, False)
            if self.addtree:
                self.lucas_lut_transrules(nat2for[z, :, :].data.T, rcm_lsm, pft_help, self.pfts_forest, self.pfts_shrubs, self.pfts_grass, 0, self.nr_forest, self.nr_shrubs, self.nr_grass, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 2, False)
            self.pft_frac_ts[:, :, :, z] = self.pft_help[:, :, :]
        print('land use change finished')
        # NORMALIZE TO GET A SUM OF 1 AND SET SEA POINTS TO MISSING VALUE
        self.recalc_pft_frac_ts(rcm_lsm)
        print("Finished normalizing")
        print('LAND USE CHANGE FINISHED')
        if self.irri:
            print_section_heading('IRRIGATION')
            self.lucas_lut_irrigation(rcm_lsm)
            self.recalc_null_pft_frac_ts(rcm_lsm)
        if self.mcgrath:
            print_section_heading('MCGRATH')
            self.lucas_lut_mcgrath(rcm_lsm, mcgrath_frac_help)
            self.recalc_null_pft_frac_ts(rcm_lsm)
        if self.plot:
            print_section_heading('PLOTTING')
            self.plot_pft_frac_ts()
            print('PLOTTING FINISHED')

    def lucas_lut_backward(self):
        rcm_lsm, pft_help, grass_backgr_help, crops_backgr_help, forest_backgr_help, shrubs_backgr_help, \
        shrubs_grass_backgr, nfv2cro, cro2nfv, for2cro, cro2for, ran2cro, cro2ran, pas2cro, cro2pas, cro2urb, nfv2urb, \
        for2urb, ran2urb, pas2urb, for2pas, pas2for, nfv2pas, ran2pas, pas2nfv, for2ran, ran2for, for2nfv, nfv2for, mcgrath_frac_help = self.lucas_lut_input()
        for z in range(self.years):
            zz = self.years - z -1
            print('year', zz)
            pft_help = self.lucas_lut_transrules(nfv2cro[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_shrubs, self.pfts_crops, 0, 0, self.nr_shrubs, self.nr_crops, 1, 1, self.pft_shrubs_default, shrubs_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(cro2nfv[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_crops, self.pfts_shrubs, self.pfts_grass, 0, self.nr_crops, self.nr_shrubs, self.nr_grass, 1, self.pft_crops_default, crops_backgr_help, self.backgrd, 2, False)
            pft_help = self.lucas_lut_transrules(for2cro[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_forest, self.pfts_crops, 0, 0, self.nr_forest, self.nr_crops, 1, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 1, self.mcgrath, mcgfrac=mcgrath_frac_help[:, :, :, zz])
            pft_help = self.lucas_lut_transrules(cro2for[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_crops, self.pfts_forest, self.pfts_shrubs, 0, self.nr_crops, self.nr_forest, self.nr_shrubs, 1, self.pft_crops_default, crops_backgr_help, False, 2, False)
            pft_help = self.lucas_lut_transrules(ran2cro[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_grass, self.pfts_crops, 0, 0, self.nr_grass, self.nr_crops, 1, 1, self.pft_grass_default, grass_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(cro2ran[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_crops, self.pfts_grass, 0, 0, self.nr_crops, self.nr_grass, 1, 1, self.pft_crops_default, crops_backgr_help, False, 1, False)
            pft_help = self.lucas_lut_transrules(pas2cro[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_grass, self.pfts_crops, 0, 0, self.nr_grass, self.nr_crops, 1, 1, self.pft_grass_default, grass_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(cro2pas[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_crops, self.pfts_grass, 0, 0, self.nr_crops, self.nr_grass, 1, 1, self.pft_crops_default, crops_backgr_help, False, 1, False)
            pft_help = self.lucas_lut_transrules(cro2urb[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_crops, self.pfts_urban, 0, 0, self.nr_crops, self.nr_urban, 1, 1, self.pft_crops_default, crops_backgr_help, False, 1, False)
            pft_help = self.lucas_lut_transrules(nfv2urb[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_shrubs, self.pfts_urban, 0, 0, self.nr_shrubs, self.nr_urban, 1, 1, self.pft_shrubs_default, shrubs_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(for2urb[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_forest, self.pfts_urban, 0, 0, self.nr_forest, self.nr_urban, 1, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 1, self.mcgrath, mcgfrac=mcgrath_frac_help[:, :, :, zz])
            pft_help = self.lucas_lut_transrules(ran2urb[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_shrubs_grass, self.pfts_urban, 0, 0, self.nr_shrubs+self.nr_grass, self.nr_urban, 1, 1, self.pft_shrubs_default, shrubs_grass_backgr, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(pas2urb[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_grass, self.pfts_urban, 0, 0, self.nr_grass, self.nr_urban, 1, 1, self.pft_grass_default, grass_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(for2pas[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_forest, self.pfts_grass, 0, 0, self.nr_forest, self.nr_grass, 1, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 1, self.mcgrath, mcgfrac=mcgrath_frac_help[:, :, :, zz])
            pft_help = self.lucas_lut_transrules(pas2for[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_grass, self.pfts_forest, 0, 0, self.nr_grass, self.nr_forest, 1, 1, self.pft_grass_default, grass_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(nfv2pas[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_shrubs, self.pfts_grass, 0, 0, self.nr_shrubs, self.nr_grass, 1, 1, self.pft_shrubs_default, shrubs_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(ran2pas[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_shrubs, self.pfts_grass, 0, 0, self.nr_shrubs, self.nr_grass, 1, 1, self.pft_shrubs_default, shrubs_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(pas2nfv[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_grass, self.pfts_shrubs, 0, 0, self.nr_grass, self.nr_shrubs, 1, 1, self.pft_grass_default, grass_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(for2ran[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_forest, self.pfts_shrubs, self.pfts_grass, 0, self.nr_forest, self.nr_shrubs, self.nr_grass, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 2, self.mcgrath, mcgfrac=mcgrath_frac_help[:, :, :, zz])
            pft_help = self.lucas_lut_transrules(ran2for[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_shrubs_grass, self.pfts_forest, 0, 0, self.nr_grass+self.nr_shrubs, self.nr_forest, 1, 1, self.pft_shrubs_default, shrubs_grass_backgr, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(for2nfv[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_forest, self.pfts_shrubs, self.pfts_grass, 0, self.nr_forest, self.nr_shrubs, self.nr_grass, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 2, self.mcgrath, mcgfrac=mcgrath_frac_help[:, :, :, zz])
            pft_help = self.lucas_lut_transrules(nfv2for[zz, :, :].data.T, rcm_lsm, pft_help, self.pfts_shrubs_grass, self.pfts_forest, 0, 0, self.nr_grass+self.nr_shrubs, self.nr_forest, 1, 1, self.pft_shrubs_default, shrubs_grass_backgr, self.backgrd, 1, False)
            self.pft_frac_ts[:, :, :, zz] = pft_help[:, :, :]
        # NORMALIZE TO GET A SUM OF 1 AND SET SEA POINTS TO MISSING VALUE
        self.recalc_pft_frac_ts(rcm_lsm)
        if self.irri:
            print_section_heading('IRRIGATION')
            self.lucas_lut_irrigation(rcm_lsm)
            self.recalc_null_pft_frac_ts(rcm_lsm)
            self.recalc_pft_frac_ts(rcm_lsm)
        if self.mcgrath:
            print_section_heading('MCGRATH')
            self.lucas_lut_mcgrath(rcm_lsm, mcgrath_frac_help)
            self.recalc_null_pft_frac_ts(rcm_lsm)
            self.recalc_pft_frac_ts(rcm_lsm)
        if self.plot:
            print_section_heading('PLOTTING')
            self.plot_pft_frac_ts()

    def lucas_lut_transrules(self, trans, rcm_lsm, inpfts, pfts_1, pfts_2, pfts_3, pfts_4,
                             nr_pfts_1, nr_pfts_2, nr_pfts_3, nr_pfts_4, defaultpft,
                             backgrdpfts, backgrd, rule, mcgrath, mcgfrac=np.zeros((3), dtype="float32")):

        # SUM VARIABLES
        pft_1_sum = sum(inpfts[:, :, pfts_1[n]-1] for n in range(nr_pfts_1))
        pft_2_sum = sum(inpfts[:, :, pfts_2[n]-1] for n in range(nr_pfts_2))
        pft_3_sum = sum(inpfts[:, :, pfts_3[n]-1] for n in range(nr_pfts_3)) if rule >= 2 else 0
        pft_4_sum = sum(inpfts[:, :, pfts_4[n]-1] for n in range(nr_pfts_4)) if rule >= 3 else 0
        mcg_sum = sum(mcgfrac[:, :, n] for n in range(3)) if mcgrath else 0
        # HELP VARIABLES
        helper = np.zeros((self.xsize, self.ysize), dtype="float32")
        helper_2 = np.zeros((self.xsize, self.ysize), dtype="float32")
        helper_3 = np.zeros((self.xsize, self.ysize), dtype="float32")
        # STATIC VALUES TO RESTORE LATER
        inpfts_static_values = np.where(rcm_lsm[..., np.newaxis] > 0.0, -998, inpfts)
        inpfts_trans_static_values = np.where(trans[..., np.newaxis] > 0.0, -999, inpfts)
        # RULE 1 : subtract one group and increas one group
        if rule == 1:
            var1 = trans - np.maximum(0.0, trans + pft_1_sum - 1)
            var2 = trans - np.maximum(0.0, trans - pft_2_sum)
            # limit the transition so that TRANS+PFT1 is less equal 1 and TRANS-PFT2 is greater equal 0
            trans = np.where(var1 < var2, var1, var2)
            mask_trans = trans > 0.0
            # filtering pft_2_sum and trans values to avoid division by zero
            filtered_pft_2_sum = np.where((mask_trans) & (pft_2_sum > 0.0), pft_2_sum, 1.0)
            filtered_trans = np.where((mask_trans) & (pft_2_sum > 0.0), trans, 0.0)
            
            # subtracting from PFT group 2
            for ipft in range(nr_pfts_2):
                inpfts[:, :, pfts_2[ipft]-1] = np.where((mask_trans) & (pft_2_sum > 0.0), inpfts[:, :, pfts_2[ipft]-1] - (inpfts[:, :, pfts_2[ipft]-1] / filtered_pft_2_sum * filtered_trans), inpfts[:, :, pfts_2[ipft]-1])
                # if more fraction is removed than available set fraction to zero
                helper = np.where((inpfts[:, :, pfts_2[ipft]-1] < 0.0) & (pft_2_sum > 0.0) & (mask_trans), helper + inpfts[:, :, pfts_2[ipft]-1], helper)
                inpfts[:, :, pfts_2[ipft]-1] = np.where((inpfts[:, :, pfts_2[ipft]-1] < 0.0) & (pft_2_sum > 0.0) & (mask_trans), 0.0, inpfts[:, :, pfts_2[ipft]-1])
            
            # adjust transition if needed (not enough of PFT group 2), helper is always <= 0
            trans = np.where((mask_trans), np.maximum(0.0, trans + helper), trans)
            # if McGrath forest data should be used
            mask_mcg_sum = mcg_sum == 1
            if mcgrath:
                for ipft in range(2, 5):
                    inpfts[:, :, pfts_1[ipft]-1] = np.where((mcgrath) & (mask_trans) & (mask_mcg_sum), inpfts[:, :, pfts_1[ipft]-1] + (mcgfrac[:, :, ipft-2] * trans), inpfts[:, :, pfts_1[ipft]-1])
            # just use the relative fractions
            filtered_pft_1_sum = np.where((mask_trans) & (pft_1_sum > 0.0) & ((~mcgrath) | (~mask_mcg_sum)), pft_1_sum, 1)
            filtered_trans = np.where((mask_trans) & (pft_1_sum > 0.0) & ((~mcgrath) | (~mask_mcg_sum)), trans, 0)
            for ipft in range(nr_pfts_1):
                inpfts[:, :, pfts_1[ipft]-1] = np.where((mask_trans) & (pft_1_sum > 0.0) & ((~mcgrath) | (~mask_mcg_sum)), inpfts[:, :, pfts_1[ipft]-1] + (inpfts[:, :, pfts_1[ipft]-1] / filtered_pft_1_sum * filtered_trans), inpfts[:, :, pfts_1[ipft]-1])
            if backgrd:
                for ipft in range(nr_pfts_1):
                    inpfts[:, :, pfts_1[ipft]-1] = np.where((pft_1_sum <= 0.0) & (mask_trans) & ((~mcgrath) | (~mask_mcg_sum)), inpfts[:, :, pfts_1[ipft]-1] + (backgrdpfts[:, :, ipft] * trans), inpfts[:, :, pfts_1[ipft]-1])        
            else:
                inpfts[:, :, defaultpft-1] = np.where((pft_1_sum <= 0.0) & (mask_trans) & ((~mcgrath) | (~mask_mcg_sum)), trans, inpfts[:, :, defaultpft-1])
        # RULE 2
        elif rule == 2:
            # limit the transition so that TRANS+PFT1 is less equal 1 and TRANS-(PFT2+PFT3) is greater equal 0
            var1 = trans - np.maximum(0.0, pft_1_sum + trans - 1)
            var2 = trans - np.maximum(0.0, trans - (pft_2_sum + pft_3_sum))
            trans = np.where(var1 < var2, var1, var2)
            mask_trans = trans > 0.0
            # filtering pft_2_sum and trans values to avoid division by zero
            filtered_pft_2_sum = np.where((mask_trans) & (pft_2_sum > 0.0), pft_2_sum, 1.0)
            filtered_trans = np.where((mask_trans) & (pft_2_sum > 0.0), trans, 0.0)

            # subtracting from PFT group 2
            for ipft in range(nr_pfts_2):
                inpfts[:, :, pfts_2[ipft]-1] = np.where((mask_trans) & (pft_2_sum > 0.0), inpfts[:, :, pfts_2[ipft]-1] - (inpfts[:, :, pfts_2[ipft]-1] / filtered_pft_2_sum * filtered_trans), inpfts[:, :, pfts_2[ipft]-1])
                
                # if more fraction is removed than available set fraction to zero
                helper = np.where((inpfts[:, :, pfts_2[ipft]-1] < 0.0) & (pft_2_sum > 0.0) & (mask_trans), helper - inpfts[:, :, pfts_2[ipft]-1], helper)
                inpfts[:, :, pfts_2[ipft]-1] = np.where((inpfts[:, :, pfts_2[ipft]-1] < 0.0) & (pft_2_sum > 0.0) & (mask_trans), 0.0, inpfts[:, :, pfts_2[ipft]-1])
                
            helper = np.where((mask_trans) & (pft_2_sum <= 0.0), trans, helper)
            mask_helper = helper > 0.0

            # filtering pft_3_sum and helper values to avoid division by zero
            filtered_pft_3_sum = np.where((mask_trans) & (mask_helper) & (pft_3_sum > 0.0), pft_3_sum, 1.0)
            filtered_helper = np.where((mask_trans) & (mask_helper) & (pft_3_sum > 0.0), helper, 0.0)
            for ipft in range(nr_pfts_3):
                inpfts[:, :, pfts_3[ipft]-1] = np.where((mask_trans) & (mask_helper) & (pft_3_sum > 0.0), inpfts[:, :, pfts_3[ipft]-1] - (inpfts[:, :, pfts_3[ipft]-1] / filtered_pft_3_sum * filtered_helper), inpfts[:, :, pfts_3[ipft]-1])
                helper_2 = np.where((inpfts[:, :, pfts_3[ipft]-1] < 0.0) & (pft_3_sum > 0.0) & (mask_trans) & (mask_helper), helper_2 + inpfts[:, :, pfts_3[ipft]-1], helper_2)
                inpfts[:, :, pfts_3[ipft]-1] = np.where((inpfts[:, :, pfts_3[ipft]-1] < 0.0) & (pft_3_sum > 0.0) & (mask_trans) & (mask_helper), 0.0, inpfts[:, :, pfts_3[ipft]-1])
            
            trans = np.where(mask_trans, np.maximum(0.0, trans + helper_2, dtype="float32"), trans)
            # if McGrath forest data should be used
            mask_mcg_sum = mcg_sum == 1
            if mcgrath:
                for ipft in range(2, 5):
                    inpfts[:, :, pfts_1[ipft]-1] = np.where((mask_trans) & (mask_mcg_sum) & (mcgrath), inpfts[:, :, pfts_1[ipft]-1] + (mcgfrac[:, :, ipft-2] * trans), inpfts[:, :, pfts_1[ipft]-1])
            # just use the relative fractions
            filtered_pft_1_sum = np.where((mask_trans) & (pft_1_sum > 0.0) & ((~mcgrath) | (~mask_mcg_sum)), pft_1_sum, 1)
            filtered_trans = np.where((mask_trans) & (pft_1_sum > 0.0) & ((~mcgrath) | (~mask_mcg_sum)), trans, 0)
            for ipft in range(nr_pfts_1):
                inpfts[:, :, pfts_1[ipft]-1] = np.where((mask_trans) & (pft_1_sum > 0.0) & ((~mcgrath) | (~mask_mcg_sum)), inpfts[:, :, pfts_1[ipft]-1] + (inpfts[:, :, pfts_1[ipft]-1] / filtered_pft_1_sum * filtered_trans), inpfts[:, :, pfts_1[ipft]-1])
            if backgrd:
                for ipft in range(nr_pfts_1):
                    inpfts[:, :, pfts_1[ipft]-1] = np.where((pft_1_sum <= 0.0) & (mask_trans) & ((~mcgrath) | (~mask_mcg_sum)), inpfts[:, :, pfts_1[ipft]-1] + (backgrdpfts[:, :, ipft] * trans), inpfts[:, :, pfts_1[ipft]-1])
            else:
                inpfts[:, :, defaultpft-1] = np.where((pft_1_sum <= 0.0) & (mask_trans) & ((~mcgrath) | (~mask_mcg_sum)), trans, inpfts[:, :, defaultpft-1])
        # RULE 3
        elif rule == 3:
            # limit the transition so that TRANS+PFT1 is less equal 1 and TRANS-(PFT2+PFT3) is greater equal 0
            var1 = trans - np.maximum(0., pft_1_sum + trans - 1.0)
            var2 = trans - np.maximum(0., trans - (pft_2_sum + pft_3_sum + pft_4_sum))
            trans = np.where(var1 < var2, var1, var2)
            mask_trans = trans > 0.0

            # filtering to avoid division by zero
            filtered_pft_2_sum = np.where((mask_trans) & (pft_2_sum > 0.0), pft_2_sum, 1.0)
            filtered_trans = np.where((mask_trans) & (pft_2_sum > 0.0), trans, 0.0)
            filtered_pft_2_sum = np.where(filtered_trans == 0.0, 1.0, filtered_pft_2_sum)
            
            # subtracting from PFT group 2
            for ipft in range(nr_pfts_2):
                inpfts[:, :, pfts_2[ipft-1]] -= (inpfts[:, :, pfts_2[ipft-1]]/filtered_pft_2_sum*filtered_trans)
                inpfts_pfts_2 = inpfts[:, :, pfts_2[ipft-1]]
                
                # if more fraction is removed than available set fraction to zero
                helper = np.where((mask_trans) & (inpfts[:, :, pfts_2[ipft-1]] < 0.0) & (pft_2_sum > 0.0), helper - inpfts[:, :, pfts_2[ipft-1]], helper)
                inpfts[:, :, pfts_2[ipft-1]] = np.where((mask_trans) & (inpfts[:, :, pfts_2[ipft-1]] < 0.0) & (pft_2_sum > 0.0), 0.0, inpfts[:, :, pfts_2[ipft-1]])
            
            helper = np.where((mask_trans) & (pft_2_sum <= 0.0), trans, helper)
            mask_helper = helper > 0.0

            # filtering to avoid division by zero
            filtered_pft_3_sum = np.where((mask_trans) & (mask_helper) & (pft_3_sum > 0.0), pft_3_sum, 1.0)
            filtered_helper = np.where((mask_trans) & (mask_helper) & (pft_3_sum > 0.0), helper, 0.0)
            filtered_pft_3_sum = np.where(filtered_helper == 0.0, 1.0, filtered_pft_3_sum)
            
            # subtracting from PFT group 3
            for ipft in range(nr_pfts_3):
                inpfts[:, :, pfts_3[ipft]-1] -= (inpfts[:, :, pfts_3[ipft]-1]/filtered_pft_3_sum*filtered_helper)

                # if more fraction is removed than available set fraction to zero
                helper_2 = np.where((mask_trans) & (mask_helper) & (pft_3_sum > 0.0) & (inpfts[:, :, pfts_3[ipft]] < 0.0), helper_2 - inpfts[:, :, pfts_3[ipft]], helper_2)
                inpfts[:, :, pfts_3[ipft]-1] = np.where((mask_trans) & (mask_helper) & (pft_3_sum > 0.0) & (inpfts[:, :, pfts_3[ipft]] < 0.0), 0.0, inpfts[:, :, pfts_3[ipft]])
            
            helper_2 = np.where((mask_trans) & (mask_helper) & (pft_3_sum <= 0.0), helper, helper_2)
            mask_helper_2 = helper_2 > 0.0

            # subtracting from PFT group 3
            filtered_pft_4_sum = np.where((mask_trans) & (mask_helper_2) & (pft_4_sum > 0.0), pft_4_sum, 1.0)
            filtered_helper_2 = np.where((mask_trans) & (mask_helper_2) & (pft_4_sum > 0.0), helper_2, 0.0)
            filtered_pft_4_sum = np.where(filtered_helper_2 == 0.0, 1.0, filtered_pft_4_sum)
            
            # adding from PFT group 4
            for ipft in range(nr_pfts_4):
                inpfts[:, :, pfts_4[ipft]-1] -= (inpfts[:, :, pfts_4[ipft]-1]/filtered_pft_4_sum*filtered_helper_2)

                # if more fraction is removed than available set fraction to zero
                helper_3 = np.where((mask_trans) & (mask_helper_2) & (pft_4_sum > 0.0) & (inpfts[:, :, pfts_4[ipft]-1] < 0.0), helper_3 + inpfts[:, :, pfts_4[ipft]-1], helper_3)
                inpfts[:, :, pfts_4[ipft]-1] = np.where((mask_trans) & (mask_helper_2) & (pft_4_sum > 0.0) & (inpfts[:, :, pfts_4[ipft]-1] < 0.0), 0.0, inpfts[:, :, pfts_4[ipft]-1])

            trans = np.where((mask_trans) & (pft_4_sum <= 0.0), np.maximum(0.0, trans + helper_3), trans)
            
            # if McGrath forest data should be used
            mask_mcg_sum = mcg_sum == 1
            if mcgrath:
                for ipft in range(2, 5):
                    inpfts[:, :, pfts_1[ipft]-1] = np.where((mask_trans) & (mask_mcg_sum) & (mcgrath), inpfts[:, :, pfts_1[ipft]-1] + (mcgfrac[:, :, ipft-2] * trans), inpfts[:, :, pfts_1[ipft]-1])
            filtered_pft_1_sum = np.where((mask_trans) & (pft_1_sum > 0.0) & ((~mcgrath) | (~mask_mcg_sum)), pft_1_sum, 1.0)
            filtered_trans = np.where((mask_trans) & (pft_1_sum > 0.0) & ((~mcgrath) | (~mask_mcg_sum)), trans, 0.0)
            filtered_pft_1_sum = np.where(filtered_trans == 0.0, 1.0, filtered_pft_1_sum)
            for ipft in range(nr_pfts_1):
                inpfts[:, :, pfts_1[ipft]-1] += (inpfts[:, :, pfts_1[ipft]-1]/filtered_pft_1_sum*filtered_trans)
            if backgrd:
                for ipft in range(nr_pfts_1):
                    inpfts[:, :, pfts_1[ipft]-1] = np.where((mask_trans) & (pft_1_sum > 0.0) & ((~mcgrath) | (~mask_mcg_sum)), inpfts[:, :, pfts_1[ipft]-1], inpfts[:, :, pfts_1[ipft]-1] + (backgrdpfts[:, :, ipft] * trans))
            else:
                inpfts[:, :, defaultpft-1] = np.where((mask_trans) & (pft_1_sum > 0.0) & ((~mcgrath) | (~mask_mcg_sum)), inpfts[:, :, defaultpft-1], trans)
        
        # RESTORE STATIC VALUES
        inpfts = np.where(inpfts_trans_static_values == -999, inpfts, inpfts_trans_static_values)
        inpfts = np.where(inpfts_static_values == -998, inpfts, inpfts_static_values)
        return inpfts

    def lucas_lut_help(self):
        data_help = np.zeros((self.xsize, self.ysize, self.npfts, self.years+1))
        data = xr.open_dataset("data/LUCAS_LUC7_ESACCI_LUH2_historical_1950_2015_reg01_Germany_mcgrath.nc")
        for i in range(self.npfts):
            num = i + 1
            if len(str(num)) > 1:
                data_help[:, :, i, :] = data[f"var8{num}"].data.T
            else:
                data_help[:, :, i, :] = data[f"var80{num}"].data.T
        return data_help

    def lucas_lut_input(self):
        pft_help = np.zeros((self.xsize, self.ysize, self.npfts), dtype='float32')
        self.pft_frac = xr.open_dataset(self.namelist["F_LC_IN"]).isel(time=0)
        for i in range(self.npfts):
            num = i + 1
            if len(str(num)) > 1:
                pft_help[:, :, i] = self.pft_frac[f"var8{num}"].data.T
                if self.forward:
                    self.pft_frac_ts[:, :, i, 0] = self.pft_frac[f"var8{num}"].data.T
                else:
                    self.pft_frac_ts[:, :, i, self.years] = self.pft_frac[f"var8{num}"].data.T
            else:
                pft_help[:, :, i] = self.pft_frac[f"var80{num}"].data.T
                if self.forward:
                    self.pft_frac_ts[:, :, i, 0] = self.pft_frac[f"var80{num}"].data.T
                else:
                    self.pft_frac_ts[:, :, i, self.years] = self.pft_frac[f"var80{num}"].data.T
        if self.forward:
            urban_backgr = xr.open_dataset(self.namelist["F_BACKURB"])
            urban_backgr_help = np.zeros((self.xsize, self.ysize, self.nr_urban), dtype='float32')
            for i in range(self.nr_urban):
                num = i + 15
                urban_backgr_help[:, :, i] = urban_backgr[f"var8{num}"].data[0, :, :].T
        grass_backgr = xr.open_dataset(self.namelist["F_BACKGRA"])
        grass_backgr_help = np.zeros((self.xsize, self.ysize, self.nr_grass), dtype='float32')
        for i in range(self.nr_grass):
            num = i + 9
            if len(str(num)) > 1:
                grass_backgr_help[:, :, i] = grass_backgr[f"var8{num}"].data[0, :, :].T
            else:
                grass_backgr_help[:, :, i] = grass_backgr[f"var80{num}"].data[0, :, :].T

        crops_backgr = xr.open_dataset(self.namelist["F_BACKCRO"])
        crops_backgr_help = np.zeros((self.xsize, self.ysize, self.nr_crops), dtype='float32')
        for i in range(self.nr_crops):
            num = i + 13
            if len(str(num)) > 1:
                crops_backgr_help[:, :, i] = crops_backgr[f"var8{num}"].data[0, :, :].T
            else:
                crops_backgr_help[:, :, i] = crops_backgr[f"var80{num}"].data[0, :, :].T
        forest_backgr = xr.open_dataset(self.namelist["F_BACKFOR"])
        forest_backgr_help = np.zeros((self.xsize, self.ysize, self.nr_forest), dtype='float32')
        for i in range(self.nr_forest):
            num = i + 1
            if len(str(num)) > 1:
                forest_backgr_help[:, :, i] = forest_backgr[f"var8{num}"].data[0, :, :].T
            else:
                forest_backgr_help[:, :, i] = forest_backgr[f"var80{num}"].data[0, :, :].T

        shrubs_backgr = xr.open_dataset(self.namelist["F_BACKSHR"])
        shrubs_backgr_help = np.zeros((self.xsize, self.ysize, self.nr_shrubs), dtype='float32')
        for i in range(self.nr_shrubs):
            num = i + 7
            if len(str(num)) > 1:
                shrubs_backgr_help[:, :, i] = shrubs_backgr[f"var8{num}"].data[0, :, :].T
            else:
                shrubs_backgr_help[:, :, i] = shrubs_backgr[f"var80{num}"].data[0, :, :].T
    
        shrubs_grass_backgr = np.zeros((self.xsize, self.ysize, self.nr_shrubs+self.nr_grass), dtype='float32')
        shrubs_grass_backgr[:, :, :self.nr_shrubs] = shrubs_backgr_help / 2.
        shrubs_grass_backgr[:, :, self.nr_shrubs:(self.nr_grass+self.nr_shrubs)] = grass_backgr_help / 2.
        # Mcfrac dataset
        mcgrath_frac = xr.open_dataset(self.namelist["F_MCGRATH"], decode_times=False)
        mcgrath_frac_help = np.zeros((self.xsize, self.ysize, 3, self.years+1), dtype='float32')
        for i in range(1, 4):
            mcgrath_frac_help[:, :, i-1, :] = mcgrath_frac[f"var80{i+2}"][:, :, :].data.T
        # RCM LSM
        rcm_lsm = xr.open_dataset(self.namelist["F_RCM_LSM_IN"]).var210.values.T
        # Transformation datasets
        nfv2cro = xr.open_dataset(self.namelist["F_NFV2CRO"], decode_times=False)["nfv2cro"]
        cro2nfv = xr.open_dataset(self.namelist["F_CRO2NFV"], decode_times=False)["cro2nfv"]
        for2cro = xr.open_dataset(self.namelist["F_FOR2CRO"], decode_times=False)["for2cro"]
        cro2for = xr.open_dataset(self.namelist["F_CRO2FOR"], decode_times=False)["cro2for"]
        ran2cro = xr.open_dataset(self.namelist["F_RAN2CRO"], decode_times=False)["ran2cro"]
        cro2ran = xr.open_dataset(self.namelist["F_CRO2RAN"], decode_times=False)["cro2ran"]
        pas2cro = xr.open_dataset(self.namelist["F_PAS2CRO"], decode_times=False)["pas2cro"]
        cro2pas = xr.open_dataset(self.namelist["F_CRO2PAS"], decode_times=False)["cro2pas"]
        cro2urb = xr.open_dataset(self.namelist["F_CRO2URB"], decode_times=False)["cro2urb"]
        nfv2urb = xr.open_dataset(self.namelist["F_NFV2URB"], decode_times=False)["nfv2urb"]
        for2urb = xr.open_dataset(self.namelist["F_FOR2URB"], decode_times=False)["for2urb"]
        ran2urb = xr.open_dataset(self.namelist["F_RAN2URB"], decode_times=False)["ran2urb"]
        pas2urb = xr.open_dataset(self.namelist["F_PAS2URB"], decode_times=False)["pas2urb"]
        for2pas = xr.open_dataset(self.namelist["F_FOR2PAS"], decode_times=False)["for2pas"]
        pas2for = xr.open_dataset(self.namelist["F_PAS2FOR"], decode_times=False)["pas2for"]
        nfv2pas = xr.open_dataset(self.namelist["F_NFV2PAS"], decode_times=False)["nfv2pas"]
        ran2pas = xr.open_dataset(self.namelist["F_RAN2PAS"], decode_times=False)["ran2pas"]
        # CHECK variable name
        pas2nfv = xr.open_dataset(self.namelist["F_PAS2NFV"], decode_times=False)["pas2nfv"]
        for2ran = xr.open_dataset(self.namelist["F_FOR2RAN"], decode_times=False)["for2ran"]
        ran2for = xr.open_dataset(self.namelist["F_RAN2FOR"], decode_times=False)["ran2for"]
        for2nfv = xr.open_dataset(self.namelist["F_FOR2NFV"], decode_times=False)["for2nfv"]
        nfv2for = xr.open_dataset(self.namelist["F_NFV2FOR"], decode_times=False)["nfv2for"]
        if self.forward:
            urb2cro = xr.open_dataset(self.namelist["F_URB2CRO"], decode_times=False)["urb2cro"]
            urb2nfv = xr.open_dataset(self.namelist["F_URB2NFV"], decode_times=False)["urb2nfv"]
            urb2for = xr.open_dataset(self.namelist["F_URB2FOR"], decode_times=False)["urb2for"]
            urb2ran = xr.open_dataset(self.namelist["F_URB2RAN"], decode_times=False)["urb2ran"]
            urb2pas = xr.open_dataset(self.namelist["F_URB2PAS"], decode_times=False)["urb2pas"]
            nfv2ran = xr.open_dataset(self.namelist["F_NFV2RAN"], decode_times=False)["nfv2ran"]
            if self.addtree:
                nat2for = xr.open_dataset(self.namelist["F_NAT2FOR"], decode_times=False)["nat2for"]
                return rcm_lsm, pft_help, grass_backgr_help, crops_backgr_help, forest_backgr_help, shrubs_backgr_help, \
                       shrubs_grass_backgr, nfv2cro, cro2nfv, for2cro, cro2for, ran2cro, cro2ran, pas2cro, cro2pas, cro2urb, nfv2urb, \
                       for2urb, ran2urb, pas2urb, for2pas, pas2for, nfv2pas, ran2pas, pas2nfv, for2ran, ran2for, for2nfv, nfv2for, mcgrath_frac_help, \
                       urb2cro, urb2nfv, urb2for, urb2ran, urb2pas, nfv2ran, nat2for, urban_backgr_help
            
            return rcm_lsm, pft_help, grass_backgr_help, crops_backgr_help, forest_backgr_help, shrubs_backgr_help, \
               shrubs_grass_backgr, nfv2cro, cro2nfv, for2cro, cro2for, ran2cro, cro2ran, pas2cro, cro2pas, cro2urb, nfv2urb, \
               for2urb, ran2urb, pas2urb, for2pas, pas2for, nfv2pas, ran2pas, pas2nfv, for2ran, ran2for, for2nfv, nfv2for, mcgrath_frac_help, \
                urb2cro, urb2nfv, urb2for, urb2ran, urb2pas, nfv2ran, urban_backgr_help

        return rcm_lsm, pft_help, grass_backgr_help, crops_backgr_help, forest_backgr_help, shrubs_backgr_help, \
               shrubs_grass_backgr, nfv2cro, cro2nfv, for2cro, cro2for, ran2cro, cro2ran, pas2cro, cro2pas, cro2urb, nfv2urb, \
               for2urb, ran2urb, pas2urb, for2pas, pas2for, nfv2pas, ran2pas, pas2nfv, for2ran, ran2for, for2nfv, nfv2for, mcgrath_frac_help

    def lucas_lut_irrigation(self, rcm_lsm):
        irri_frac = xr.open_dataset(self.namelist["F_IRRI_IN"], decode_times=False)["irrig_frac"].data.T
        for z in range(self.years+1):
            mask = (rcm_lsm > 0.0) & ((self.pft_frac_ts[:, :, 12, z] + self.pft_frac_ts[:, :, 13, z]) > 0.0)
            sum_crops = self.pft_frac_ts[:, :, 12, z] + self.pft_frac_ts[:, :, 13, z]
            irri_mask = irri_frac[:, :, z] > 0.0
            self.pft_frac_ts[mask & irri_mask, 12, z] = (1.0 - irri_frac[mask & irri_mask, z]) * sum_crops[mask & irri_mask]
            self.pft_frac_ts[mask & irri_mask, 13, z] = irri_frac[mask & irri_mask, z] * sum_crops[mask & irri_mask]
            self.pft_frac_ts[mask & ~irri_mask, 12, z] = sum_crops[mask & ~irri_mask]
            mask = rcm_lsm > 0.0
            self.pft_frac_ts[:, :, :, z] = np.where(self.pft_frac_ts[:, :, :, z] < 0.0, 0.0, self.pft_frac_ts[:, :, :, z])
            pft_sum = np.sum(self.pft_frac_ts[:, :, :, z], axis=2)
            pft_sum_mask = pft_sum > 0.0
            self.pft_frac_ts[mask & pft_sum_mask, :, z] /= pft_sum[mask & pft_sum_mask, np.newaxis]

    def recalc_pft_frac_ts(self, rcm_lsm):
        print("normalize to get a sum of 1 and set sea points to missing value")
        for z in range(self.years+1):
            mask = rcm_lsm > 0.0
            pft_sum = np.sum(self.pft_frac_ts[:, :, :, z], axis=2)
            pft_sum_mask = pft_sum > 0.0
            self.pft_frac_ts[mask & pft_sum_mask, :, z] /= pft_sum[mask & pft_sum_mask, np.newaxis]
            self.pft_frac_ts[~mask, :, z] = -999.0

    def recalc_null_pft_frac_ts(self, rcm_lsm):
        print('remove zeros (they still occur due to rounding issues)')
        for z in range(self.years+1):
            mask = rcm_lsm > 0.0
            mask_pft_frac_ts = self.pft_frac_ts[:, :, :, z] < 0.0
            self.pft_frac_ts[mask_pft_frac_ts, z] = 0.0
            pft_sum = np.sum(self.pft_frac_ts[:, :, :, z], axis=2)
            pft_sum_mask = pft_sum > 0.0
            self.pft_frac_ts[mask & pft_sum_mask, :, z] /= pft_sum[mask & pft_sum_mask, np.newaxis]

    def lucas_lut_mcgrath(self, rcm_lsm, mcgrath_frac):
        mask_rcm_lsm = rcm_lsm > 0.0
        for z in range(self.years):
            sum_forest = np.zeros((self.xsize, self.ysize), dtype="float32")
            sum_forest_p1 = np.zeros((self.xsize, self.ysize), dtype="float32")
            sum_mcg = np.zeros((self.xsize, self.ysize), dtype="float32")
            d_mcg_frac = np.zeros((self.xsize, self.ysize, 3), dtype="float32")
            d_frac = np.zeros((self.xsize, self.ysize, 3), dtype="float32")
            d_forest = np.zeros((self.xsize, self.ysize, 3), dtype="float32")
            abssum = np.zeros((self.xsize, self.ysize), dtype="float32")
            zz = self.years - 1 - z
            # compute sun of the three forest types for t and t+1 as well as sum of the forest fractions given by mcgrath
            for k in range(2, 5):
                sum_forest = np.where(mask_rcm_lsm, sum_forest + self.pft_frac_ts[:, :, k, zz], sum_forest)
                sum_forest_p1 = np.where(mask_rcm_lsm, sum_forest_p1 + self.pft_frac_ts[:, :, k, zz+1], sum_forest_p1)
                sum_mcg = np.where(mask_rcm_lsm, sum_mcg + mcgrath_frac[:, :, k-2, zz], sum_mcg)

            mask_sum_forest_p1 = sum_forest_p1 > 0.0
            mask_sum_mcg = sum_mcg > 0.0
            mask_sum_forest = sum_forest > 0.0
            # eliminate missing values
            for k in range(2, 5):
                sum_mcg = np.where(mask_rcm_lsm & ((mcgrath_frac[:, :, k-2, zz] < -998.0) | (mcgrath_frac[:, :, k-2, zz+1] < -998.0)), 0.0, sum_mcg)
            
            # Handle problem condition
            #if sum_forest_p1 < 0.0 or sum_forest < 0.0:
            #    pass
            
            # compute adjustment to mcgrath forest fractions only if sum of forest PFTs and sum of mcgrath data is greater zero
            for k in range(3):
                # conpute difference in forest type fractions (=Delta_mcg) compute absolute sum of Delta_mcg
                d_mcg_frac[:, :, k] = np.where(mask_rcm_lsm & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest, mcgrath_frac[:, :, k, zz+1] - mcgrath_frac[:, :, k, zz], d_mcg_frac[:, :, k])
                abssum = np.where(mask_rcm_lsm & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest, abssum + np.absolute(d_mcg_frac[:, :, k]), abssum)
            
            
            for k in range(2, 5):    
                # compute difference in forest PFTs from t to t+1 (delta_forest_frac_p1=forest_pft_frac/(sum_forest(t+1)*(sum_forest(t+1)-sum_forest(t)))
                filtered_sum_forest_p1 = np.where(mask_rcm_lsm & (sum_forest_p1 > 0.0), sum_forest_p1, 1.0)
                rest = (sum_forest_p1 - sum_forest)
                rest = np.where(mask_rcm_lsm & (sum_forest_p1 > 0.0), rest, 0.0)
                d_forest[:, :, k-2] = np.where(mask_rcm_lsm & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest, self.pft_frac_ts[:, :, k, zz+1] / filtered_sum_forest_p1 * rest, d_forest[:, :, k-2])
            # continue only if there is a change in mcgrath data (i.e. Delta_mcg neq 0)
            mask_abssum = abssum > 0.0
            
            # compute change in fraction due to mcgrath data for each pft Delta_forest_frac_mcg=sum_forest(t+1)*Delta_mcg
            for k in range(2, 5):
                d_frac[:, :, k-2] = np.where(mask_rcm_lsm & mask_abssum & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest, sum_forest_p1 * d_mcg_frac[:, :, k-2], d_frac[:, :, k-2])
            for k in range(2, 5):
                self.pft_frac_ts[:, :, k, zz] = np.where(mask_rcm_lsm & mask_abssum & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest, self.pft_frac_ts[:, :, k, zz+1] - d_frac[:, :, k-2] - d_forest[:, :, k-2], self.pft_frac_ts[:, :, k, zz])
            lhelp = np.ones((self.xsize, self.ysize, 3), dtype="float32")
            helper = np.zeros((self.xsize, self.ysize), dtype="float32")
            # check if negative fractions occured, if so distrute negative fractions to fraction
            for k in range(2, 5):
                helper = np.where(mask_rcm_lsm & mask_abssum & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest & (self.pft_frac_ts[:, :, k, zz] < 0.0), helper + self.pft_frac_ts[:, :, k, zz], helper)
                lhelp[:, :, k-2] = np.where(mask_rcm_lsm & mask_abssum & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest & (self.pft_frac_ts[:, :, k, zz] < 0.0), 0.0, lhelp[:, :, k-2])
                self.pft_frac_ts[:, :, k, zz] = np.where(mask_rcm_lsm & mask_abssum & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest & (self.pft_frac_ts[:, :, k, zz] < 0.0), 0.0, self.pft_frac_ts[:, :, k, zz])
            sum_helper = np.zeros((self.xsize, self.ysize), dtype="float32")
            sum_lhelp = np.zeros((self.xsize, self.ysize), dtype="float32")
            for k in range(3):
                sum_lhelp = np.where(mask_rcm_lsm & mask_abssum & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest, sum_lhelp + lhelp[:, :, k], sum_lhelp)
            
            # if at least one fraction is below zero
            mask_sum_lhelp = sum_lhelp < 3
            # sum of fractions that are larger than zero
            mask_lhelp_1 = lhelp == 1
            for k in range(2, 5):
                sum_helper = np.where(mask_rcm_lsm & mask_abssum & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest & mask_sum_lhelp & mask_lhelp_1[:, :, k-2], sum_helper + self.pft_frac_ts[:, :, k, zz], sum_helper)
            # if there is are fractions which are larger than zero add helper to these fractions
            mask_sum_helper = sum_helper > 0.0
            for k in range(2, 5):
                filtered_sum_helper = np.where(mask_rcm_lsm & mask_abssum & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest & mask_sum_lhelp & mask_sum_helper & mask_lhelp_1[:, :, k-2], sum_helper, 1.0)
                filtered_helper = np.where(mask_rcm_lsm & mask_abssum & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest & mask_sum_lhelp & mask_sum_helper & mask_lhelp_1[:, :, k-2], helper, 0.0)
                self.pft_frac_ts[:, :, k, zz] = np.where(mask_rcm_lsm & mask_abssum & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest & mask_sum_lhelp & mask_sum_helper & mask_lhelp_1[:, :, k-2], self.pft_frac_ts[:, :, k, zz] + (self.pft_frac_ts[:, :, k, zz] / filtered_sum_helper * filtered_helper), self.pft_frac_ts[:, :, k, zz])
                self.pft_frac_ts[:, :, k, zz] = np.where(mask_rcm_lsm & mask_abssum & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest & mask_sum_lhelp & ~mask_sum_helper, self.pft_frac_ts[:, :, k, zz] - ((1 / 3) * helper), self.pft_frac_ts[:, :, k, zz])
            # abssum <= 0.0
            for k in range(2, 5):
                self.pft_frac_ts[:, :, k, zz] = np.where(mask_rcm_lsm & ~mask_abssum & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest, self.pft_frac_ts[:, :, k, zz+1]-d_forest[:, :, k-2], self.pft_frac_ts[:, :, k, zz])
            lhelp = np.ones((self.xsize, self.ysize, 3), dtype="float32")
            helper = np.zeros((self.xsize, self.ysize), dtype="float32")
            for k in range(2, 5):
                helper = np.where(mask_rcm_lsm & ~mask_abssum & (self.pft_frac_ts[:, :, k, zz] < 0.0) & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest, helper + self.pft_frac_ts[:, :, k, zz], helper)
                lhelp[:, :, k-2] = np.where(mask_rcm_lsm & ~mask_abssum & (self.pft_frac_ts[:, :, k, zz] < 0.0) & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest, 0.0, lhelp[:, :, k-2])
                self.pft_frac_ts[:, :, k, zz] = np.where(mask_rcm_lsm & ~mask_abssum & (self.pft_frac_ts[:, :, k, zz] < 0.0) & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest, 0.0, self.pft_frac_ts[:, :, k, zz])
            sum_helper = np.zeros((self.xsize, self.ysize), dtype="float32")
            sum_lhelp = np.zeros((self.xsize, self.ysize), dtype="float32")
            for k in range(3):
                sum_lhelp = np.where(mask_rcm_lsm & ~mask_abssum & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest, sum_lhelp + lhelp[:, :, k], sum_lhelp)
            mask_sum_lhelp = sum_lhelp < 3
            mask_lhelp_1 = lhelp == 1
            for k in range(2, 5):
                sum_helper = np.where(mask_rcm_lsm & ~mask_abssum & mask_sum_lhelp & mask_lhelp_1[:, :, k-2] & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest, sum_helper + self.pft_frac_ts[:, :, k, zz], sum_helper)
            mask_sum_helper = sum_helper > 0.0
            for k in range(2, 5):
                filter_sum_helper = np.where(mask_rcm_lsm & ~mask_abssum & mask_sum_lhelp & mask_sum_helper & mask_lhelp_1[:, :, k-2] & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest, sum_helper, 1.0)
                filtered_helper = np.where(mask_rcm_lsm & ~mask_abssum & mask_sum_lhelp & mask_sum_helper & mask_lhelp_1[:, :, k-2] & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest, helper, 0.0)
                self.pft_frac_ts[:, :, k, zz] = np.where(mask_rcm_lsm & ~mask_abssum & mask_sum_lhelp & mask_sum_helper & mask_lhelp_1[:, :, k-2] & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest, self.pft_frac_ts[:, :, k, zz] + (self.pft_frac_ts[:, :, k, zz] / filtered_sum_helper * filtered_helper), self.pft_frac_ts[:, :, k, zz])
                self.pft_frac_ts[:, :, k, zz] = np.where(mask_rcm_lsm & ~mask_abssum & mask_sum_lhelp & ~mask_sum_helper & mask_sum_forest_p1 & mask_sum_mcg & mask_sum_forest, self.pft_frac_ts[:, :, k, zz] - ((1 / 3) * helper), self.pft_frac_ts[:, :, k, zz])


    def lucas_lut_output(self):
        coords = self.reg.split(",")
        lon = np.linspace(float(coords[0]), float(coords[1]), self.xsize)
        lat = np.linspace(float(coords[2]), float(coords[3]), self.ysize)
        time = np.linspace(self.syear, self.eyear, self.years+1)
        data_array = xr.DataArray(
            self.pft_frac_ts,
            dims=("x", "y", "npft", "time"),
            coords={"x": lon, "y": lat, "time": time}
        )
        # Convert to Dataset and specify the variable name
        dataset = data_array.to_dataset(name="pft_frac")
        # Save the DataArray to a NetCDF file
        dataset.to_netcdf(self.namelist["F_LC_OUT"])

    def generate_namelist(self):
        if self.res == 250:
            ext = "NINT"
        else:
            if self.remap == "bilinear":
                ext = "BIL"
            else:
                ext = ""

        # Select period and self.region
        # CHECK
        if self.region == "Europe":
            if self.res == 100:
                xsize = 1400
                ysize = 630
            elif self.res == 250:
                xsize = 560
                ysize = 252
        elif self.region == "Global":
            xsize = 3600
            ysize = 1800
        elif self.region == "Australasia":
            xsize = 1160
            ysize = 570
        elif self.region == "NorthAmerica":
            xsize = 1900
            ysize = 850
        elif self.region == "GAR011":
            xsize = 145
            ysize = 129
        elif self.region == "Germany":
            if self.res ==  25:
                xsize = 371
                ysize = 351
            elif self.res == 100:
                xsize = 95
                ysize = 91
            elif self.res == 250:
                xsize = 38
                ysize = 36
            elif self.res == 500:
                xsize = 19
                ysize = 18

        if self.scenario == "historical":
            sdir = f"{luhdir}/historic/{self.region}/{self.grid}"
        elif self.scenario == "historical_high":
            sdir = f"{luhdir}/historic_high/{self.region}/{self.grid}"
        elif self.scenario == "historical_low":
            sdir = f"{luhdir}/historic_low/{self.region}/{self.grid}"
        else:
            sdir = f"{luhdir}/scenarios/{self.scenario}/{self.region}/{self.grid}"

        if self.irri:
            ofile = f"{oname}_{self.scenario}_{self.syear}_{self.eyear}_{self.grid}_irri"
        else:
            ofile = f"{oname}_{self.scenario}_{self.syear}_{self.eyear}_{self.grid}"

        if self.mcgback and self.scenario in ["historical", "historical_high", "historical_low"]:
            ofile = f"{ofile}_mcg2"

        if self.addtree and self.scenario not in ["historical", "historical_high", "historical_low"]:
            ofile = f"{ofile}_addtr"
        # Creating directories if they do not exist
        Path(os.path.join(sdir)).mkdir(parents=True, exist_ok=True)
        Path(glcdir).mkdir(parents=True, exist_ok=True)
        Path(pftdir).mkdir(parents=True, exist_ok=True)
        Path(odir).mkdir(parents=True, exist_ok=True)
        namelist_dict = {
            # FILES
            "F_RCM_LSM_IN": f"{pftdir}/{self.glc_lsm}_LSM_{self.grid}.nc", # lsmfile
            "F_LC_IN": f"{pftdir}/PFTS_{self.glc}_{self.grid}_v11.nc", # pftfile
            "F_BACKGRA": f"{pftdir}/GRAB_{self.glc_lsm}_{self.grid}_v11.nc", # grabfile
            "F_BACKSHR": f"{pftdir}/SHRB_{self.glc_lsm}_{self.grid}_v11.nc", # shrbfile
            "F_BACKFOR": f"{pftdir}/FORB_{self.glc_lsm}_{self.grid}_v11.nc", # forbfile
            "F_BACKCRO": f"{pftdir}/CROB_{self.glc_lsm}_{self.grid}_v11.nc", # crobfile
            "F_BACKURB": f"{pftdir}/URBB_{self.glc_lsm}_{self.grid}_v11.nc", # urbbfile
            #"F_MCGRATH": f"{mcgdir}/{mcg}_{self.syear}_{self.eyear}_ForestBckgrdMcGrath_{self.grid}.nc", # mcgfile
            "F_MCGRATH": f"data/{mcg}_{self.syear}_{self.eyear}_ForestBckgrdMcGrath_{self.grid}.nc",
            "F_IRRI_IN": f"{sdir}/irrigation_{self.syear}_{self.eyear}_{self.grid}.nc", # irrfile
            "F_LC_OUT": f"{odir}/{ofile}.nc", # outfile
            "F_FOR2CRO": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_for2cro_{ext}_{self.grid}.nc", # for2cro
            "F_CRO2FOR": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_cro2for_{ext}_{self.grid}.nc", # cro2for
            "F_FOR2RAN": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_for2ran_{ext}_{self.grid}.nc", # for2ran
            "F_RAN2FOR": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_ran2for_{ext}_{self.grid}.nc", # ran2for
            "F_FOR2PAS": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_for2pas_{ext}_{self.grid}.nc", # for2pas
            "F_PAS2FOR": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_pas2for_{ext}_{self.grid}.nc", # pas2for
            "F_FOR2URB": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_for2urb_{ext}_{self.grid}.nc", # for2urb
            "F_NFV2CRO": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_nfv2cro_{ext}_{self.grid}.nc", # nfv2cro
            "F_CRO2NFV": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_cro2nfv_{ext}_{self.grid}.nc", # cro2nfv
            "F_NFV2RAN": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_nfv2ran_{ext}_{self.grid}.nc", # nfv2ran
            "F_RAN2NFV": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_ran2nfv_{ext}_{self.grid}.nc", # ran2nfv
            "F_NFV2URB": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_nfv2urb_{ext}_{self.grid}.nc", # nfv2urb
            "F_RAN2CRO": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_ran2cro_{ext}_{self.grid}.nc", # ran2cro
            "F_CRO2RAN": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_cro2ran_{ext}_{self.grid}.nc", # cro2ran
            "F_RAN2URB": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_ran2urb_{ext}_{self.grid}.nc", # ran2urb
            "F_PAS2CRO": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_pas2cro_{ext}_{self.grid}.nc", # pas2cro
            "F_CRO2PAS": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_cro2pas_{ext}_{self.grid}.nc", # cro2pas
            "F_PAS2URB": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_pas2urb_{ext}_{self.grid}.nc", # pas2urb
            "F_CRO2URB": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_cro2urb_{ext}_{self.grid}.nc", # cro2urb
            "F_NFV2PAS": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_nfv2pas_{ext}_{self.grid}.nc", # nfv2pas
            "F_PAS2NFV": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_pas2nfv_{ext}_{self.grid}.nc", # pas2nfv
            "F_RAN2PAS": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_ran2pas_{ext}_{self.grid}.nc", # ran2pas
            "F_FOR2NFV": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_for2nfv_{ext}_{self.grid}.nc", # for2nfv
            "F_NFV2FOR": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_nfv2for_{ext}_{self.grid}.nc", # nfv2for
            "F_URB2FOR": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_urb2for_{ext}_{self.grid}.nc", # urb2for
            "F_URB2NFV": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_urb2nfv_{ext}_{self.grid}.nc", # urb2nfv
            "F_URB2CRO": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_urb2cro_{ext}_{self.grid}.nc", # urb2cro
            "F_URB2PAS": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_urb2pas_{ext}_{self.grid}.nc", # urb2pas
            "F_URB2RAN": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_urb2ran_{ext}_{self.grid}.nc", # urb2ran
            "F_ADDTREE": f"{sdir}/addtree_{self.syear}_{self.eyear}_{self.grid}.nc", # nat2for
            "F_GRID": f"{scriptsdir}/grid_{self.grid}", # grid
            "XSIZE": xsize,
            "YSIZE": ysize,
        }
        return namelist_dict

    def func_prepare_luh2_data(self):
        if self.scenario == "historical":
            sdir="historic"
            sfile="states"
            tfile="transitions"
            mfile="management"
        elif self.scenario == "historical_high":
            sdir="historic_high"
            sfile="states"
            tfile="transitions"
            mfile="management"
        elif self.scenario == "historical_low":
            sdir="historic_low"
            sfile="states"
            tfile="transitions"
            mfile="management"
        elif self.scenario in scenario_dict.keys():
            sdir=f"scenarios/{self.scenario}"
            afile=f"added_tree_cover_input4MIPs_landState_ScenarioMIP_UofMD-{scenario_dict[self.scenario]}-2-1-f_gn_2015-2100"
            sfile=f"multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-{scenario_dict[self.scenario]}-2-1-f_gn_2015-2100"
            tfile=f"multiple-transitions_input4MIPs_landState_ScenarioMIP_UofMD-{scenario_dict[self.scenario]}-2-1-f_gn_2015-2100"
            mfile=f"multiple-management_input4MIPs_landState_ScenarioMIP_UofMD-{scenario_dict[self.scenario]}-2-1-f_gn_2015-2100"

        # interpolation
        if self.grid == "reg025_Europe":
            ext="NINT"
            remap_com="invertlat"
            cutting=''
        else:
            if self.remap == "bilinear":
                ext="BIL"
                remap_com=f"remapbil"
                if self.grid == "EUR-011":
                    cutting = "-selindexbox,2,434,2,434"
                else:
                    cutting = ""
            elif self.remap == "con2":
                ext = "CON2"
                remap_com = f"remapcon2"
                if self.grid == 'EUR-011':
                    cutting = "-selindexbox,2,434,2,434"
                else:
                    cutting = ''

        path_region = os.path.join(luhdir, sdir, self.region)
        path_sdir = os.path.join(luhdir, sdir)
        if self.trans:
            print_section_heading(f"Selecting variables for transitions")
            ofile=f"{tfile}_{self.syear}_{self.eyear}_{self.region}.nc"
            cdo.selyear(f"{self.syear}/{self.eyear}", input=f"{datadir}/{tfile}.nc", output=f"{path_sdir}/{self.region}/tmp_{ofile}")
            cdo.sellonlatbox(self.reg, input=f"-selvar,{vars_trans} {path_sdir}/{self.region}/tmp_{ofile}", output=f"{path_sdir}/{self.region}/{ofile}")
            fromto_array = [
                {"varn": "for2urb", "for_1": FOR, "for_2": URB, "outvar_condition": None},
                {"varn": "urb2for", "for_1": URB, "for_2": FOR, "outvar_condition": "primf"},
                {"varn": "for2nfv", "for_1": FOR, "for_2": NFV, "outvar_condition": "primn"},
                {"varn": "for2cro", "for_1": FOR, "for_2": CRO, "outvar_condition": None},
                {"varn": "nfv2for", "for_1": NFV, "for_2": FOR, "outvar_condition": "primf"},
                {"varn": "cro2for", "for_1": CRO, "for_2": FOR, "outvar_condition": "primf"},
                {"varn": "cro2urb", "for_1": CRO, "for_2": URB, "outvar_condition": None},
                {"varn": "urb2cro", "for_1": URB, "for_2": CRO, "outvar_condition": None},
                {"varn": "cro2nfv", "for_1": CRO, "for_2": NFV, "outvar_condition": "primn"},
                {"varn": "nfv2cro", "for_1": NFV, "for_2": CRO, "outvar_condition": None},
                {"varn": "nfv2urb", "for_1": NFV, "for_2": URB, "outvar_condition": None},
                {"varn": "urb2nfv", "for_1": URB, "for_2": NFV, "outvar_condition": "primn"},
                {"varn": "ran2nfv", "for_1": RAN, "for_2": NFV, "outvar_condition": "primn"},
                {"varn": "nfv2ran", "for_1": NFV, "for_2": RAN, "outvar_condition": None},
                {"varn": "ran2for", "for_1": RAN, "for_2": FOR, "outvar_condition": "primf"},
                {"varn": "for2ran", "for_1": FOR, "for_2": RAN, "outvar_condition": None},
                {"varn": "ran2cro", "for_1": RAN, "for_2": CRO, "outvar_condition": None},
                {"varn": "cro2ran", "for_1": CRO, "for_2": RAN, "outvar_condition": None},
                {"varn": "ran2urb", "for_1": RAN, "for_2": URB, "outvar_condition": None},
                {"varn": "urb2ran", "for_1": URB, "for_2": RAN, "outvar_condition": None},
                {"varn": "pas2nfv", "for_1": PAS, "for_2": NFV, "outvar_condition": "primn"},
                {"varn": "nfv2pas", "for_1": NFV, "for_2": PAS, "outvar_condition": None},
                {"varn": "pas2for", "for_1": PAS, "for_2": FOR, "outvar_condition": "primf"},
                {"varn": "for2pas", "for_1": FOR, "for_2": PAS, "outvar_condition": None},
                {"varn": "pas2cro", "for_1": PAS, "for_2": CRO, "outvar_condition": None},
                {"varn": "cro2pas", "for_1": CRO, "for_2": PAS, "outvar_condition": None},
                {"varn": "pas2urb", "for_1": PAS, "for_2": URB, "outvar_condition": None},
                {"varn": "ran2pas", "for_1": RAN, "for_2": PAS, "outvar_condition": None},
                {"varn": "urb2pas", "for_1": URB, "for_2": PAS, "outvar_condition": None}

            ]
            for data in fromto_array:
                self.fromto(data["varn"], data["for_1"], data["for_2"], tfile, ext, cutting, path_region, remap_com, data["outvar_condition"])

        if self.state:
            print_section_heading(f"Selecting variables for states")
            cdo.selyear(f"{self.syear}/{self.eyear}", input=f"{datadir}/{sfile}.nc", output=f"{path_sdir}/{self.region}/tmp_{sfile}_{self.syear}_{self.eyear}.nc")
            cdo.sellonlatbox(self.reg, input=f"-selvar,{vars_state} {path_sdir}/{self.region}/tmp_{sfile}_{self.syear}_{self.eyear}.nc", output=f"{path_region}/states_{self.syear}_{self.eyear}_{self.region}.nc")
            if remap_com == "invertlat":
                cdo.invertlat(input=f"{path_region}/states_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/states_{self.syear}_{self.eyear}_{self.grid}.nc")
            elif remap_com == "remapbil":
                cdo.remapbil(f"{scriptsdir}/grid_{self.grid}", input=f"{path_region}/states_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/states_{self.syear}_{self.eyear}_{self.grid}.nc")
            elif remap_com == "remapcon2":
                cdo.remapcon2(f"{scriptsdir}/grid_{self.grid}", input=f"{path_region}/states_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/states_{self.syear}_{self.eyear}_{self.grid}.nc")
        if self.addtree:
            print_section_heading(f"Selecting variables for added tree cover")
            cdo.sellonlatbox(self.reg, input=f"-selyear,{self.syear}/{self.eyear} -selvar,added_tree_cover {path_sdir}/{afile}.nc", output=f"{path_region}/addtree_{self.syear}_{self.eyear}_{self.region}.nc")
            if remap_com == "invertlat":
                cdo.invertlat(input=f"{path_region}/addtree_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/{self.grid}/addtree_{self.syear}_{self.eyear}_{self.grid}.nc")
            elif remap_com == "remapbil":
                cdo.remapbil(f"{scriptsdir}/grid_{self.grid}", input=f"{path_region}/addtree_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/{self.grid}/addtree_{self.syear}_{self.eyear}_{self.grid}.nc")
            elif remap_com == "remapcon2":
                cdo.remapcon2(f"{scriptsdir}/grid_{self.grid}", input=f"{path_region}/addtree_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/{self.grid}/addtree_{self.syear}_{self.eyear}_{self.grid}.nc")
            cdo.copy(input=f'{path_region}/{self.grid}/addtree_{self.syear}_{self.eyear}_{self.grid}.nc', output=f"{path_region}/{self.grid}/addtree_{self.syear}_{self.eyear}_{self.grid}.nc")
        # compute irragtion fraction 
        if self.irri:
            print_section_heading(f"Selecting variables for irrigation")
            cdo.sellonlatbox(self.reg, input=f"-selyear,{self.syear}/{self.eyear} -selvar,{vars_irrig} {datadir}/{mfile}.nc", output=f"{path_region}/irri_vars_{self.syear}_{self.eyear}_{self.region}.nc")

            if self.scenario in ["historical", "historical_low", "historical_high"]:
                if remap_com == "invertlat":
                    cdo.invertlat(input=f"-chname,{IRR[0]},irrig_frac -selvar,{IRR[0]} {path_region}/irri_vars_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/{self.grid}/irrigation_{self.syear}_{self.eyear}_{self.grid}.nc")
                elif remap_com == "remapbil":
                    cdo.remapbil(f"{scriptsdir}/grid_{self.grid}", input=f"-chname,{IRR[0]},irrig_frac -selvar,{IRR[0]} {path_region}/irri_vars_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/{self.grid}/irrigation_{self.syear}_{self.eyear}_{self.grid}.nc")
                elif remap_com == "remapcon2":
                    cdo.remapcon2(f"{scriptsdir}/grid_{self.grid}", input=f"-chname,{IRR[0]},irrig_frac -selvar,{IRR[0]} {path_region}/irri_vars_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/{self.grid}/irrigation_{self.syear}_{self.eyear}_{self.grid}.nc")
                cdo.copy(options="-setmisstoc,-999", input=f'{path_region}/{self.grid}/irrigation_{self.syear}_{self.eyear}_{self.grid}.nc', output=f"{path_region}/{self.grid}/irrigation_{self.syear}_{self.eyear}_$grid.nc")
            else:
                if remap_com in ["remapbil", "remapcon2"]:
                    cdo.setmisstoc(input=f"-999 -{remap_com},{scriptsdir}/grid_{self.grid} -varssum -selvar,{vars_crops} {path_region}/states_{self.syear}_{self.eyear}_{self.grid}.nc", output=f"{path_region}/sum_crop_frac.nc")
                else:
                    cdo.setmisstoc(input=f"-999 -{remap_com} -varssum -selvar,{vars_crops} {path_region}/states_{self.syear}_{self.eyear}_{self.grid}.nc", output=f"{path_region}/sum_crop_frac.nc")

                # Change variables names

                for n in range(5):
                    if remap_com in ["remapbil", "remapcon2"]:
                        cdo.mul(input=f"-{remap_com}, grid_{self.grid} -selvar,{IRR[n]} {path_region}/irri_vars_{self.syear}_{self.eyear}_{self.region}.nc -{remap_com}, grid_{self.grid} -selvar,{ICR[n]} {path_region}/states_{self.syear}_{self.eyear}_{self.grid}.nc", output=f"{sdir}/{self.region}/dummy.nc")
                    else:
                        cdo.mul(input=f"-{remap_com} -selvar,{IRR[n]} {path_region}/irri_vars_{self.syear}_{self.eyear}_{self.region}.nc -{remap_com} -selvar,{ICR[n]} {path_region}/states_{self.syear}_{self.eyear}_{self.grid}.nc", output=f"{sdir}/{self.region}/dummy.nc")
                    if n == 0:
                        cdo.chname(input=f"{IRR[0]},irrig_frac -selvar,{IRR[0]} {path_region}/dummy.nc", output=f"{path_region}/sum_irri_frac_{self.syear}_{self.eyear}_{self.grid}.nc")
                    else:
                        cdo.add(input=f"-chname,{IRR[n]},irrig_frac -selvar,{IRR[n]} {path_region}/dummy.nc {path_region}/sum_irri_frac_{self.syear}_{self.eyear}_{self.grid}.nc")   
                        shutil.move(f"{path_region}/dummy2.nc", f"{path_region}/sum_irri_frac_{self.syear}_{self.eyear}_{self.grid}.nc")
                    os.remove(f"{path_region}/dummy.nc")
                cdo.div(input=f"{path_region}/sum_irri_frac_{self.syear}_{self.eyear}_{self.grid}.nc {path_region}/sum_crop_frac.nc", output=f"{path_region}/{self.grid}/irrigation_{self.syear}_{self.eyear}_{self.grid}.nc")
                cdo.copy(input=f"-setmisstoc,-999 {path_region}/{self.grid}/irrigation_{self.syear}_{self.eyear}_{self.grid}.nc", output=f"{path_region}/{self.grid}/irrigation_{self.syear}_{self.eyear}_{self.grid}.nc")

    def func_prepare_mcgrath(self):
        # commands for interpolation to given grid
        if self.remap == "bilinear":
            ext = "BIL"
            remap_com = f"remapbil, grid_{self.grid}"
            if self.grid == "EUR-011":
                cutting = "-selindexbox,2,434,2,434"
            else:
                cutting = ""
        else:
            remap_com = ""

        ifile = f"{datadir}/{self.mcg}_{self.syear}_{self.mcgrath_eyear}.nc"

        # compute background for LUT classes using zonal mean
        cdo.chname(f"maxvegetfrac,{PFT_TeBrEv}", input=f"-vertsum -sellevel,{TeBrEv} {ifile}", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_TeBrEv.nc")
        cdo.chname(f"maxvegetfrac,{PFT_TeBrDec}", input=f"-vertsum -sellevel,{TeBrDec} {ifile}", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_TeBrDec.nc")
        cdo.chname(f"maxvegetfrac,{PFT_ConEv}", input=f"-vertsum -sellevel,{EvCon} {ifile}", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_ConEv.nc")
        cdo.chname(f"maxvegetfrac,forest", input=f"-vertsum -sellevel,{TeBrEv},{TeBrDec},{EvCon} {ifile}", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_FOR.nc")
        cdo.merge(input=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_TeBrEv.nc {glcdir}/{self.lcd}_{self.syear}_{self.eyear}_TeBrDec.nc {glcdir}/{self.lcd}_{self.syear}_{self.eyear}_ConEv.nc", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_dummy.nc")
        if self.reg:
            cdo.setmisstoc(-999, input=f"-remapbil,{scriptsdir}/grid_{self.grid} -sellonlatbox,{self.reg} -div {glcdir}/{self.lcd}_{self.syear}_{self.eyear}_dummy.nc -varssum {glcdir}/{self.lcd}_{self.syear}_{self.eyear}_dummy.nc", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_ForestBckgrdMcGrath_{self.grid}.nc")
        else:
            cdo.setmisstoc(-999, input=f"-remapbil,{scriptsdir}/grid_{self.grid} -div {glcdir}/{self.lcd}_{self.syear}_{self.eyear}_dummy.nc -varssum {glcdir}/{self.lcd}_{self.syear}_{self.eyear}_dummy.nc", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_ForestBckgrdMcGrath_{self.grid}.nc")
        if self.mcgrath_eyear:
            if self.mcgrath_eyear < self.eyear:
                for year in range(self.mcgrath_eyear, self.eyear+1):
                    cdo.setdate(f"{year}-06-15", input=f"-selyear,2010 {glcdir}/{self.lcd}_{self.syear}_{self.eyear}_ForestBckgrdMcGrath_{self.grid}.nc", output=f"{glcdir}/dummy_{year}.nc")
        cdo.mergetime(input=f"{glcdir}/dummy_????.nc", output=f"{glcdir}/{self.lcd}_{self.mcgrath_eyear}_{self.eyear}_ForestBckgrdMcGrath_{self.grid}.nc")
        cdo.mergetime(input=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_ForestBckgrdMcGrath_{self.grid}.nc {glcdir}/{self.lcd}_{self.mcgrath_eyear}_{self.eyear}_ForestBckgrdMcGrath_{self.grid}.nc", output=f"{glcdir}/{self.lcd}_{self.syear}_2015_ForestBckgrdMcGrath_{self.grid}.nc")

    def fromto(self, varn, for_1, for_2, tfile, ext, cutting, path, remap_com, outvar_condition=None):
        odir = self.grid
        print(f"{for_1} to {for_2}")
        # combine land-use changes using reclassifcation
        ifile=f"{tfile}_{self.syear}_{self.eyear}_{self.region}.nc"
        logfile=f"{tfile}_{self.syear}_{self.eyear}_{self.region}.log"
        ofile=f"transitions_{self.syear}_{self.eyear}_{self.region}_{varn}"
        cdo.mulc(0, input=f"-selvar,primf_to_urban {path}/{ifile}", output=f"{path}/dummy.nc")
        for inivar in for_1:
            for outvar in for_2:
                if not outvar_condition or outvar != outvar_condition:
                    cdo.add(input=f"-selvar,{inivar}_to_{outvar} {path}/{ifile} {path}/dummy.nc", output=f"{path}/{ofile}.nc")
                    cdo.chname(f"{inivar}_to_{outvar},{varn}", input=f"{path}/{ofile}.nc", output=f"{path}/dummy.nc")
        shutil.move(f"{path}/dummy.nc", f"{path}/{ofile}.nc")
        if remap_com == "invertlat":
            cdo.invertlat(input=f"{path}/{ofile}.nc", output=f"{path}/{self.grid}/{ofile}_{ext}_{self.grid}_2.nc")
        elif remap_com == "remapbil":
            cdo.remapbil(f"{scriptsdir}/grid_{self.grid}", input=f"{path}/{ofile}.nc", output=f"{path}/{self.grid}/{ofile}_{ext}_{self.grid}_2.nc")
        elif remap_com == "remapcon2":
            cdo.remapcon2(f"{scriptsdir}/grid_{self.grid}", input=f"{path}/{ofile}.nc", output=f"{path}/{self.grid}/{ofile}_{ext}_{self.grid}_2.nc")
        cdo.copy(input=f"{cutting} -setmisstoc,-999. {path}/{self.grid}/{ofile}_{ext}_{self.grid}_2.nc", output=f"{path}/{self.grid}/{ofile}_{ext}_{self.grid}.nc", options="-f nc")

    def plot_pft_frac_ts(self):
        year = self.plot_year
        Path(os.path.join(plotdir, self.region, str(self.syear+year))).mkdir(parents=True, exist_ok=True)
        data = self.lucas_lut_help()
        if self.plot_npft:
            self.main_plot(self.plot_npft, year)
            self.diff_plot(self.plot_npft, year)
            self.diff_plot(self.plot_npft, year, data=data, directory=f"lucas_diff")
            self.diff_2_plot(self.plot_npft, year, self.pft_frac_ts, data, directory=f"lucas_python_fortran_diff")
        else:
            for inpft in range(self.npfts):
                self.main_plot(inpft+1, year)
                self.diff_plot(inpft+1, year)
                self.diff_plot(inpft+1, year, data=data, directory=f"lucas_diff")
                self.diff_2_plot(inpft+1, year, self.pft_frac_ts, data, directory=f"lucas_python_fortran_diff")

    def diff_2_plot(self, npft, year, data_1, data_2, directory):
        Path(os.path.join(plotdir, self.region, directory)).mkdir(parents=True, exist_ok=True)
        coords = self.reg.split(",")
        lon = self.pft_frac.lon.values
        lat = self.pft_frac.lat.values
        # Create a meshgrid for longitude and latitude
        lon_grid, lat_grid = np.meshgrid(lon, lat)
        # Plotting
        fig = plt.figure(figsize=(12, 10))
        ax = plt.axes(projection=ccrs.PlateCarree())

        # Add natural earth features for better visualization
        ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=1, edgecolor='black')
        ax.coastlines(resolution='10m')
        ax.add_feature(cfeature.LAND, facecolor="white")
        ax.add_feature(cfeature.LAKES, facecolor='blue')
        ax.add_feature(cfeature.RIVERS, edgecolor='blue')
        # Plot the data
        levels = np.linspace(-0.4, 0.4, 100)  # 30 intervals
        #values = (data_1[:, :, npft-1, year].T-data_1[:, :, npft-1, 0].T) - (data_2[:, :, npft-1, year].T-data_2[:, :, npft-1, 0].T)
        values = (data_1[:, :, npft-1, year].T) - (data_2[:, :, npft-1, year].T)
        # Create a custom colormap with white at the center
        colors = [(0, 0, 1), (1, 1, 1), (1, 0, 0)]  # Blue -> White -> Red
        n_bins = 100  # Discretize the color range into 100 bins
        cmap_name = 'blue_white_red'
        cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

        # Create the filled contour plot
        contour = ax.contourf(lon_grid, lat_grid, values, levels=levels, cmap=cm, transform=ccrs.PlateCarree(), alpha=1)

        # Add colorbar with label
        cbar = plt.colorbar(contour, orientation='vertical', pad=0.22, aspect=50)
        cbar.set_label('PFT Fraction')

        # Set the extent to focus on Germany (coordinates provided by self.reg)
        ax.set_extent([float(coords[0]), float(coords[1]), float(coords[2]), float(coords[3])], crs=ccrs.PlateCarree())

        # Add gridlines with labels
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False

        # Set the title of the plot
        plt.title(f'Average {npft} PFT Fraction Diff (PYTHON-FORTRAN) over Germany for {self.syear + year}')

        # Save the plot to the specified directory
        plt.savefig(f'{os.path.join(plotdir, self.region, directory)}/pft_frac_{self.region}_npft_{npft}_diff_{self.syear + year}.png')

        # Close the plot to free memory
        plt.close()

    def diff_plot(self, npft, year, data=None, directory=None):
        directory = f"diff_{self.syear+year}" if directory is None else directory
        data = self.pft_frac_ts if data is None else data
        Path(os.path.join(plotdir, self.region, directory)).mkdir(parents=True, exist_ok=True)
        coords = self.reg.split(",")
        lon = self.pft_frac.lon.values
        lat = self.pft_frac.lat.values
        # Create a meshgrid for longitude and latitude
        lon_grid, lat_grid = np.meshgrid(lon, lat)
        # Plotting
        fig = plt.figure(figsize=(12, 10))
        ax = plt.axes(projection=ccrs.PlateCarree())

        # Add natural earth features for better visualization
        ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=1, edgecolor='black')
        ax.coastlines(resolution='10m')
        ax.add_feature(cfeature.LAND, facecolor="white")
        ax.add_feature(cfeature.LAKES, facecolor='blue')
        ax.add_feature(cfeature.RIVERS, edgecolor='blue')
        # Plot the data
        levels = np.linspace(-0.4, 0.3, 30)  # 30 intervals
        contour = ax.contourf(lon_grid, lat_grid, data[:, :, npft-1, self.years].T-data[:, :, npft-1, 0].T, levels=levels, transform=ccrs.PlateCarree(), alpha=1)
        # Add colorbar
        cbar = plt.colorbar(contour, orientation='vertical', pad=0.22, aspect=50)
        cbar.set_label('PFT Fraction')
        # Set the extent to focus on Germany
        ax.set_extent([float(coords[0]), float(coords[1]), float(coords[2]), float(coords[3])], crs=ccrs.PlateCarree())
        # Add gridlines
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
        # Title
        plt.title(f'Average {npft} PFT Fraction Diff over Germany for {self.syear}-{self.eyear}')
        # Save the plot
        plt.savefig(f'{os.path.join(plotdir, self.region, directory)}/pft_frac_{self.region}_npft_{npft}_diff_{self.syear}_{self.eyear}.png')
        plt.close()

    def main_plot(self, npft, year):
        coords = self.reg.split(",")
        lon = self.pft_frac.lon.values
        lat = self.pft_frac.lat.values
        # Create a meshgrid for longitude and latitude
        lon_grid, lat_grid = np.meshgrid(lon, lat)
        # Plotting
        fig = plt.figure(figsize=(12, 10))
        ax = plt.axes(projection=ccrs.PlateCarree())

        # Add natural earth features for better visualization
        ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=1, edgecolor='black')
        ax.coastlines(resolution='10m')
        #ax.add_feature(cfeature.LAND, facecolor="white")
        ax.add_feature(cfeature.LAKES, facecolor='blue')
        ax.add_feature(cfeature.RIVERS, edgecolor='blue')
        # Plot the data
        levels = np.linspace(0.0, 1.0, 30)  # 30 intervals
        contour = ax.contourf(lon_grid, lat_grid, self.pft_frac_ts[:, :, npft-1, year].T, levels=levels, transform=ccrs.PlateCarree(), alpha=1)
        # Add colorbar
        cbar = plt.colorbar(contour, orientation='vertical', pad=0.22, aspect=50)
        cbar.set_label('PFT Fraction')
        # Set the extent to focus on Germany
        ax.set_extent([float(coords[0]), float(coords[1]), float(coords[2]), float(coords[3])], crs=ccrs.PlateCarree())
        # Add gridlines
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
        # Title
        plt.title(f'Average {npft} PFT Fraction over Germany for {self.syear+year}')
        # Save the plot
        plt.savefig(f'{os.path.join(plotdir, self.region, str(self.syear+year))}/pft_frac_{self.region}_npft_{npft}_{self.syear+year}.png')
        plt.close()

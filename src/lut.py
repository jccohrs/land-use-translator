from pathlib import Path
import shutil
from cdo import *
from src.config import *
import os
import xarray as xr
import numpy as np
import netCDF4
from src.utils import print_section_heading

import matplotlib.pyplot as plt
from scipy.stats import linregress
from matplotlib.patches import Patch
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import ListedColormap, BoundaryNorm

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
        self.pft_grass_default = self.namelist["GRADEF"]
        self.pft_crops_default = self.namelist["CRODEF"]
        self.pft_shrubs_default = self.namelist["SHRDEF"]
        self.pft_forest_default = self.namelist["FORDEF"]
        self.pft_urban_default = self.namelist["URBDEF"]
        self.pft_crops_default = self.namelist["CRODEF"]
        self.xsize = self.namelist["XSIZE"]
        self.ysize = self.namelist["YSIZE"]
        self.pft_frac_ts = np.zeros((self.xsize, self.ysize, self.npfts, self.years+1))
        self.pfts_shrubs_grass = self.pfts_grass + self.pfts_shrubs

    def lucas_lut_forward(self):
        if self.addtree:
            rcm_lsm, rcm_lsm_help, pft_help, grass_backgr_help, crops_backgr_help, forest_backgr_help, shrubs_backgr_help, \
            shrubs_grass_backgr, nfv2cro, cro2nfv, for2cro, cro2for, ran2cro, cro2ran, pas2cro, cro2pas, cro2urb, nfv2urb, \
            for2urb, ran2urb, pas2urb, for2pas, pas2for, nfv2pas, ran2pas, pas2nfv, for2ran, ran2for, for2nfv, nfv2for, mcgrath_frac_help, \
            urb2cro, urb2nfv, urb2for, urb2ran, urb2pas, nfv2ran, nat2for, urban_backgr_help = self.lucas_lut_input()
        else:
            rcm_lsm, rcm_lsm_help, pft_help, grass_backgr_help, crops_backgr_help, forest_backgr_help, shrubs_backgr_help, \
            shrubs_grass_backgr, nfv2cro, cro2nfv, for2cro, cro2for, ran2cro, cro2ran, pas2cro, cro2pas, cro2urb, nfv2urb, \
            for2urb, ran2urb, pas2urb, for2pas, pas2for, nfv2pas, ran2pas, pas2nfv, for2ran, ran2for, for2nfv, nfv2for, mcgrath_frac_help, \
            urb2cro, urb2nfv, urb2for, urb2ran, urb2pas, nfv2ran, urban_backgr_help = self.lucas_lut_input()
        for z in range(self.years):
                    
            pft_help = self.lucas_lut_transrules(for2cro[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_crops, self.pfts_forest, self.pfts_shrubs, self.pfts_grass, self.nr_crops, self.nr_forest, self.nr_shrubs, self.nr_grass, self.pft_crops_default, crops_backgr_help, False, 3, False)
            pft_help = self.lucas_lut_transrules(nfv2cro[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_crops, self.pfts_shrubs, self.pfts_grass, 0, self.nr_crops, self.nr_shrubs, self.nr_grass, 1, self.pft_crops_default, crops_backgr_help, False, 2, False)
            pft_help = self.lucas_lut_transrules(ran2cro[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_crops, self.pfts_shrubs, self.pfts_grass, 0, self.nr_crops, self.nr_shrubs, self.nr_grass, 1, self.pft_crops_default, crops_backgr_help, False, 2, False)
            pft_help = self.lucas_lut_transrules(pas2cro[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_crops, self.pfts_grass, 0, 0, self.nr_crops, self.nr_grass, 1, 1, self.pft_crops_default, crops_backgr_help, False, 1, False)
            pft_help = self.lucas_lut_transrules(cro2for[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_forest, self.pfts_crops, 0, 0, self.nr_forest, self.nr_crops, 1, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(cro2nfv[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_shrubs_grass, self.pfts_crops, 0, 0, self.nr_shrubs + self.nr_grass, self.nr_crops, 1, 1, self.pft_shrubs_default, shrubs_grass_backgr, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(cro2ran[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_grass, self.pfts_crops, 0, 0, self.nr_grass, self.nr_crops, 1, 1, self.pft_grass_default, grass_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(cro2pas[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_grass, self.pfts_crops, 0, 0, self.nr_grass, self.nr_crops, 1, 1, self.pft_grass_default, grass_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(cro2urb[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_urban, self.pfts_crops, 0, 0, self.nr_urban, self.nr_crops, 1, 1, self.pft_urban_default, urban_backgr_help, False, 1, False)
            pft_help = self.lucas_lut_transrules(for2urb[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_urban, self.pfts_forest, self.pfts_shrubs, self.pfts_grass, self.nr_urban, self.nr_forest, self.nr_shrubs, self.nr_grass, self.pft_urban_default, urban_backgr_help, False, 3, False)
            pft_help = self.lucas_lut_transrules(nfv2urb[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_urban, self.pfts_shrubs, self.pfts_grass, 0, self.nr_urban, self.nr_shrubs, self.nr_grass, 1, self.pft_urban_default, self.backgrd, 2, False)
            pft_help = self.lucas_lut_transrules(ran2urb[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_urban, self.pfts_shrubs_grass, 0, 0, self.nr_urban, self.nr_shrubs + self.nr_grass, 1, 1, self.pft_urban_default, urban_backgr_help, False, 1, False)
            pft_help = self.lucas_lut_transrules(pas2urb[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_urban, self.pfts_grass, 0, 0, self.nr_urban, self.nr_grass, 1, 1, self.pft_urban_default, urban_backgr_help, False, 1, False)
            pft_help = self.lucas_lut_transrules(urb2cro[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_crops, self.pfts_urban, 0, 0, self.nr_crops, self.nr_urban, 1, 1, self.pft_crops_default, crops_backgr_help, False, 1, False)
            pft_help = self.lucas_lut_transrules(urb2nfv[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_shrubs_grass, self.pfts_urban, 0, 0, self.nr_shrubs + self.nr_grass, self.nr_urban, 1, 1, self.pft_shrubs_default, shrubs_grass_backgr, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(urb2for[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_forest, self.pfts_urban, 0, 0, self.nr_forest, self.nr_urban, 1, 1, self.pft_forest_default, forest_backgr, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(urb2ran[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_shrubs_grass, self.pfts_urban, 0, 0, self.nr_shrubs + self.nr_grass, self.nr_urban, 1, 1, self.pft_shrubs_default, shrubs_grass_backgr, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(urb2pas[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_grass, self.pfts_urban, 0, 0, self.nr_grass, self.nr_urban, 1, 1, pft_grass_default, grass_backgr[i, j, :], self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(for2pas[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_grass, self.pfts_forest, 0, 0, self.nr_grass, self.nr_forest, 1, 1, self.pft_grass_default, grass_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(nfv2pas[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_grass, self.pfts_shrubs, 0, 0, self.nr_grass, self.nr_shrubs, 1, 1, self.pft_grass_default, grass_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(ran2pas[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_grass, self.pfts_shrubs, 0, 0, self.nr_grass, self.nr_shrubs, 1, 1, self.pft_grass_default, grass_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(pas2for[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_forest, self.pfts_grass, 0, 0, self.nr_forest, self.nr_grass, 1, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(pas2nfv[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_shrubs, self.pfts_grass, 0, 0, self.nr_shrubs, self.nr_grass, 1, 1, self.pft_shrubs_default, shrubs_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(for2ran[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_shrubs_grass, self.pfts_forest, 0, 0, self.nr_shrubs + self.nr_grass, self.nr_forest, 1, 1, self.pft_shrubs_default, shrubs_grass_backgr, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(nfv2ran[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_grass, self.pfts_shrubs, 0, 0, self.nr_grass, self.nr_shrubs, 1, 1, self.pft_grass_default, grass_backgr, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(ran2for[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_forest, self.pfts_shrubs, self.pfts_grass, 0, self.nr_forest, self.nr_shrubs, self.nr_grass, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 2, False)
            pft_help = self.lucas_lut_transrules(for2nfv[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_shrubs_grass, self.pfts_forest, 0, 0, self.nr_shrubs + self.nr_grass, self.nr_forest, 1, 1, self.pft_shrubs_default, shrubs_grass_backgr, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(nfv2for[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_forest, self.pfts_shrubs, self.pfts_grass, 0, self.nr_forest, self.nr_shrubs, self.nr_grass, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 2, False)
            if self.addtree:
                self.lucas_lut_transrules(nat2for[z, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_forest, self.pfts_shrubs, self.pfts_grass, 0, self.nr_forest, self.nr_shrubs, self.nr_grass, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 2, False)
            self.pft_frac_ts[:, :, :, z] = self.pft_help[:, :, :]
        print('land use change finished')
        # NORMALIZE TO GET A SUM OF 1 AND SET SEA POINTS TO MISSING VALUE
        self.recalc_pft_frac_ts(rcm_lsm)
        if self.irri:
            self.lucas_lut_irrigation(rcm_lsm)
        if self.mcgrath:
            self.lucas_lut_mcgrath(rcm_lsm)

    def lucas_lut_backward(self):
        rcm_lsm, rcm_lsm_help, pft_help, grass_backgr_help, crops_backgr_help, forest_backgr_help, shrubs_backgr_help, \
        shrubs_grass_backgr, nfv2cro, cro2nfv, for2cro, cro2for, ran2cro, cro2ran, pas2cro, cro2pas, cro2urb, nfv2urb, \
        for2urb, ran2urb, pas2urb, for2pas, pas2for, nfv2pas, ran2pas, pas2nfv, for2ran, ran2for, for2nfv, nfv2for, mcgrath_frac_help = self.lucas_lut_input()
        for z in range(self.years):
            zz = self.years - z - 1
            print('year', zz)
            pft_help = self.lucas_lut_transrules(nfv2cro[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_shrubs, self.pfts_crops, 0, 0, self.nr_shrubs, self.nr_crops, 1, 1, self.pft_shrubs_default, shrubs_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(cro2nfv[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_crops, self.pfts_shrubs, self.pfts_grass, 0, self.nr_crops, self.nr_shrubs, self.nr_grass, 1, self.pft_crops_default, crops_backgr_help, self.backgrd, 2, False)
            pft_help = self.lucas_lut_transrules(for2cro[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_forest, self.pfts_crops, 0, 0, self.nr_forest, self.nr_crops, 1, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 1, self.mcgrath, mcgfrac=mcgrath_frac_help[zz, :, :, :])
            pft_help = self.lucas_lut_transrules(cro2for[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_crops, self.pfts_forest, self.pfts_shrubs, 0, self.nr_crops, self.nr_forest, self.nr_shrubs, 1, self.pft_crops_default, crops_backgr_help, False, 2, False)
            pft_help = self.lucas_lut_transrules(ran2cro[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_grass, self.pfts_crops, 0, 0, self.nr_grass, self.nr_crops, 1, 1, self.pft_grass_default, grass_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(cro2ran[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_crops, self.pfts_grass, 0, 0, self.nr_crops, self.nr_grass, 1, 1, self.pft_crops_default, crops_backgr_help, False, 1, False)
            pft_help = self.lucas_lut_transrules(pas2cro[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_grass, self.pfts_crops, 0, 0, self.nr_grass, self.nr_crops, 1, 1, self.pft_grass_default, grass_backgr_help, False, 1, False)
            pft_help = self.lucas_lut_transrules(cro2pas[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_crops, self.pfts_grass, 0, 0, self.nr_crops, self.nr_grass, 1, 1, self.pft_crops_default, crops_backgr_help, False, 1, False)
            pft_help = self.lucas_lut_transrules(cro2urb[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_crops, self.pfts_urban, 0, 0, self.nr_crops, self.nr_urban, 1, 1, self.pft_crops_default, crops_backgr_help, False, 1, False)
            pft_help = self.lucas_lut_transrules(nfv2urb[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_shrubs, self.pfts_urban, 0, 0, self.nr_shrubs, self.nr_urban, 1, 1, self.pft_shrubs_default, shrubs_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(for2urb[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_forest, self.pfts_urban, 0, 0, self.nr_forest, self.nr_urban, 1, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(ran2urb[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_shrubs_grass, self.pfts_urban, 0, 0, self.nr_shrubs+self.nr_grass, self.nr_urban, 1, 1, self.pft_shrubs_default, shrubs_grass_backgr, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(pas2urb[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_grass, self.pfts_urban, 0, 0, self.nr_grass, self.nr_urban, 1, 1, self.pft_grass_default, grass_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(for2pas[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_forest, self.pfts_grass, 0, 0, self.nr_forest, self.nr_grass, 1, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 1, self.mcgrath, mcgfrac=mcgrath_frac_help[zz, :, :, :])
            pft_help = self.lucas_lut_transrules(pas2for[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_grass, self.pfts_forest, 0, 0, self.nr_grass, self.nr_forest, 1, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(nfv2pas[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_shrubs, self.pfts_grass, 0, 0, self.nr_shrubs, self.nr_grass, 1, 1, self.pft_shrubs_default, shrubs_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(ran2pas[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_shrubs, self.pfts_grass, 0, 0, self.nr_shrubs, self.nr_grass, 1, 1, self.pft_shrubs_default, shrubs_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(pas2nfv[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_grass, self.pfts_shrubs, 0, 0, self.nr_grass, self.nr_shrubs, 1, 1, self.pft_grass_default, grass_backgr_help, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(for2ran[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_forest, self.pfts_shrubs, self.pfts_grass, 0, self.nr_forest, self.nr_shrubs, self.nr_grass, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 2, self.mcgrath, mcgfrac=mcgrath_frac_help[zz, :, :, :])
            pft_help = self.lucas_lut_transrules(ran2for[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_shrubs_grass, self.pfts_forest, 0, 0, self.nr_grass+self.nr_shrubs, self.nr_shrubs, 1, 1, self.pft_shrubs_default, shrubs_grass_backgr, self.backgrd, 1, False)
            pft_help = self.lucas_lut_transrules(for2nfv[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_forest, self.pfts_shrubs, self.pfts_grass, 0, self.nr_forest, self.nr_shrubs, self.nr_grass, 1, self.pft_forest_default, forest_backgr_help, self.backgrd, 2, self.mcgrath, mcgfrac=mcgrath_frac_help[zz, :, :, :])
            pft_help = self.lucas_lut_transrules(nfv2for[zz, :, :].data.T, rcm_lsm_help, pft_help, self.pfts_shrubs_grass, self.pfts_forest, 0, 0, self.nr_grass+self.nr_shrubs, self.nr_forest, 1, 1, self.pft_shrubs_default, shrubs_grass_backgr, self.backgrd, 1, False)
            self.pft_frac_ts[:, :, :, zz] = pft_help[:, :, :]
        if self.plot:
            self.plot_pft_frac_ts()
        print('LAND USE CHANGE FINISHED')
        # NORMALIZE TO GET A SUM OF 1 AND SET SEA POINTS TO MISSING VALUE
        self.recalc_pft_frac_ts(rcm_lsm)
        if self.irri:
            self.lucas_lut_irrigation(rcm_lsm)
        if self.mcgrath:
            self.lucas_lut_mcgrath(rcm_lsm)

    def lucas_lut_transrules(self, trans, rcm_lsm, inpfts, pfts_1, pfts_2, pfts_3, pfts_4,
                             nr_pfts_1, nr_pfts_2, nr_pfts_3, nr_pfts_4, defaultpft,
                             backgrdpfts, backgrd, rule, mcgrath, mcgfrac=np.zeros(3)):
        pft_1_sum = sum(inpfts[:, :, pfts_1[n]] for n in range(nr_pfts_1))
        pft_2_sum = sum(inpfts[:, :, pfts_2[n]] for n in range(nr_pfts_2))
        pft_3_sum = sum(inpfts[:, :, pfts_3[n]] for n in range(nr_pfts_3)) if rule >= 2 else 0
        pft_4_sum = sum(inpfts[:, :, pfts_4[n]] for n in range(nr_pfts_4)) if rule >= 3 else 0
        mcg_sum = sum(mcgfrac[n] for n in range(3)) if mcgrath else 0
        helper = np.zeros((self.xsize, self.ysize))
        helper_2 = np.zeros((self.xsize, self.ysize))
        helper_3 = np.zeros((self.xsize, self.ysize))
        inpfts_static_values = np.where(rcm_lsm > 0.0, False, inpfts)
        if rule == 1:
            selected_pft_1_sum = pft_1_sum[pft_1_sum > 0.0]
            selected_pft_2_sum = pft_2_sum[pft_2_sum > 0.0]
            trans = np.minimum(trans, trans - np.maximum(0.0, pft_1_sum + trans - 1), trans - np.maximum(0.0, trans - pft_2_sum))
            selected_trans = trans[trans > 0.0]
            if selected_trans.size > 0:
                if selected_pft_2_sum.size > 0:
                    filtered_pft_2_sum = np.where(pft_2_sum > 0.0, pft_2_sum, 1)
                    filtered_trans = np.where(trans > 0.0, trans, 1)
                    filtered_pft_2_sum = np.where(filtered_trans == 1, 1, filtered_pft_2_sum)
                    filtered_trans = np.where(filtered_pft_2_sum == 1, 1, filtered_trans)
                    for ipft in range(nr_pfts_2):
                        inpfts[:, :, pfts_2[ipft]] -= (inpfts[:, :, pfts_2[ipft]] / filtered_pft_2_sum * filtered_trans)
                        inpfts_pfts_2 = inpfts[:, :, pfts_2[ipft]]
                        selected_inpfts = inpfts_pfts_2[inpfts_pfts_2 < 0.0]
                        if selected_inpfts.size > 0:
                            filtered_inpfts = np.where(inpfts_pfts_2 < 0.0, inpfts_pfts_2, 0)
                            helper += filtered_inpfts
                            inpfts[:, :, pfts_2[ipft]] = np.where(inpfts_pfts_2 < 0.0, 0.0, inpfts_pfts_2)
                trans = np.maximum(0.0, trans + helper)
                if mcgrath and mcg_sum == 1:
                    filtered_trans = np.where(trans > 0, trans, 0)
                    for ipft in range(2, 5):
                        inpfts[:, :, pfts_1[ipft]] += (mcgfrac[:, :, ipft-2].T * filtered_trans)
                else:
                    if selected_pft_1_sum.size > 0:
                        filtered_pft_1_sum = np.where(pft_1_sum > 0.0, pft_1_sum, 1)
                        filtered_trans = np.where(trans > 0, trans, 1)
                        filtered_pft_1_sum = np.where(filtered_trans == 1, 1, filtered_pft_1_sum)
                        filtered_trans = np.where(filtered_pft_1_sum == 1, 1, filtered_trans)
                        for ipft in range(nr_pfts_1):
                            inpfts[:, :, pfts_1[ipft]] += (inpfts[:, :, pfts_1[ipft]] / filtered_pft_1_sum * filtered_trans)
                    else:
                        if backgrd:
                            filtered_trans = np.where(trans > 0, trans, 0)
                            for ipft in range(nr_pfts_1):
                                inpfts[:, :, pfts_1[ipft]] += (backgrdpfts[:, :, ipft] * filtered_trans)
                            
                        else:
                            inpfts[:, :, defaultpft] = np.where(trans > 0, trans, inpfts[:, :, defaultpft])
        elif rule == 2:
            selected_pft_1_sum = pft_1_sum[pft_1_sum > 0.0]
            selected_pft_2_sum = pft_2_sum[pft_2_sum > 0.0]
            selected_pft_3_sum = pft_3_sum[pft_3_sum > 0.0]
            trans = np.minimum(trans, trans - np.maximum(0.0, pft_1_sum + trans - 1), trans - np.maximum(0.0, trans - pft_2_sum + pft_3_sum))
            selected_trans = trans[trans > 0.0]
            
            if selected_trans.size > 0:
                if selected_pft_2_sum.size > 0:
                    filtered_pft_2_sum = np.where(pft_2_sum > 0.0, pft_2_sum, 1)
                    filtered_trans = np.where(trans > 0, trans, 1)
                    filtered_pft_2_sum = np.where(filtered_trans == 1, 1, filtered_pft_2_sum)
                    filtered_trans = np.where(filtered_pft_2_sum == 1, 1, filtered_trans)
                    for ipft in range(nr_pfts_2):
                        inpfts[:, :, pfts_2[ipft]] -= (inpfts[:, :, pfts_2[ipft]] / filtered_pft_2_sum * filtered_trans)
                        inpfts_pfts_2 = inpfts[:, :, pfts_2[ipft]]
                        selected_inpfts = inpfts_pfts_2[inpfts_pfts_2 < 0.0]
                        if selected_inpfts.size > 0:
                            filtered_inpfts = np.where(inpfts_pfts_2 < 0.0, inpfts_pfts_2, 0)
                            helper += filtered_inpfts
                            inpfts[:, :, pfts_2[ipft]] = np.where(inpfts_pfts_2 < 0.0, 0, inpfts_pfts_2)
                else:
                    helper = np.where(trans > 0, trans, helper)
                selected_helper = helper[helper > 0.0]
                if selected_helper.size > 0 and selected_pft_3_sum.size > 0:
                    filtered_pft_3_sum = np.where(pft_3_sum > 0.0, pft_3_sum, 1)
                    filtered_helper = np.where(helper > 0, helper, 1)
                    filtered_pft_3_sum = np.where(filtered_helper == 1, 1, filtered_pft_3_sum)
                    filtered_helper = np.where(filtered_pft_3_sum == 1, 1, filtered_helper)
                    for ipft in range(nr_pfts_3):
                        inpfts[:, :, pfts_3[ipft]] -= (inpfts[:, :, pfts_3[ipft]] / filtered_pft_3_sum * filtered_helper)
                        inpfts_pfts_3 = inpfts[:, :, pfts_3[ipft]]
                        selected_inpfts = inpfts_pfts_3[inpfts_pfts_3 < 0.0]
                        if selected_inpfts.size > 0:
                            filtered_inpfts = np.where(inpfts_pfts_3 < 0.0, inpfts_pfts_3, 0)
                            helper_2 += filtered_inpfts
                            inpfts[:, :, pfts_3[ipft]] = np.where(inpfts_pfts_3 < 0.0, 0, inpfts_pfts_3)
                            
                trans = np.maximum(0.0, trans + helper_2)
                filtered_trans = np.where(trans > 0, trans, 0)
                if mcgrath and mcg_sum == 1:
                    for ipft in range(2, 5):
                        inpfts[:, :, pfts_1[ipft]] += (mcgfrac[:, :, ipft - 2].T * filtered_trans)
                else:
                    if selected_pft_1_sum.size > 0:
                        filtered_pft_1_sum = np.where(pft_1_sum > 0.0, pft_1_sum, 1)
                        filtered_trans = np.where(trans > 0, trans, 1)
                        filtered_pft_1_sum = np.where(filtered_trans == 1, 1, filtered_pft_1_sum)
                        filtered_trans = np.where(filtered_pft_1_sum == 1, 1, filtered_trans)
                        for ipft in range(nr_pfts_1):
                            inpfts[:, :, pfts_1[ipft]] += (inpfts[:, :, pfts_1[ipft]] / filtered_pft_1_sum * filtered_trans)
                    else:
                        if backgrd:
                            filtered_trans = np.where(trans > 0, trans, 0)
                            for ipft in range(nr_pfts_1):
                                inpfts[:, :, pfts_1[ipft]] += (backgrdpfts[:, :, ipft] * filtered_trans)
                        else:
                            inpfts[:, :, defaultpft] = np.where(trans > 0, trans, inpfts[:, :, defaultpft])
        elif rule == 3:
            selected_pft_1_sum = pft_1_sum[pft_1_sum > 0.0]
            selected_pft_2_sum = pft_2_sum[pft_2_sum > 0.0]
            selected_pft_3_sum = pft_3_sum[pft_3_sum > 0.0]
            selected_pft_4_sum = pft_4_sum[pft_4_sum > 0.0]
            
            trans = np.minimum(trans, trans-np.maximum(0., pft_1_sum+trans-1), trans-np.maximum(0., trans-(pft_2_sum+pft_3_sum+pft_4_sum)))
            selected_trans = trans[trans > 0.0]
            if selected_trans.size > 0:
                if selected_pft_2_sum.size > 0:
                    filtered_pft_2_sum = np.where(pft_2_sum > 0.0, pft_2_sum, 1)
                    filtered_trans = np.where(trans > 0, trans, 1)
                    filtered_pft_2_sum = np.where(filtered_trans == 1, 1, filtered_pft_2_sum)
                    filtered_trans = np.where(filtered_pft_2_sum == 1, 1, filtered_trans)
                    for ipft in range(nr_pfts_2):
                        inpfts[:, :, pfts_2[ipft]] -= (inpfts[:, :, pfts_2[ipft]]/filtered_pft_2_sum*filtered_trans)
                        inpfts_pfts_2 = inpfts[:, :, pfts_2[ipft-1]]
                        selected_inpfts = inpfts_pfts_2[inpfts_pfts_2 < 0.0]
                        if selected_inpfts.size > 0:
                            filtered_inpfts = np.where(inpfts_pfts_2 < 0.0, inpfts_pfts_2, 0)
                            helper -= filtered_inpfts
                            inpfts[:, :, pfts_2[ipft]] = np.where(inpfts_pfts_2 < 0.0, 0.0, inpfts_pfts_2)
                else:
                    helper = np.where(trans > 0, trans, helper)
                selected_helper = helper[helper > 0.0]
                if selected_helper.size > 0:
                    if selected_pft_3_sum.size > 0:
                        for ipft in range(nr_pfts_3):
                            filtered_pft_3_sum = np.where(pft_3_sum > 0.0, pft_3_sum, 1)
                            filtered_helper = np.where(helper > 0, helper, 1)
                            filtered_pft_3_sum = np.where(filtered_helper == 1, 1, filtered_pft_3_sum)
                            filtered_helper = np.where(filtered_pft_3_sum == 1, 1, filtered_helper)
                            inpfts[:, :, pfts_3[ipft]] -= (inpfts[:, :, pfts_3[ipft]]/filtered_pft_3_sum*filtered_helper)
                            inpfts_pfts_3 = inpfts[:, :, pfts_3[ipft]]
                            selected_inpfts = inpfts_pfts_3[inpfts_pfts_3 < 0.0]
                            if selected_inpfts.size > 0:
                                filtered_inpfts = np.where(inpfts_pfts_3 < 0.0, inpfts_pfts_3, 0)
                                helper_2 -= filtered_inpfts
                                inpfts[:, :, pfts_3[ipft]] = np.where(inpfts_pfts_3 < 0.0, 0.0, inpfts_pfts_3)
                    else:
                        helper_2 = np.where(helper > 0, helper, helper_2)
                selected_helper_2 = helper_2[helper_2 > 0.0] 
                if selected_helper_2.size > 0:
                    if selected_pft_4_sum.size > 0:
                        filtered_pft_4_sum = np.where(pft_4_sum > 0.0, pft_4_sum, 1)
                        filtered_helper_2 = np.where(helper_2 > 0, helper_2, 1)
                        filtered_pft_4_sum = np.where(filtered_helper_2 == 1, 1, filtered_pft_4_sum)
                        filtered_helper_2 = np.where(filtered_pft_4_sum == 1, 1, filtered_helper_2)
                        for ipft in range(nr_pfts_4):
                            inpfts[:, :, pfts_4[ipft]] -= (inpfts[:, :, pfts_4[ipft]]/filtered_pft_4_sum*filtered_helper_2)
                            inpfts_pfts_4 = inpfts[:, :, pfts_4[ipft]]
                            selected_inpfts = inpfts_pfts_4[inpfts_pfts_4 < 0.0]
                            if selected_inpfts.size > 0:
                                filtered_inpfts = np.where(inpfts_pfts_4 < 0.0, inpfts_pfts_4, 0)
                                helper_3 += filtered_inpfts
                                inpfts[:, :, pfts_4[ipft]] = np.where(inpfts_pfts_4 < 0.0, 0.0, inpfts_pfts_4)
            # adjust transition if needed (not enough of pft group 2), helper is always <= 0
            trans = np.maximum(0., trans+helper_3)
            # add to pft group 2
            if mcgrath and mcg_sum == 1:  # if mcgrath forest data should be used (hard coded for
                for ipft in range(2, 5):
                    inpfts[:, :, pfts_1[ipft]] += (mcgfrac[:, :, ipft-2].T*trans)
            else:  # just use the relative fractions
                if selected_pft_1_sum.size > 0:
                    for ipft in range(nr_pfts_1):
                        filtered_pft_1_sum = np.where(pft_1_sum > 0.0, pft_1_sum, 1)
                        filtered_trans = np.where(trans > 0, trans, 1)
                        filtered_pft_1_sum = np.where(filtered_trans == 1, 1, filtered_pft_1_sum)
                        filtered_trans = np.where(filtered_pft_1_sum == 1, 1, filtered_trans)
                        inpfts[:, :, pfts_1[ipft]] += (inpfts[:, :, pfts_1[ipft]]/filtered_pft_1_sum*filtered_trans)
                if backgrd:
                    filtered_trans = np.where(trans > 0, trans, 0)
                    for ipft in range(nr_pfts_1):
                        inpfts[:, :, pfts_1[ipft]] += (backgrdpfts[:, :, ipft] * filtered_trans)
                else:
                    inpfts[:, :, defaultpft] = np.where(trans > 0, trans, inpfts[:, :, defaultpft])
        inpfts = np.where(inpfts_static_values == False, inpfts, inpfts_static_values)
        return inpfts

    def lucas_lut_input(self):
        pft_help = np.zeros((self.xsize, self.ysize, self.npfts))
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
            urban_backgr_help = np.zeros((self.xsize, self.ysize, self.nr_urban))
            for i in range(self.nr_urban):
                num = i + 15
                urban_backgr_help[:, :, i] = urban_backgr[f"var8{num}"].data[0, :, :].T
        grass_backgr = xr.open_dataset(self.namelist["F_BACKGRA"])
        grass_backgr_help = np.zeros((self.xsize, self.ysize, self.nr_grass))
        for i in range(self.nr_grass):
            num = i + 9
            if len(str(num)) > 1:
                grass_backgr_help[:, :, i] = grass_backgr[f"var8{num}"].data[0, :, :].T
            else:
                grass_backgr_help[:, :, i] = grass_backgr[f"var80{num}"].data[0, :, :].T

        crops_backgr = xr.open_dataset(self.namelist["F_BACKCRO"])
        crops_backgr_help = np.zeros((self.xsize, self.ysize, self.nr_crops))
        for i in range(self.nr_crops):
            num = i + 13
            if len(str(num)) > 1:
                crops_backgr_help[:, :, i] = crops_backgr[f"var8{num}"].data[0, :, :].T
            else:
                crops_backgr_help[:, :, i] = crops_backgr[f"var80{num}"].data[0, :, :].T

        forest_backgr = xr.open_dataset(self.namelist["F_BACKFOR"])
        forest_backgr_help = np.zeros((self.xsize, self.ysize, self.nr_forest))
        for i in range(self.nr_forest):
            num = i + 1
            if len(str(num)) > 1:
                forest_backgr_help[:, :, i] = forest_backgr[f"var8{num}"].data[0, :, :].T
            else:
                forest_backgr_help[:, :, i] = forest_backgr[f"var80{num}"].data[0, :, :].T

        shrubs_backgr = xr.open_dataset(self.namelist["F_BACKSHR"])
        shrubs_backgr_help = np.zeros((self.xsize, self.ysize, self.nr_shrubs))
        for i in range(self.nr_shrubs):
            num = i + 7
            if len(str(num)) > 1:
                shrubs_backgr_help[:, :, i] = shrubs_backgr[f"var8{num}"].data[0, :, :].T
            else:
                shrubs_backgr_help[:, :, i] = shrubs_backgr[f"var80{num}"].data[0, :, :].T
    
        shrubs_grass_backgr = np.zeros((self.xsize, self.ysize, self.nr_shrubs+self.nr_grass))
        shrubs_grass_backgr[:, :, :self.nr_shrubs] = shrubs_backgr_help / 2.
        shrubs_grass_backgr[:, :, self.nr_shrubs:(self.nr_grass+self.nr_shrubs)] = grass_backgr_help / 2.

        # Mcfrac dataset
        mcgrath_frac = netCDF4.Dataset(self.namelist["F_MCGRATH"], decode_times=False)
        mcgrath_frac_help = np.zeros((self.years+1, self.ysize, self.xsize, 3))
        for i in range(1, 4):
            mcgrath_frac_help[:, :, :, i-1] = mcgrath_frac[f"var80{i+2}"][:, :, :].data

        # RCM LSM
        rcm_lsm = xr.open_dataset(self.namelist["F_RCM_LSM_IN"]).var210.values.T
        rcm_lsm_help = np.zeros((self.xsize, self.ysize, self.npfts))

        # Transformation datasets
        nfv2cro = netCDF4.Dataset(self.namelist["F_NFV2CRO"], decode_times=False).variables["nfv2cro"]
        cro2nfv = netCDF4.Dataset(self.namelist["F_CRO2NFV"], decode_times=False).variables["cro2nfv"]
        for2cro = netCDF4.Dataset(self.namelist["F_FOR2CRO"], decode_times=False).variables["for2cro"]
        cro2for = netCDF4.Dataset(self.namelist["F_CRO2FOR"], decode_times=False).variables["cro2for"]
        ran2cro = netCDF4.Dataset(self.namelist["F_RAN2CRO"], decode_times=False).variables["ran2cro"]
        cro2ran = netCDF4.Dataset(self.namelist["F_CRO2RAN"], decode_times=False).variables["cro2ran"]
        pas2cro = netCDF4.Dataset(self.namelist["F_PAS2CRO"], decode_times=False).variables["pas2cro"]
        cro2pas = netCDF4.Dataset(self.namelist["F_CRO2PAS"], decode_times=False).variables["cro2pas"]
        cro2urb = netCDF4.Dataset(self.namelist["F_CRO2URB"], decode_times=False).variables["cro2urb"]
        nfv2urb = netCDF4.Dataset(self.namelist["F_NFV2URB"], decode_times=False).variables["nfv2urb"]
        for2urb = netCDF4.Dataset(self.namelist["F_FOR2URB"], decode_times=False).variables["for2urb"]
        ran2urb = netCDF4.Dataset(self.namelist["F_RAN2URB"], decode_times=False).variables["ran2urb"]
        pas2urb = netCDF4.Dataset(self.namelist["F_PAS2URB"], decode_times=False).variables["pas2urb"]
        for2pas = netCDF4.Dataset(self.namelist["F_FOR2PAS"], decode_times=False).variables["for2pas"]
        pas2for = netCDF4.Dataset(self.namelist["F_PAS2FOR"], decode_times=False).variables["pas2for"]
        nfv2pas = netCDF4.Dataset(self.namelist["F_NFV2PAS"], decode_times=False).variables["nfv2pas"]
        ran2pas = netCDF4.Dataset(self.namelist["F_RAN2PAS"], decode_times=False).variables["ran2pas"]
        pas2nfv = netCDF4.Dataset(self.namelist["F_PAS2NFV"], decode_times=False).variables["pas2nfv"]
        for2ran = netCDF4.Dataset(self.namelist["F_FOR2RAN"], decode_times=False).variables["for2ran"]
        ran2for = netCDF4.Dataset(self.namelist["F_RAN2FOR"], decode_times=False).variables["ran2for"]
        for2nfv = netCDF4.Dataset(self.namelist["F_FOR2NFV"], decode_times=False).variables["for2nfv"]
        nfv2for = netCDF4.Dataset(self.namelist["F_NFV2FOR"], decode_times=False).variables["nfv2for"]

        if self.forward:
            urb2cro = netCDF4.Dataset(self.namelist["F_URB2CRO"], decode_times=False).variables["urb2cro"]
            urb2nfv = netCDF4.Dataset(self.namelist["F_URB2NFV"], decode_times=False).variables["urb2nfv"]
            urb2for = netCDF4.Dataset(self.namelist["F_URB2FOR"], decode_times=False).variables["urb2for"]
            urb2ran = netCDF4.Dataset(self.namelist["F_URB2RAN"], decode_times=False).variables["urb2ran"]
            urb2pas = netCDF4.Dataset(self.namelist["F_URB2PAS"], decode_times=False).variables["urb2pas"]
            nfv2ran = netCDF4.Dataset(self.namelist["F_NFV2RAN"], decode_times=False).variables["nfv2ran"]
            if self.addtree:
                nat2for = netCDF4.Dataset(self.namelist["F_NAT2FOR"], decode_times=False).variables["nat2for"]
                return rcm_lsm, rcm_lsm_help, pft_help, grass_backgr_help, crops_backgr_help, forest_backgr_help, shrubs_backgr_help, \
                       shrubs_grass_backgr, nfv2cro, cro2nfv, for2cro, cro2for, ran2cro, cro2ran, pas2cro, cro2pas, cro2urb, nfv2urb, \
                       for2urb, ran2urb, pas2urb, for2pas, pas2for, nfv2pas, ran2pas, pas2nfv, for2ran, ran2for, for2nfv, nfv2for, mcgrath_frac_help, \
                       urb2cro, urb2nfv, urb2for, urb2ran, urb2pas, nfv2ran, nat2for, urban_backgr_help
            
            return rcm_lsm, rcm_lsm_help, pft_help, grass_backgr_help, crops_backgr_help, forest_backgr_help, shrubs_backgr_help, \
               shrubs_grass_backgr, nfv2cro, cro2nfv, for2cro, cro2for, ran2cro, cro2ran, pas2cro, cro2pas, cro2urb, nfv2urb, \
               for2urb, ran2urb, pas2urb, for2pas, pas2for, nfv2pas, ran2pas, pas2nfv, for2ran, ran2for, for2nfv, nfv2for, mcgrath_frac_help, \
                urb2cro, urb2nfv, urb2for, urb2ran, urb2pas, nfv2ran, urban_backgr_help

        return rcm_lsm, rcm_lsm_help, pft_help, grass_backgr_help, crops_backgr_help, forest_backgr_help, shrubs_backgr_help, \
               shrubs_grass_backgr, nfv2cro, cro2nfv, for2cro, cro2for, ran2cro, cro2ran, pas2cro, cro2pas, cro2urb, nfv2urb, \
               for2urb, ran2urb, pas2urb, for2pas, pas2for, nfv2pas, ran2pas, pas2nfv, for2ran, ran2for, for2nfv, nfv2for, mcgrath_frac_help

    def lucas_lut_irrigation(self, rcm_lsm):
        irri_frac = xr.open_dataset(self.namelist["F_IRRI_IN"])
        mask = (rcm_lsm > 0.0) & ((self.pft_frac_ts[:, :, 12, :] + self.pft_frac_ts[:, :, 13, :]) > 0.0)
        sum_crops = self.pft_frac_ts[:, :, 12, :] + self.pft_frac_ts[:, :, 13, :]
        irri_mask = irri_frac > 0.0
        self.pft_frac_ts[mask & irri_mask, 12, :] = (1.0 - irri_frac[mask & irri_mask]) * sum_crops[mask & irri_mask]
        self.pft_frac_ts[mask & irri_mask, 13, :] = irri_frac[mask & irri_mask] * sum_crops[mask & irri_mask]
        self.pft_frac_ts[mask & ~irri_mask, 12, :] = sum_crops[mask & ~irri_mask]
        mask = rcm_lsm > 0.0
        self.pft_frac_ts = np.where(self.pft_frac_ts < 0.0, 0.0, self.pft_frac_ts)
        pft_sum = np.sum(self.pft_frac_ts, axis=2)
        pft_sum_mask = pft_sum > 0.0
        self.pft_frac_ts[mask & pft_sum_mask] /= pft_sum[mask & pft_sum_mask, np.newaxis]

    def recalc_pft_frac_ts(self, rcm_lsm):
        print_section_heading('NORMALIZE TO GET A SUM OF 1 AND SET SEA POINTS TO MISSING VALUE')
        for z in range(self.years+1):
            mask = rcm_lsm > 0.0
            pft_sum = np.sum(self.pft_frac_ts[:, :, :, z], axis=2)
            pft_sum_mask = pft_sum > 0.0
            self.pft_frac_ts[mask & pft_sum_mask, :, z] /= pft_sum[mask & pft_sum_mask, np.newaxis]
            self.pft_frac_ts[~mask, :, z] = -999.

    def lucas_lut_mcgrath(self, rcm_lsm):
        mcgrath_frac = xr.open_dataset(self.namelist["F_MCGRATH"])
        for z in range(self.years+1):
            zz = years + 1 - z
            d_mcg_frac = [0.0, 0.0, 0.0]
            d_frac = [0.0, 0.0, 0.0]
            d_forest = [0.0, 0.0, 0.0]
            sum_forest = 0.0
            sum_forest_p1 = 0.0
            sum_helper = 0.0
            sum_mcg = 0.0
            for k in range(2, 5):
                sum_forest += pft_frac_ts[:, :, k, zz]
                sum_forest_p1 += pft_frac_ts[:, :, k, zz]
                sum_mcg += mcgrath_frac[:, :, k-2, zz]
            for k in range(2, 5):
                if mcgrath_frac[:, :, k-2, :] < -998.0 or mcgrath_frac[i-1][j-1][k-2][zz] < -998.0:
                    sum_mcg = 0.0
            if sum_forest_p1 < 0.0 or sum_forest < 0.0:
                pass  # Handle problem condition
            if sum_forest_p1 > 0.0 and sum_mcg > 0.0 and sum_forest > 0.0:
                abssum = 0.0

                for k in range(3):
                    d_mcg_frac[k] = mcgrath_frac[i-1][j-1][k][zz] - mcgrath_frac[i-1][j-1][k][zz-1]
                    abssum += abs(d_mcg_frac[k])

                for k in range(3):
                    d_forest[k] = pft_frac_ts[i-1][j-1][k+2][zz] / sum_forest_p1 * (sum_forest_p1 - sum_forest)

                if abssum > 0.0:
                    for k in range(3):
                        d_frac[k] = sum_forest_p1 * d_mcg_frac[k]

                    for k in range(3):
                        pft_frac_ts[i-1][j-1][k+2][zz-1] = pft_frac_ts[i-1][j-1][k+2][zz] - d_frac[k] - d_forest[k]
                    lhelp = [1, 1, 1]
                    helper = 0.0

                    for k in range(3):
                        if pft_frac_ts[i-1][j-1][k+2][zz-1] < 0.0:
                            helper += pft_frac_ts[i-1][j-1][k+2][zz-1]
                            pft_frac_ts[i-1][j-1][k+2][zz-1] = 0.0
                            lhelp[k] = 0

                    sum_helper = 0.0
                    sum_lhelp = 0.0
                    for k in range(3):
                        sum_helper += lhelp[k]

                    if sum_helper < 3:

                        for k in range(3):
                            if lhelp[k] == 1:
                                sum_helper += pft_frac_ts[i-1][j-1][k+2][zz-1]

                        for k in range(3):
                            if sum_helper > 0.0 and lhelp[k] == 1:
                                pft_frac_ts[i-1][j-1][k+2][zz-1] += pft_frac_ts[i-1][j-1][k+2][zz-1] / sum_helper * helper
                            else:
                                pft_frac_ts[i-1][j-1][k+2][zz-1] -= (1 / 3) * helper
                else:
                    for k in range(3):
                        pft_frac_ts[i-1][j-1][k+2][zz-1] = pft_frac_ts[i-1][j-1][k+2][zz] - d_forest[k]

                    lhelp = [1, 1, 1]
                    helper = 0.0

                    for k in range(3):
                        if pft_frac_ts[i-1][j-1][k+2][zz-1] < 0.0:
                            helper += pft_frac_ts[i-1][j-1][k+2][zz-1]
                            pft_frac_ts[i-1][j-1][k+2][zz-1] = 0.0
                            lhelp[k] = 0

                    sum_helper = 0.0
                    sum_lhelp = 0.0
                    for k in range(3):
                        sum_helper += lhelp[k]

                    if sum_helper < 3:
                        for k in range(3):
                            if lhelp[k] == 1:
                                sum_helper += pft_frac_ts[i-1][j-1][k+2][zz-1]

                        for k in range(3):
                            if sum_helper > 0.0 and lhelp[k] == 1:
                                pft_frac_ts[i-1][j-1][k+2][zz-1] += pft_frac_ts[i-1][j-1][k+2][zz-1] / sum_helper * helper
                            else:
                                pft_frac_ts[i-1][j-1][k+2][zz-1] -= (1 / 3) * helper

        # remove zeros (they still occur due to rounding issues)
        

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

        lmcg = False
        if self.scenario == "historical":
            sdir = f"{luhdir}/historic/{self.region}/{self.grid}"
            if self.mcgback:
                lmcg = True
        elif self.scenario == "historical_high":
            sdir = f"{luhdir}/historic_high/{self.region}/{self.grid}"
            if self.mcgback:
                lmcg = True
        elif self.scenario == "historical_low":
            sdir = f"{luhdir}/historic_low/{self.region}/{self.grid}"
            if self.mcgback:
                lmcg = True
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
            "F_MCGRATH": f"{mcgdir}/{mcg}_{self.syear}_{self.eyear}_ForestBckgrdMcGrath_{self.grid}.nc", # mcgfile
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

            # CONTROL
            "XSIZE": xsize,
            "YSIZE": ysize,
            "MCGRATH": lmcg,
            "FORWARD": self.forward,
            "BACKGRD": self.backgrd,
            "ADDTREE": self.addtree,
            "IRRI": self.irri,
            "NPFTS" : 16,
            "GRADEF" : 9,
            "CRODEF": 13,
            "SHRDEF": 8,
            "FORDEF": 4,
            "URBDEF": 15,
            "FORPFTS": "1, 2, 3, 4, 5, 6, 0, 0, 0, 0",
            "SHRPFTS": "7, 8, 0, 0, 0, 0, 0, 0, 0, 0",
            "GRAPFTS": "9, 10, 11, 0, 0, 0, 0, 0, 0, 0",
            "CROPFTS": "13, 14, 0, 0, 0, 0, 0, 0, 0, 0",
            "URBPFTS": "15, 0, 0, 0, 0, 0, 0, 0, 0, 0",
            "CONPFTS": "12, 16, 0, 0, 0, 0, 0, 0, 0, 0"
        }
        cdo.setgrid(os.path.join(scriptsdir, f"grid_{self.grid}"), input=f"{odir}/{ofile}.srv", output=f"{odir}/{ofile}.nc", options=f"-f nc -settime,12:00:00 -setctomiss,-999") 
        return namelist_dict

    def prepare_luh2_data(self):
        """
        V2.0
        """
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
        print(f"Selecting time period {self.syear}/{self.eyear}")
        cdo.selyear(f"{self.syear}/{self.eyear}", input=f"{datadir}/{tfile}.nc", output=f"{datadir}/{tfile}_2.nc")
        print(f"Done")
        # still have to fix this part to make it lighter. 
        if self.trans:
            print(f"Selecting variables for transitions")
            ofile=f"{tfile}_{self.syear}_{self.eyear}_{self.region}.nc"
            cdo.sellonlatbox(self.reg, input=f"-selvar,{vars_trans} {datadir}/{tfile}_2.nc", output=f"{path_sdir}/{self.region}/{ofile}")
            print("Done")
        if self.state:
            print(f"Selecting variables for states")
            cdo.sellonlatbox(self.reg, input=f"-selyear,{self.syear}/{self.eyear} -selvar,{vars_state} {datadir}/{sfile}.nc", output=f"{path_region}/states_{self.syear}_{self.eyear}_{self.region}.nc")
            if remap_com == "invertlat":
                cdo.invertlat(input=f"{path_region}/states_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/states_{self.syear}_{self.eyear}_{self.grid}.nc")
            elif remap_com == "remapbil":
                cdo.remapbil(f"{scriptsdir}/grid_{self.grid}", input=f"{path_region}/states_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/states_{self.syear}_{self.eyear}_{self.grid}.nc")
            elif remap_com == "remapcon2":
                cdo.remapcon2(f"{scriptsdir}/grid_{self.grid}", input=f"{path_region}/states_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/states_{self.syear}_{self.eyear}_{self.grid}.nc")
            print("Done")
        if self.addtree:
            print(f"Selecting variables for added tree cover")
            cdo.sellonlatbox(self.reg, input=f"-selyear,{self.syear}/{self.eyear} -selvar,added_tree_cover {path_sdir}/{afile}.nc", output=f"{path_region}/addtree_{self.syear}_{self.eyear}_{self.region}.nc")
            if remap_com == "invertlat":
                cdo.invertlat(input=f"{path_region}/addtree_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/{self.grid}/addtree_{self.syear}_{self.eyear}_{self.grid}.nc")
            elif remap_com == "remapbil":
                cdo.remapbil(f"{scriptsdir}/grid_{self.grid}", input=f"{path_region}/addtree_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/{self.grid}/addtree_{self.syear}_{self.eyear}_{self.grid}.nc")
            elif remap_com == "remapcon2":
                cdo.remapcon2(f"{scriptsdir}/grid_{self.grid}", input=f"{path_region}/addtree_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/{self.grid}/addtree_{self.syear}_{self.eyear}_{self.grid}.nc")
            cdo.copy(input=f'{path_region}/{self.grid}/addtree_{self.syear}_{self.eyear}_{self.grid}.nc', output=f"{path_region}/{self.grid}/addtree_{self.syear}_{self.eyear}_{self.grid}.srv")
            print("Done")
        # compute irragtion fraction 
        if self.irri:
            print(f"Selecting variables for irrigation")
            cdo.sellonlatbox(self.reg, input=f"-selyear,{self.syear}/{self.eyear} -selvar,{vars_irrig} {path_sdir}/{mfile}.nc", output=f"{path_region}/irri_vars_{self.syear}_{self.eyear}_{self.region}.nc")

            if self.scenario in ["historical", "historical_low", "historical_high"]:
                if remap_com == "invertlat":
                    cdo.invertlat(input=f"-chname,{IRR[0]},irrig_frac -selvar,{IRR[0]} {path_region}/irri_vars_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/{grid}/irrigation_{syear}_{eyear}_{grid}.nc")
                elif remap_com == "remapbil":
                    cdo.remapbil(f"{scriptsdir}/grid_{grid}", input=f"-chname,{IRR[0]},irrig_frac -selvar,{IRR[0]} {path_region}/irri_vars_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/{self.grid}/irrigation_{syear}_{eyear}_{grid}.nc")
                elif remap_com == "remapcon2":
                    cdo.remapcon2(f"{scriptsdir}/grid_{grid}", input=f"-chname,{IRR[0]},irrig_frac -selvar,{IRR[0]} {path_region}/irri_vars_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/{self.grid}/irrigation_{syear}_{eyear}_{grid}.nc")
                cdo.copy(options="-setmisstoc,-999", input=f'{path_region}/{self.grid}/irrigation_{self.syear}_{self.eyear}_{self.grid}.nc', output=f"{path_region}/{self.grid}/irrigation_{self.syear}_{self.eyear}_$grid.srv")
            else:
                if remap_com in ["remapbil", "remapcon2"]:
                    cdo.setmisstoc(input=f"-999 -{remap_com},{scriptsdir}/grid_{self.grid} -varssum -selvar,{vars_crops} {path_region}/states_{self.syear}_{self.eyear}_{self.grid}.nc", output=f"{path_region}/sum_crop_frac.nc")
                else:
                    cdo.setmisstoc(input=f"-999 -{remap_com} -varssum -selvar,{vars_crops} {path_region}/states_{self.syear}_{self.eyear}_{self.grid}.nc", output=f"{path_region}/sum_crop_frac.nc")

                # Change variables names

                for n in range(5):
                    if remap_com in ["remapbil", "remapcon2"]:
                        cdo.mul(input=f"-{remap_com}, grid_{grid} -selvar,{IRR[n]} {path_region}/irri_vars_{self.syear}_{self.eyear}_{self.region}.nc -{remap_com}, grid_{self.grid} -selvar,{ICR[n]} {path_region}/states_{self.syear}_{self.eyear}_{self.grid}.nc", output=f"{sdir}/{self.region}/dummy.nc")
                    else:
                        cdo.mul(input=f"-{remap_com} -selvar,{IRR[n]} {path_region}/irri_vars_{self.syear}_{self.eyear}_{self.region}.nc -{remap_com} -selvar,{ICR[n]} {path_region}/states_{syear}_{eyear}_{grid}.nc", output=f"{sdir}/{self.region}/dummy.nc")
                    if n == 0:
                        cdo.chname(input=f"{IRR[0]},irrig_frac -selvar,{IRR[0]} {path_region}/dummy.nc", output=f"{path_region}/sum_irri_frac_{self.syear}_{self.eyear}_{self.grid}.nc")
                    else:
                        cdo.add(input=f"-chname,{IRR[n]},irrig_frac -selvar,{IRR[n]} {path_region}/dummy.nc {path_region}/sum_irri_frac_{self.syear}_{self.eyear}_{self.grid}.nc")   
                        shutil.move(f"{path_region}/dummy2.nc", f"{path_region}/sum_irri_frac_{self.syear}_{self.eyear}_{self.grid}.nc")
                    os.remove(f"{path_region}/dummy.nc")
                cdo.div(input=f"{path_region}/sum_irri_frac_{self.syear}_{self.eyear}_{self.grid}.nc {path_region}/sum_crop_frac.nc", output=f"{path_region}/{self.grid}/irrigation_{self.syear}_{self.eyear}_{self.grid}.nc")
                cdo.copy(input=f"-setmisstoc,-999 {path_region}/{self.grid}/irrigation_{self.syear}_{self.eyear}_{self.grid}.nc", output=f"{path_region}/{self.grid}/irrigation_{self.syear}_{self.eyear}_{self.grid}.srv")
            print("Done")
        if self.trans:
            # New classification
            print(f"Selecting variables for transitions")
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
            print("Done")

    def prepare_mcgrath(self):
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
        #
        # compute background for LUT classes using zonal mean
        #
        cdo.chname(f"maxvegetfrac,{PFT_TeBrEv}", input=f"-vertsum -sellevel,{TeBrEv} {ifile}", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_TeBrEv.nc")
        cdo.chname(f"maxvegetfrac,{PFT_TeBrDec}", input=f"-vertsum -sellevel,{TeBrDec} {ifile}", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_TeBrDec.nc")
        cdo.chname(f"maxvegetfrac,{PFT_ConEv}", input=f"-vertsum -sellevel,{EvCon} {ifile}", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_ConEv.nc")
        cdo.chname(f"maxvegetfrac,forest", input=f"-vertsum -sellevel,{TeBrEv},{TeBrDec},{EvCon} {ifile}", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_FOR.nc")
        cdo.merge(input=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_TeBrEv.nc {glcdir}/{self.lcd}_{self.syear}_{self.eyear}_TeBrDec.nc {glcdir}/{self.lcd}_{self.syear}_{self.eyear}_ConEv.nc", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_dummy.nc")
        if self.reg:
            cdo.setmisstoc(-999, input=f"-remapbil,{scriptsdir}/grid_{self.grid} -sellonlatbox,{self.reg} -div {glcdir}/{self.lcd}_{self.syear}_{self.eyear}_dummy.nc -varssum {glcdir}/{self.lcd}_{self.syear}_{self.eyear}_dummy.nc", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_ForestBckgrdMcGrath_{self.grid}.nc")
        else:
            cdo.setmisstoc(-999, input=f"-remapbil,{scriptsdir}/grid_{self.grid} -div {glcdir}/{self.lcd}_{self.syear}_{self.eyear}_dummy.nc -varssum {glcdir}/{self.lcd}_{self.syear}_{self.eyear}_dummy.nc", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_ForestBckgrdMcGrath_{self.grid}.nc")
        if self.mcgrath_eyear < self.eyear:
            for year in range(self.mcgrath_eyear, self.eyear+1):
                cdo.setdate(f"{year}-06-15", input=f"-selyear,2010 {glcdir}/{self.lcd}_{self.syear}_{self.eyear}_ForestBckgrdMcGrath_{self.grid}.nc", output=f"{glcdir}/dummy_{year}.nc")
        cdo.mergetime(input=f"{glcdir}/dummy_????.nc", output=f"{glcdir}/{self.lcd}_2011_2015_ForestBckgrdMcGrath_{self.grid}.nc")
        cdo.mergetime(input=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_ForestBckgrdMcGrath_{self.grid}.nc {glcdir}/{self.lcd}_2011_2015_ForestBckgrdMcGrath_{self.grid}.nc", output=f"{glcdir}/{self.lcd}_{self.syear}_2015_ForestBckgrdMcGrath_{self.grid}.nc")
        cdo.copy(input=f"{glcdir}/{self.lcd}_{self.syear}_2015_ForestBckgrdMcGrath_{self.grid}.nc", output=f"{glcdir}/{self.lcd}_{self.syear}_2015_ForestBckgrdMcGrath_{self.grid}.srv")

    def fromto(self, varn, for_1, for_2, tfile, ext, cutting, path, remap_com, outvar_condition=None):
        odir = self.grid
        # combine land-use changes using reclassifcation
        ifile=f"{tfile}_{self.syear}_{self.eyear}_{self.region}.nc"
        logfile=f"{tfile}_{self.syear}_{self.eyear}_{self.region}.log"
        ofile=f"transitions_{self.syear}_{self.eyear}_{self.region}_{varn}"
        cdo.mulc(0, input=f"-selvar,primf_to_urban {path}/{ifile}", output=f"{path}/dummy.nc")
        for inivar in for_1:
            for outvar in for_2:
                if not outvar_condition:
                    cdo.add(input=f"-selvar,{inivar}_to_{outvar} {path}/{ifile} {path}/dummy.nc", output=f"{path}/{ofile}.nc")
                    cdo.chname(f"{inivar}_to_{outvar},{varn}", input=f"{path}/{ofile}.nc", output=f"{path}/dummy.nc")
                if outvar != outvar_condition:
                    cdo.add(input=f"-selvar,{inivar}_to_{outvar} {path}/{ifile} {path}/dummy.nc", output=f"{path}/{ofile}.nc")
                    cdo.chname(f"{inivar}_to_{outvar},{varn}", input=f"{path}/{ofile}.nc", output=f"{path}/dummy.nc")
        shutil.move(f"{path}/dummy.nc", f"{path}/{ofile}.nc")
        if remap_com == "invertlat":
            cdo.invertlat(input=f"{path}/{ofile}.nc", output=f"{path}/{self.grid}/{ofile}_{ext}_{self.grid}.nc")
        elif remap_com == "remapbil":
            cdo.remapbil(f"{scriptsdir}/grid_{self.grid}", input=f"{path}/{ofile}.nc", output=f"{path}/{self.grid}/{ofile}_{ext}_{self.grid}.nc")
        elif remap_com == "remapcon2":
            cdo.remapcon2(f"{scriptsdir}/grid_{self.grid}", input=f"{path}/{ofile}.nc", output=f"{path}/{self.grid}/{ofile}_{ext}_{self.grid}.nc")
        cdo.copy(options=f"{cutting} -setmisstoc,-999.", input=f'{path}/{self.grid}/{ofile}_{ext}_{self.grid}.nc', output=f"{path}/{self.grid}/{ofile}_{ext}_{self.grid}.srv")

    def plot_pft_frac_ts(self):
        year = self.plot_year if self.plot_year else 1
        Path(os.path.join(plotdir, self.region, str(self.syear+year))).mkdir(parents=True, exist_ok=True)
        if self.plot_npft:
            self.main_plot(self.plot_npft, year)
        else:
            for inpft in range(self.npfts):
                self.main_plot(inpft+1, year)

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
        ax.add_feature(cfeature.LAND, facecolor="white")
        ax.add_feature(cfeature.LAKES, facecolor='blue')
        ax.add_feature(cfeature.RIVERS, edgecolor='blue')
        # Plot the data
        min_value = np.min(self.pft_frac_ts[:, :, npft-1, year-1])
        max_value = np.max(self.pft_frac_ts[:, :, npft-1, year-1])
        #levels = np.linspace(min_value, max_value, 11)
        levels = np.linspace(0, 1, 11)  # 10 intervals
        colors = ['#ffffff', '#f7fcf5', '#e5f5e0', '#c7e9c0', '#a1d99b', '#74c476', '#41ab5d', '#238b45', '#006d2c', '#00441b']
        cmap = ListedColormap(colors)
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        contour = ax.contourf(lon_grid, lat_grid, self.pft_frac_ts[:, :, npft-1, year-1].T, levels=levels, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm, alpha=1)
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
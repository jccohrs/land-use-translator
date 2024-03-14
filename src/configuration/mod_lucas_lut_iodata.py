# Arrays for input and output data
import numpy as np
from conf.conf import NPFTS, NR_GRASS, NR_CROPS, NR_SHRUBS, NR_FOREST, NR_URBAN

def return_iodata_arrays(XSIZE, YSIZE, YEARS):
    rcm_lsm = np.zeros((XSIZE, YSIZE)) # input: labd-sea mask
    pft_frac = np.zeros((XSIZE, YSIZE, NPFTS)) # input: land cover fractions
    pft_frac_ts = np.zeros((XSIZE, YSIZE, NPFTS, YEARS+1)) # output: time series of land cover fractions
    grass_backgr = np.zeros((XSIZE, YSIZE, NR_GRASS)) # background map for grass 
    crops_backgr = np.zeros((XSIZE, YSIZE, NR_CROPS)) # background map for crops
    shrubs_backgr = np.zeros((XSIZE, YSIZE, NR_SHRUBS)) # background map for shrubs
    forest_backgr = np.zeros((XSIZE, YSIZE, NR_FOREST)) # background map for forest
    urban_backgr = np.zeros((XSIZE, YSIZE, NR_URBAN)) # background map for forest
    irri_frac = np.zeros((XSIZE, YSIZE, YEARS+1)) # background map for forest
    mcgrath_frac = np.zeros((XSIZE, YSIZE, 3, YEARS+1)) # forest data from McGrath

    return rcm_lsm, pft_frac, pft_frac_ts, grass_backgr, crops_backgr, shrubs_backgr, forest_backgr, urban_backgr, irri_frac, mcgrath_frac

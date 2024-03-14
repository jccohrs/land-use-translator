# Arrays of land cover transition time series
import numpy as np
from conf.conf import FORPFTS, SHRPFTS, GRAPFTS, CROPFTS, URBPFTS, CONPFTS

def return_trans_arrays(XSIZE, YSIZE, YEARS):
    NR_GRASS = 0
    NR_CROPS = 0
    NR_SHRUBS = 0
    NR_FOREST = 0
    NR_URBAN = 1
    PFTS_GRASS = list(filter(0, GRAPFTS))
    PFTS_CROPS = list(filter(0, CROPFTS))
    PFTS_SHRUBS = list(filter(0, SHRPFTS))
    PFTS_FOREST = list(filter(0, FORPFTS))
    PFTS_URBAN = list(filter(0, URBPFTS))
    return NR_GRASS, NR_CROPS, NR_SHRUBS, NR_FOREST, NR_URBAN, PFTS_GRASS, PFTS_CROPS, PFTS_SHRUBS, PFTS_FOREST, PFTS_URBAN

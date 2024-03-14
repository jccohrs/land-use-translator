# Arrays of land cover transition time series
import numpy as np

def return_trans_arrays(XSIZE, YSIZE, YEARS):
    for2cro = np.zeros((XSIZE, YSIZE, YEARS)) # forest to crops
    cro2for = np.zeros((XSIZE, YSIZE, YEARS)) # crops to forest
    for2ran = np.zeros((XSIZE, YSIZE, YEARS)) # crops to urban
    ran2for = np.zeros((XSIZE, YSIZE, YEARS)) # non-grass vegetation to grass
    for2pas = np.zeros((XSIZE, YSIZE, YEARS)) # non-grass vegetation to urban
    pas2for = np.zeros((XSIZE, YSIZE, YEARS)) # non-grass vegetation to crops
    for2urb = np.zeros((XSIZE, YSIZE, YEARS)) # grass to crops
    nfv2cro = np.zeros((XSIZE, YSIZE, YEARS)) # grass to non-grass vegetation
    cro2nfv = np.zeros((XSIZE, YSIZE, YEARS)) # grass to urban
    nfv2ran = np.zeros((XSIZE, YSIZE, YEARS)) # urban to grass
    ran2nfv = np.zeros((XSIZE, YSIZE, YEARS)) # urban to ucrops
    nfv2pas = np.zeros((XSIZE, YSIZE, YEARS)) # urban to non-grass vegetation
    pas2nfv = np.zeros((XSIZE, YSIZE, YEARS)) # forest to crops
    nfv2urb = np.zeros((XSIZE, YSIZE, YEARS)) # crops to forest
    ran2cro = np.zeros((XSIZE, YSIZE, YEARS)) # non-grass vegetation to grass
    ran2urb = np.zeros((XSIZE, YSIZE, YEARS)) # non-grass vegetation to urban
    pas2cro = np.zeros((XSIZE, YSIZE, YEARS)) # non-grass vegetation to crops
    cro2pas = np.zeros((XSIZE, YSIZE, YEARS)) # grass to crops
    pas2urb = np.zeros((XSIZE, YSIZE, YEARS)) # grass to non-grass vegetation
    cro2urb = np.zeros((XSIZE, YSIZE, YEARS)) # grass to urban
    ran2pas = np.zeros((XSIZE, YSIZE, YEARS)) # rangeland to pasture
    for2nfv = np.zeros((XSIZE, YSIZE, YEARS)) # forest to non-forest vegetation
    nfv2for = np.zeros((XSIZE, YSIZE, YEARS)) # non-forest to forest
    urb2for = np.zeros((XSIZE, YSIZE, YEARS)) # urban to forest
    urb2nfv = np.zeros((XSIZE, YSIZE, YEARS)) # urban to non-forest vegetation
    urb2cro = np.zeros((XSIZE, YSIZE, YEARS)) # urban to crops
    urb2pas = np.zeros((XSIZE, YSIZE, YEARS)) # urban to pasture
    urb2ran = np.zeros((XSIZE, YSIZE, YEARS)) # urban to rangeland
    nat2for = np.zeros((XSIZE, YSIZE, YEARS)) # natural vegetation to forest
    return for2cro, cro2for, for2ran, ran2for, for2pas, pas2for, \
           for2urb, nfv2cro, cro2nfv, nfv2ran, ran2nfv, nfv2pas, \
           pas2nfv, nfv2urb, ran2cro, ran2urb, pas2cro, cro2pas, \
           pas2urb, cro2urb, ran2pas, for2nfv, nfv2for, urb2for, \
           urb2nfv, urb2cro, urb2pas, urb2ran, nat2for

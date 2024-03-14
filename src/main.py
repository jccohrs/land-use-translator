from conf.mod_lucas_lut_iodata import return_iodata_arrays
from conf.mod_lucas_lut_trans import return_trans_arrays
from conf.conf import XSIZE, YSIZE, SYEAR, EYEAR


class LULCC:

    def __init__(self):
        # import arrays for input and output data
        x_size = int(input("Please enter the X size (default 443): ") or XSIZE)
        y_size = int(input("Please enter the Y size (default 443): ") or YSIZE)
        forward = bool(input("Please enter forward (default True): ") or True)
        years = EYEAR-SYEAR if forward else SYEAR-EYEAR
        rcm_lsm, pft_frac, pft_frac_ts, grass_backgr, crops_backgr, shrubs_backgr, \
        forest_backgr, urban_backgr, irri_frac, mcgrath_frac = return_iodata_arrays(XSIZE=x_size, YSIZE=y_size, YEARS=years)
        for2cro, cro2for, for2ran, ran2for, for2pas, pas2for, for2urb, nfv2cro, \
        cro2nfv, nfv2ran, ran2nfv, nfv2pas, pas2nfv, nfv2urb, ran2cro, ran2urb, \
        pas2cro, cro2pas, pas2urb, cro2urb, ran2pas, for2nfv, nfv2for, urb2for, \
        urb2nfv, urb2cro, urb2pas, urb2ran, nat2for = return_trans_arrays(XSIZE=x_size, YSIZE=y_size, YEARS=years)


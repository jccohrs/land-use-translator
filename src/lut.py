from pathlib import Path
import shutil
from cdo import *
from src.config import *
import os
import xarray as xr

cdo = Cdo()

class LUT:

    def __init__(self, config):
        for key, value in config.items():
            setattr(self, key, value)
        self.grid = f"reg01_{self.region}"
        self.glc = f"{config.lcd}-{config.esayear}-{config.vers}"
        self.glc_lsm = f"{config.lcd}-2015-{config.vers}"
        self.namelist = self.generate_namelist()
        if config.forward:
            self.years = config.eyear - config.syear
        else:
            self.years = config.syear - config.eyear
        self.nr_grass = nr_grass
        self.nr_crops = nr_crops
        self.nr_shrubs = nr_shrubs
        self.nr_forest = nr_forest
        self.nr_urban = nr_urban
        for i in range(0, 10):
            if CROPFTS[i] > 0:
                self.nr_crops += 1
            if FORPFTS[i] > 0:
                self.nr_forest += 1
            if GRAPFTS[i] > 0:
                self.nr_grass += 1
            if SHRPFTS[i] > 0:
                self.nr_shrubs += 1

        self.pfts_grass = GRAPFTS[0:self.nr_grass]
        self.pfts_crops = CROPFTS[0:self.nr_crops]
        self.pfts_shrubs = SHRPFTS[0:self.nr_shrubs]
        self.pfts_forest = FORPFTS[0:self.nr_forest]
        self.pfts_urban = URBPFTS[0:self.nr_urban]
        self.pfts_forest_shrubs = self.pfts_forest
        self.pfts_forest_shrubs_grass = self.pfts_forest
        self.pfts_shrubs_grass = self.pfts_shrubs
        self.pfts_forest_shrubs += self.pfts_shrubs
        self.pfts_forest_shrubs_grass += self.pfts_shrubs
        self.pfts_shrubs_grass += self.pfts_grass
        self.pfts_forest_shrubs_grass += self.pfts_grass
        self.pft_frac = xr.open_dataset(self.namelist["F_LC_IN"])

    def lucas_lut_forward(self):
        print(f"GRASS CLASSES: {self.pfts_grass}")
        print(f"CROP CLASSES: {self.pfts_crops}")
        print(f"SHRUB CLASSES: {self.pfts_shrubs}")
        print(f"FOREST CLASSES: {self.pfts_forest}")
        dummy = np.zeros(3)
        verb = False
        # list of pfts that consists of two or three PFT groups
        pfts_forest_shrubs = self.pfts_forest
        pfts_forest_shrubs_grass = self.pfts_forest
        pfts_shrubs_grass = self.pfts_shrubs
        pfts_forest_shrubs += self.pfts_shrubs
        pfts_forest_shrubs_grass += self.pfts_shrubs
        pfts_shrubs_grass += self.pfts_grass
        pfts_forest_shrubs_grass += self.pfts_grass
        # Dataset of F_BACKFOR
        xsize = self.namelist["XSIZE"]
        ysize = self.namelist["YSIZE"]
        forest_backgr = xr.open_dataset(self.namelist["F_BACKFOR"])
        shrubs_backgr = xr.open_dataset(self.namelist["F_BACKSHR"])
        grass_backgr = xr.open_dataset(self.namelist["F_BACKGRA"])
        third_dimension_size = max([self.nr_forest, self.nr_shrubs, self.nr_grass])
        forest_shrubs_backgr = np.arange(xsize*ysize*third_dimension_size).reshape(xsize, ysize, third_dimension_size)
        nat_backgr = np.arange(xsize*ysize*third_dimension_size).reshape(xsize, ysize, third_dimension_size)
        for i in range(xsize):
            for j in range(ysize):
                for ipft in range(self.nr_forest):
                    forest_shrubs_backgr[i, j, ipft] = forest_backgr[i, j, ipft] / 2.
                    nat_backgr[i, j, ipft] = forest_backgr[i, j, ipft] / 3.
                for ipft in range(self.nr_shrubs):
                    forest_shrubs_backgr[i, j, ipft + self.nr_forest] = shrubs_backgr[i, j, ipft] / 2.
                    shrubs_grass_backgr[i, j, ipft] = shrubs_backgr[i, j, ipft] / 2.
                    nat_backgr[i, j, ipft + self.nr_forest] = shrubs_backgr[i, j, ipft] / 3.
                for ipft in range(self.nr_grass):
                    shrubs_grass_backgr[i, j, ipft + self.nr_shrubs] = grass_backgr[i, j, ipft] / 2.
                    nat_backgr[i, j, ipft + self.nr_forest + self.nr_shrubs] = grass_backgr[i, j, ipft] / 3.
        forest_backgr.close()
        shrubs_backgr.close()
        grass_backgr.close()
        pft_frac_ts[:, :, :, 0] = self.pft_frac[:, :, :self.npfts]
        pft_help[:, :, :] = pft_frac[:, :, :]
        rcm_lsm = xr.open_dataset(self.namelist["F_RCM_LSM_IN"])
        for i in range(xsize):
            for j in range(ysize):
                for z in range(self.years):
                    if rcm_lsm[i, j] > 0.0:
                        # forest to crops
                        for2cro = xr.open_dataset(self.namelist["F_FOR2CRO"])
                        self.lucas_lut_transrules(for2cro[i, j, z], pft_help[i, j, :], self.pfts_crops, self.pfts_forest, self.pfts_shrubs, self.pfts_grass, self.nr_crops, self.nr_forest, self.nr_shrubs, self.nr_grass, pft_crops_default, crops_backgr, npfts, self.xsize, self.ysize, False, 3, False, dummy, verb)
                        for2cro.close()
                        # non-grass to crops
                        nfv2cro = xr.open_dataset(self.namelist["F_NFV2CRO"])
                        self.lucas_lut_transrules(nfv2cro[i, j, z], pft_help[i, j, :], self.pfts_crops, self.pfts_shrubs, self.pfts_grass, 0, self.nr_crops, self.nr_shrubs, self.nr_grass, 1, pft_crops_default, crops_backgr, self.npfts, self.xsize, self.ysize, False, 2, False, dummy, verb)
                        nfv2cro.close()
                        # rangeland to crops
                        ran2cro = xr.open_dataset(self.namelist["F_RAN2CRO"])
                        self.lucas_lut_transrules(ran2cro[i, j, z], pft_help[i, j, :], self.pfts_crops, self.pfts_shrubs, self.pfts_grass, 0, self.nr_crops, self.nr_shrubs, self.nr_grass, 1, pft_crops_default, crops_backgr, self.npfts, self.xsize, self.ysize, False, 2, False, dummy, verb)
                        ran2cro.close()
                        # pasture to crops
                        pas2cro = xr.open_dataset(self.namelist["F_PAS2CRO"])
                        self.lucas_lut_transrules(pas2cro[i, j, z], pft_help[i, j, :], self.pfts_crops, self.pfts_grass, 0, 0, self.nr_crops, self.nr_grass, 1, 1, pft_crops_default, crops_backgr, npfts, xsize, ysize, false, 1, false, dummy, verb)
                        pas2cro.close()
                        # crops to forest
                        cro2for = xr.open_dataset(self.namelist["F_CRO2FOR"])
                        self.lucas_lut_transrules(cro2for[i, j, z], pft_help[i, j, :], self.pfts_forest, self.pfts_crops, 0, 0, self.nr_forest, self.nr_crops, 1, 1, pft_forest_default, forest_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                        cro2for.close()
                        # crops to non-forest
                        cro2nfv = xr.open_dataset(self.namelist["F_CRO2NFV"])
                        self.lucas_lut_transrules(cro2nfv[i, j, z], pft_help[i, j, :], self.pfts_shrubs_grass, self.pfts_crops, 0, 0, self.nr_shrubs + self.nr_grass, self.nr_crops, 1, 1, pft_shrubs_default, shrubs_grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                        cro2nfv.close()
                        # crops to rangeland
                        cro2ran = xr.open_dataset(self.namelist["F_CRO2RAN"])
                        self.lucas_lut_transrules(cro2ran[i, j, z], pft_help[i, j, :], self.pfts_grass, self.pfts_crops, 0, 0, self.nr_grass, self.nr_crops, 1, 1, pft_grass_default, grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                        cro2ran.close()
                        # crops to rangeland
                        cro2pas = xr.open_dataset(self.namelist["F_CRO2PAS"])
                        self.lucas_lut_transrules(cro2pas[i, j, z], pft_help[i, j, :], self.pfts_grass, self.pfts_crops, 0, 0, self.nr_grass, self.nr_crops, 1, 1, pft_grass_default, grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                        cro2pas.close()
                        # crops to urban
                        cro2urb = xr.open_dataset(self.namelist["F_CRO2URB"])
                        self.lucas_lut_transrules(cro2urb[i, j, z], pft_help[i, j, :], self.pfts_urban, self.pfts_crops, 0, 0, self.nr_urban, self.nr_crops, 1, 1, pft_urban_default, urban_backgr, npfts, xsize, ysize, false, 1, false, dummy, verb)
                        cro2urb.close()
                        # forest to urban
                        for2urb = xr.open_dataset(self.namelist["F_FOR2URB"])
                        self.lucas_lut_transrules(for2urb[i, j, z], pft_help[i, j, :], self.pfts_urban, self.pfts_forest, self.pfts_shrubs, self.pfts_grass, self.nr_urban, self.nr_forest, nr_shrubs, nr_grass, pft_urban_default, urban_backgr, xsize, ysize, npfts, false, 3, false, dummy, verb)
                        for2urb.close()
                        # non-forest to urban
                        nfv2urb = xr.open_dataset(self.namelist["F_NFV2URB"])
                        self.lucas_lut_transrules(nfv2urb[i, j, z], pft_help[i, j, :], self.pfts_urban, self.pfts_shrubs, self.pfts_grass, 0, self.nr_urban, self.nr_shrubs, self.nr_grass, 1, pft_urban_default, npfts, xsize, ysize, backgrd, 2, false, dummy, verb)
                        nfv2urb.close()
                        # rangeland to urban
                        ran2urb = xr.open_dataset(self.namelist["F_RAN2URB"])
                        self.lucas_lut_transrules(ran2urb[i, j, z], pft_help[i, j, :], self.pfts_urban, self.pfts_shrubs_grass, 0, 0, self.nr_urban, self.nr_shrubs + self.nr_grass, 1, 1, pft_urban_default, urban_backgr[i, j, :], npfts, xsize, ysize, false, 1, false, dummy, verb)
                        ran2urb.close()
                        # pasture to urban
                        pas2urb = xr.open_dataset(self.namelist["F_PAS2URB"])
                        self.lucas_lut_transrules(pas2urb[i, j, z], pft_help[i, j, :], self.pfts_urban, self.pfts_grass, 0, 0, self.nr_urban, self.nr_grass, 1, 1, pft_urban_default, urban_backgr[i, j, :], npfts, xsize, ysize, false, 1, false, dummy, verb)
                        pas2urb.close()
                        # urban to crops
                        urb2cro = xr.open_dataset(self.namelist["F_URB2CRO"])
                        self.lucas_lut_transrules(urb2cro[i, j, z], pft_help[i, j, :], self.pfts_crops, self.pfts_urban, 0, 0, self.nr_crops, self.nr_urban, 1, 1, pft_crops_default, crops_backgr[i, j, :], npfts, xsize, ysize, false, 1, false, dummy, verb)
                        urb2cro.close()
                        # urban to not-forest
                        urb2nfv = xr.open_dataset(self.namelist["F_URB2NFV"])
                        self.lucas_lut_transrules(urb2nfv[i, j, z], pft_help[i, j, :], self.pfts_shrubs_grass, self.pfts_urban, 0, 0, self.nr_shrubs + self.nr_grass, self.nr_urban, 1, 1, pft_shrubs_default, shrubs_grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                        urb2nfv.close()
                        # urban to forest
                        urb2for = xr.open_dataset(self.namelist["F_URB2FOR"])
                        self.lucas_lut_transrules(urb2for[i, j, z], pft_help[i, j, :], self.pfts_forest, self.pfts_urban, 0, 0, self.nr_forest, self.nr_urban, 1, 1, pft_forest_default, forest_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                        urb2for.close()
                        # urban to rangeland
                        urb2ran = xr.open_dataset(self.namelist["F_URB2RAN"])
                        self.lucas_lut_transrules(urb2ran[i, j, z], pft_help[i, j, :], self.pfts_shrubs_grass, self.pfts_urban, 0, 0, self.nr_shrubs + self.nr_grass, self.nr_urban, 1, 1, pft_shrubs_default, shrubs_grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                        urb2ran.close()
                        # urban to pasture
                        urb2pas = xr.open_dataset(self.namelist["F_URB2PAS"])
                        self.lucas_lut_transrules(urb2pas[i, j, z], pft_help[i, j, :], self.pfts_grass, self.pfts_urban, 0, 0, self.nr_grass, self.nr_urban, 1, 1, pft_grass_default, grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                        urb2pas.close()
                        # forest to pasture
                        for2pas = xr.open_dataset(self.namelist["F_FOR2PAS"])
                        self.lucas_lut_transrules(for2pas[i, j, z], pft_help[i, j, :], self.pfts_grass, self.pfts_forest, 0, 0, self.nr_grass, self.nr_forest, 1, 1, pft_grass_default, grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                        for2pas.close()
                        # non-forest to pasture
                        nfv2pas = xr.open_dataset(self.namelist["F_NFV2PAS"])
                        self.lucas_lut_transrules(nfv2pas[i, j, z], pft_help[i, j, :], self.pfts_grass, self.pfts_shrubs, 0, 0, self.nr_grass, self.nr_shrubs, 1, 1, pft_grass_default, grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                        nfv2pas.close()
                        # rangeland to pasture
                        ran2pas = xr.open_dataset(self.namelist["F_RAN2PAS"])
                        self.lucas_lut_transrules(ran2pas[i, j, z], pft_help[i, j, :], self.pfts_grass, self.pfts_shrubs, 0, 0, self.nr_grass, self.nr_shrubs, 1, 1, pft_grass_default, grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                        ran2pas.close()
                        # pasture to forest
                        pas2for = xr.open_dataset(self.namelist["F_PAS2FOR"])
                        self.lucas_lut_transrules(pas2for[i, j, z], pft_help[i, j, :], self.pfts_forest, self.pfts_grass, 0, 0, self.nr_forest, self.nr_grass, 1, 1, pft_forest_default, forest_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                        pas2for.close()
                        # pasture to non-forest
                        pas2nfv = xr.open_dataset(self.namelist["F_PAS2NFV"])
                        self.lucas_lut_transrules(pas2nfv[i, j, z], pft_help[i, j, :], self.pfts_shrubs, self.pfts_grass, 0, 0, self.nr_shrubs, self.nr_grass, 1, 1, pft_shrubs_default, shrubs_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                        pas2nfv.close()
                        # forest to rangeland
                        for2ran = xr.open_dataset(self.namelist["F_FOR2RAN"])
                        self.lucas_lut_transrules(for2ran[i, j, z], pft_help[i, j, :], self.pfts_shrubs_grass, self.pfts_forest, 0, 0, self.nr_shrubs + self.nr_grass, self.nr_forest, 1, 1, pft_shrubs_default, shrubs_grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                        for2ran.close()
                        # non-forest to rangeland
                        nfv2ran = xr.open_dataset(self.namelist["F_NFV2RAN"])
                        self.lucas_lut_transrules(nfv2ran[i, j, z], pft_help[i, j, :], self.pfts_grass, self.pfts_shrubs, 0, 0, self.nr_grass, self.nr_shrubs, 1, 1, pft_grass_default, grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                        nfv2ran.close()
                        # rangeland to forest
                        ran2for = xr.open_dataset(self.namelist["F_RAN2FOR"])
                        self.lucas_lut_transrules(ran2for[i, j, z], pft_help[i, j, :], self.pfts_forest, self.pfts_shrubs, self.pfts_grass, 0, self.nr_forest, self.nr_shrubs, self.nr_grass, 1, pft_forest_default, forest_backgr[i, j, :], npfts, xsize, ysize, backgrd, 2, false, dummy, verb)
                        ran2for.close()
                        # forest to non-forest
                        for2nfv = xr.open_dataset(self.namelist["F_FOR2NFV"])
                        self.lucas_lut_transrules(for2nfv[i, j, z], pft_help[i, j, :], self.pfts_shrubs_grass, self.pfts_forest, 0, 0, self.nr_shrubs + self.nr_grass, self.nr_forest, 1, 1, pft_shrubs_default, shrubs_grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                        for2nfv.close()
                        # non-forest to forest
                        nfv2for = xr.open_dataset(self.namelist["F_NFV2FOR"])
                        self.lucas_lut_transrules(nfv2for[i, j, z], pft_help[i, j, :], self.pfts_forest, self.pfts_shrubs, self.pfts_grass, 0, self.nr_forest, self.nr_shrubs, self.nr_grass, 1, pft_forest_default, forest_backgr[i, j, :], npfts, xsize, ysize, backgrd, 2, false, dummy, verb)
                        nfv2for.close()
                    # add tree cover (only if addtree = true)
                    if self.addtree:
                        self.lucas_lut_transrules(nat2for[i, j, z], pft_help[i, j, :], self.pfts_forest, self.pfts_shrubs, self.pfts_grass, 0, self.nr_forest, self.nr_shrubs, self.nr_grass, 1, pft_forest_default, forest_backgr[i, j, :], npfts, xsize, ysize, backgrd, 2, false, dummy, verb)
                    pft_frac_ts[i, j, :, z + 1] = pft_help[i, j, :]
        print('land use change finished')
        for i in range(self.xsize):
            for j in range(self.ysize):
                for z in range(self.years):
                    if rcm_lsm[i][j] > 0.0:
                        pft_sum = 0.0
                        pft_sum = np.sum(pft_frac_ts[i][j][:][z])
                        #for ipft in range(self.npfts):
                        #    pft_sum += pft_frac_ts[i-1][j-1][ipft-1][z-1]
                        if pft_sum > 0.0:
                            #for ipft in range(self.npfts):
                            #    pft_frac_ts[i-1][j-1][ipft-1][z-1] /= pft_sum
                            pft_frac_ts[i][j][:][z] /= pft_sum
                    else:
                        #for ipft in range(self.npfts):
                        #    pft_frac_ts[i-1][j-1][ipft-1][z-1] = -999.
                        pft_frac_ts[i][j][:][z] = -999.

    def lucas_lut_transrules(self, trans, inpfts, pfts_1, pfts_2, pfts_3, pfts_4,
                             nr_pfts_1, nr_pfts_2, nr_pfts_3, nr_pfts_4, defaultpft,
                             backgrdpfts, npfts, xsize, ysize, backgrd, rule, mcgrath,
                             mcgfrac, verb):
        print(trans)
        pft_1_sum = 0.
        pft_2_sum = 0.
        pft_3_sum = 0.
        pft_4_sum = 0.
        mcg_sum = 0.
        for n in range(1, nr_pfts_1 + 1):
            pft_1_sum += inpfts[pfts_1[n-1]]
        for n in range(1, nr_pfts_2 + 1):
            pft_2_sum += inpfts[pfts_2[n-1]]
        if rule >= 2:
            for n in range(1, nr_pfts_3 + 1):
                pft_3_sum += inpfts[pfts_3[n-1]]
        if rule >= 3:
            for n in range(1, nr_pfts_4 + 1):
                pft_4_sum += inpfts[pfts_4[n-1]]
        if MCGRATH:
            for n in range(1, 4):
                mcg_sum += mcgfrac[n-1]
        helper = 0.0
        helper_2 = 0.0
        helper_3 = 0.0
        if trans > 0.0:
            if rule == 1:
                if verb:
                    print('BEGINN')
                    print(f"{str(pft_1_sum)}, {str(pft_2_sum)}, {str(pft_3_sum)}, {str(pft_3_sum)}")
                    print(trans)
                trans = min(trans, trans - max(0.0, pft_1_sum + trans - 1), trans - max(0.0, trans - pft_2_sum))
                if verb:
                    print(trans)
                    print(inpfts)
                if trans > 0.0:
                    if pft_2_sum > 0.0:
                        for ipft in range(NR_pfts_2):
                            inpfts[pfts_2[ipft]] = inpfts[pfts_2[ipft]] - (inpfts[pfts_2[ipft]] / pft_2_sum * trans)
                            if inpfts[pfts_2[ipft]] < 0.0:
                                helper = helper + inpfts[pfts_2[ipft]]
                                inpfts[pfts_2[ipft]] = 0.0
                    trans = max(0.0, trans + helper)
                    if mcgrath and mcg_sum == 1:
                        for ipft in range(3, 6):
                            inpfts[pfts_1[ipft]] = inpfts[pfts_1[ipft]] + (mcgfrac[ipft - 2] * trans)
                    else:
                        if pft_1_sum > 0.0:
                            for ipft in range(NR_pfts_1):
                                inpfts[pfts_1[ipft]] = inpfts[pfts_1[ipft]] + (inpfts[pfts_1[ipft]] / pft_1_sum * trans)
                        else:
                            if backgrd:
                                for ipft in range(NR_pfts_1):
                                    print(ipft, pfts_1[ipft], backgrdpfts[ipft])
                                    inpfts[pfts_1[ipft]] = inpfts[pfts_1[ipft]] + backgrdpfts[ipft] * trans
                            else:
                                inpfts[defaultpft] = trans
                if verb:
                    print(inpfts)
                    print('END')
            elif rule == 2:
                trans = min(trans, trans - max(0.0, pft_1_sum + trans - 1), trans - max(0.0, trans - (pft_2_sum + pft_3_sum)))
                if trans > 0:
                    if pft_2_sum > 0.0:
                        for ipft in range(nr_pfts_2):
                            inpfts[pfts_2[ipft]] = inpfts[pfts_2[ipft]] - (inpfts[pfts_2[ipft]] / pft_2_sum * trans)
                            if inpfts[pfts_2[ipft]] < 0.0:
                                helper = helper - inpfts[pfts_2[ipft]]
                                inpfts[pfts_2[ipft]] = 0.0
                    else:
                        helper = trans
                    if helper > 0.0:
                        if pft_3_sum > 0.0:
                            for ipft in range(nr_pfts_3):
                                inpfts[pfts_3[ipft]] = inpfts[pfts_3[ipft]] - (inpfts[pfts_3[ipft]] / pft_3_sum * helper)
                                if inpfts[pfts_3[ipft]] < 0.0:
                                    helper_2 = helper_2 + inpfts[pfts_3[ipft]]
                                    inpfts[pfts_3[ipft]] = 0.0
                    trans = max(0.0, trans + helper_2)
                    if mcgrath and mcg_sum == 1:
                        for ipft in range(3, 6):
                            inpfts[pfts_1[ipft]] = inpfts[pfts_1[ipft]] + (mcgfrac[ipft - 2] * trans)
                    else:
                        if pft_1_sum > 0.0:
                            for ipft in range(nr_pfts_1):
                                inpfts[pfts_1[ipft]] = inpfts[pfts_1[ipft]] + (inpfts[pfts_1[ipft]] / pft_1_sum * trans)
                        else:
                            if backgrd:
                                for ipft in range(nr_pfts_1):
                                    print(f"{str(ipft)}, {str(pfts_1[ipft])}, {str(backgrdpfts[ipft])}")
                                    inpfts[pfts_1[ipft]] = inpfts[pfts_1[ipft]] + backgrdpfts[ipft] * trans
                            else:
                                inpfts[defaultpft] = trans
            elif rule == 3:
                # limit the transition so that trans+pft1 is less equal 1 and trans-(pft2+pft3) is greater equal 0
                trans = min(trans, trans-max(0., pft_1_sum+trans-1), trans-max(0., trans-(pft_2_sum+pft_3_sum+pft_4_sum)))
                # subtracting from pft group 2
                if trans > 0.0:
                    if pft_2_sum > 0.0:
                        for ipft in range(1, nr_pfts_2+1):
                            inpfts[pfts_2[ipft-1]] -= (inpfts[pfts_2[ipft-1]]/pft_2_sum*trans)
                            # if more fraction is removed than available set fraction to zero
                            if inpfts[pfts_2[ipft-1]] < 0.0:
                                helper -= inpfts[pfts_2[ipft-1]]
                                inpfts[pfts_2[ipft-1]] = 0.0
                    else:
                        helper = trans
                    if helper > 0.0:
                        if pft_3_sum > 0.0:
                            for ipft in range(1, nr_pfts_3+1):
                                inpfts[pfts_3[ipft-1]] -= (inpfts[pfts_3[ipft-1]]/pft_3_sum*helper)
                                if inpfts[pfts_3[ipft-1]] < 0.0:
                                    helper_2 -= inpfts[pfts_3[ipft-1]]
                                    inpfts[pfts_3[ipft-1]] = 0.0
                        else:
                            helper_2 = helper
                    if helper_2 > 0.0:
                        if pft_4_sum > 0.0:
                            for ipft in range(1, nr_pfts_4+1):
                                inpfts[pfts_4[ipft-1]] -= (inpfts[pfts_4[ipft-1]]/pft_4_sum*helper_2)
                                if inpfts[pfts_4[ipft-1]] < 0.0:
                                    helper_3 += inpfts[pfts_4[ipft-1]]
                                    inpfts[pfts_4[ipft-1]] = 0.0
                # adjust transition if needed (not enough of pft group 2), helper is always <= 0
                trans = max(0., trans+helper_3)
                # add to pft group 2
                if mcgrath and (mcg_sum == 1):  # if mcgrath forest data should be used (hard coded for
                    for ipft in range(3, 6):
                        inpfts[pfts_1[ipft-1]] += (mcgfrac[ipft-3]*trans)
                else:  # just use the relative fractions
                    if pft_1_sum > 0.0:
                        for ipft in range(1, nr_pfts_1+1):
                            inpfts[pfts_1[ipft-1]] += (inpfts[pfts_1[ipft-1]]/pft_1_sum*trans)
                    else:
                        if backgrd:
                            for ipft in range(1, nr_pfts_1+1):
                                print(f"{str(ipft)}, {str(pfts_1[ipft-1])}, {str(backgrdpfts[ipft-1])}")
                                inpfts[pfts_1[ipft-1]] += backgrdpfts[ipft-1]*trans
                        else:
                            inpfts[defaultpft] = trans

    def lucas_lut_irrigation(self, pft_frac_ts, irri_frac, rcm_lsm):
        for i in range(self.xsize):
            for j in range(self.ysize):
                if rcm_lsm[i][j] > 0.0:
                    for z in range(years + 1):
                        sum_crops = pft_frac_ts[i][j][13][z] + pft_frac_ts[i][j][14][z]
                        if sum_crops > 0.0:
                            if irri_frac[i][j][z] > 0.0:
                                pft_frac_ts[i][j][13][z] = (1.0 - irri_frac[i][j][z]) * sum_crops
                                pft_frac_ts[i][j][14][z] = irri_frac[i][j][z] * sum_crops
                            else:
                                pft_frac_ts[i][j][13][z] = sum_crops
        for i in range(self.xsize):
            for j in range(self.ysize):
                for z in range(self.years + 1):
                    if rcm_lsm[i][j] > 0.0:
                        pft_sum = 0.0
                        for ipft in range(self.npfts):
                            if pft_frac_ts[i][j][ipft][z] < 0.0:
                                pft_frac_ts[i][j][ipft][z] = 0.0
                            pft_sum += pft_frac_ts[i][j][ipft][z]
                        if pft_sum > 0.0:
                            for ipft in range(self.npfts):
                                pft_frac_ts[i][j][ipft][z] /= pft_sum

    def lucas_lut_backward(self):
        dummy = np.zeros_like(pft_frac)
    
        forest_backgr = xr.open_dataset(self.namelist["F_BACKFOR"])
        shrubs_backgr = xr.open_dataset(self.namelist["F_BACKSHR"])

        forest_shrubs_backgr = np.zeros((xsize, ysize, nr_forest))
        nat_backgr = np.zeros((xsize, ysize, nr_forest))
        shrubs_grass_backgr = np.zeros((xsize, ysize, nr_shrubs))
        
        for i in range(xsize):
            for j in range(ysize):
                for ipft in range(nr_forest):
                    forest_shrubs_backgr[i, j, ipft] = forest_backgr[i, j, ipft] / 2.0
                    nat_backgr[i, j, ipft] = forest_backgr[i, j, ipft] / 3.0
                
                for ipft in range(nr_shrubs):
                    forest_shrubs_backgr[i, j, ipft + nr_forest] = shrubs_backgr[i, j, ipft] / 2.0
                    shrubs_grass_backgr[i, j, ipft] = shrubs_backgr[i, j, ipft] / 2.0
                    nat_backgr[i, j, ipft + nr_forest] = shrubs_backgr[i, j, ipft] / 3.0
        
        print('finished background map')
        
        pft_frac_ts = np.zeros((xsize, ysize, npfts, years + 1))
        pft_frac_ts[:, :, :, years] = pft_frac[:, :, :]
        
        pft_help = np.copy(pft_frac)
        
        for i in range(xsize):
            for j in range(ysize):
                for z in range(years):
                    zz = (years + 1) - z
                    if rcm_lsm[i, j] > 0.0:
                        # perform land cover changes
                        if i == 500 and j == 163:
                            print(zz)
                            print(pft_help[i, j, :])
                            VERB = True
                        else:
                            VERB = False
    
                    if i == 500 and j == 163:
                        print(nfv2cro[i, j, zz], PFT_help[i, j, 13], PFTS_SHRUBS, PFTS_CROPS, NR_SHRUBS, NR_CROPS)
    
                        self.lucas_lut_transrules(nfv2cro[i, j, zz], PFT_help[i, j], PFTS_SHRUBS, PFTS_CROPS, 0, 0, NR_SHRUBS, NR_CROPS, 1, 1, PFT_SHRUBS_DEFAULT, SHRUBS_BACKGR[i, j], NPFTS, XSIZE, YSIZE, BACKGRD, 1, False, DUMMY, VERB)
    
                        if i == 500 and j == 163:
                            print(PFT_help[i, j])
    
                        # Non-forest to crops
                        if i == 500 and j == 163:
                            print(cro2nfv[i, j, zz], PFT_help[i, j, 13])
    
                        self.lucas_lut_transrules(cro2nfv[i, j, zz], PFT_help[i, j], PFTS_CROPS, PFTS_SHRUBS, PFTS_GRASS, 0, NR_CROPS, NR_SHRUBS, NR_GRASS, 1, PFT_CROPS_DEFAULT, CROPS_BACKGR[i, j], NPFTS, XSIZE, YSIZE, BACKGRD, 2, False, DUMMY, VERB)
    
                        # Non-forest to crops
                        if i == 500 and j == 163:
                            print(for2cro[i, j, zz], PFT_help[i, j, 13])
    
                        self.lucas_lut_transrules(for2cro[i, j, zz], PFT_help[i, j], PFTS_FOREST, PFTS_CROPS, 0, 0, NR_FOREST, NR_CROPS, 1, 1, PFT_FOREST_DEFAULT, FOREST_BACKGR[i, j], NPFTS, XSIZE, YSIZE, BACKGRD, 1, MCGRATH, MCGRATH_FRAC[i, j, :, zz], VERB)
    
                        # Crops to non-forest
                        if i == 500 and j == 163:
                            print(cro2for[i, j, zz], PFT_help[i, j, 13])
    
                        self.lucas_lut_transrules(cro2for[i, j, zz], PFT_help[i, j], PFTS_CROPS, PFTS_FOREST, PFTS_SHRUBS, 0, NR_CROPS, NR_FOREST, NR_SHRUBS, 1, PFT_CROPS_DEFAULT, CROPS_BACKGR[i, j], NPFTS, XSIZE, YSIZE, False, 2, False, DUMMY, VERB)
    
                        # Crops to rangeland
                        if i == 500 and j == 163:
                            print(ran2cro[i, j, zz], PFT_help[i, j, 13])
    
                        self.lucas_lut_transrules(ran2cro[i, j, zz], PFT_help[i, j], PFTS_GRASS, PFTS_CROPS, 0, 0, NR_GRASS, NR_CROPS, 1, 1, PFT_GRASS_DEFAULT, GRASS_BACKGR[i, j], NPFTS, XSIZE, YSIZE, BACKGRD, 1, False, DUMMY, VERB)
    
                        # Crops to rangeland
                        if i == 500 and j == 163:
                            print(cro2ran[i, j, zz], PFT_help[i, j, 13])
    
                        self.lucas_lut_transrules(cro2ran[i, j, zz], PFT_help[i, j], PFTS_CROPS, PFTS_GRASS, 0, 0, NR_CROPS, NR_GRASS, 1, 1, PFT_CROPS_DEFAULT, CROPS_BACKGR[i, j], NPFTS, XSIZE, YSIZE, False, 1, False, DUMMY, VERB)

    def lucas_lut_mcgrath(years, x_size, y_size, npfts, pft_frac_ts, mcgrath_frac, rcm_lsm):
        for i in range(1, x_size + 1):
            for j in range(1, y_size + 1):
                if rcm_lsm[i-1][j-1] > 0.0:
                    for z in range(1, years + 1):
                        zz = years + 1 - z
                        d_mcg_frac = [0.0, 0.0, 0.0]
                        d_frac = [0.0, 0.0, 0.0]
                        d_forest = [0.0, 0.0, 0.0]
                        sum_forest = 0.0
                        sum_forest_p1 = 0.0
                        sum_helper = 0.0
                        sum_mcg = 0.0

                        for k in range(3, 6):
                            sum_forest += pft_frac_ts[i-1][j-1][k-1][zz-1]
                            sum_forest_p1 += pft_frac_ts[i-1][j-1][k-1][zz]
                            sum_mcg += mcgrath_frac[i-1][j-1][k-3][zz-1]

                        for k in range(3, 6):
                            if mcgrath_frac[i-1][j-1][k-3][zz-1] < -998.0 or mcgrath_frac[i-1][j-1][k-3][zz] < -998.0:
                                sum_mcg = 0.0

                        if sum_forest_p1 < 0.0 or sum_forest < 0.0:
                            pass  # Handle problem condition
                        elif sum_forest_p1 > 0.0 and sum_mcg > 0.0 and sum_forest > 0.0:
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

                                for k in range(3):
                                    sum_helper += lhelp[k]

                                if sum_helper < 3:
                                    sum_helper = 0.0

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

                            for k in range(3):
                                sum_helper += lhelp[k]

                            if sum_helper < 3:
                                sum_helper = 0.0

                                for k in range(3):
                                    if lhelp[k] == 1:
                                        sum_helper += pft_frac_ts[i-1][j-1][k+2][zz-1]

                                for k in range(3):
                                    if sum_helper > 0.0 and lhelp[k] == 1:
                                        pft_frac_ts[i-1][j-1][k+2][zz-1] += pft_frac_ts[i-1][j-1][k+2][zz-1] / sum_helper * helper
                                    else:
                                        pft_frac_ts[i-1][j-1][k+2][zz-1] -= (1 / 3) * helper

        for i in range(x_size):
            for j in range(y_size):
                for z in range(years + 1):
                    if rcm_lsm[i][j] > 0.0:
                        pft_sum = 0.0
                        for ipft in range(1, npfts + 1):
                            if pft_frac_ts[i][j][ipft-1][z-1] < 0.0:
                                pft_frac_ts[i][j][ipft-1][z-1] = 0.0
                            pft_sum += pft_frac_ts[i][j][ipft-1][z-1]

                        if pft_sum > 0.0:
                            for ipft in range(1, npfts + 1):
                                pft_frac_ts[i][j][ipft-1][z-1] = pft_frac_ts[i][j][ipft-1][z-1] / pft_sum

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
            reg = [-56, 84, 16, 79]
            if self.res == 100:
                xsize = 1400
                ysize = 630
            elif self.res == 250:
                xsize = 560
                ysize = 252
        elif self.region == "Global":
            reg = [-180, 180, -90, 90]
            xsize = 3600
            ysize = 1800
        elif self.region == "Australasia":
            reg = [102, 218, -53, 4]
            xsize = 1160
            ysize = 570
        elif self.region == "NorthAmerica":
            reg = [175, 360, 0, 85]
            xsize = 1900
            ysize = 850
        elif self.region == "GAR011":
            xsize = 145
            ysize = 129
        elif self.region == "Germany":
            if self.res ==  25:
                reg = [6, 15, 5, 46.4, 55, 5]
                xsize = 371
                ysize = 351
            elif self.res == 100:
                xsize = 95
                ysize = 89
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
        #PH need to include other options when using other LC dataset
        if self.mcgback and self.scenario in ["historical", "historical_high", "historical_low"]:
            ofile = f"{ofile}_mcg2"
        if self.forward:
            lutsyear = self.syear
            luteyear = self.eyear
        else:
            lutsyear = self.eyear
            luteyear = self.syear

        if self.addtree and self.scenario not in ["historical", "historical_high", "historical_low"]:
            ofile = f"{ofile}_addtr"
        Path(luhdir).mkdir(parents=True, exist_ok=True)
        namelist_dict = {
            # FILES
            "F_RCM_LSM_IN": f"{lsmdir}/{self.glc_lsm}_LSM_{self.grid}.nc", # lsmfile
            "F_LC_IN": f"{pftdir}/PFTS_{self.glc}_CRU3_{self.grid}_v11.nc", # pftfile
            "F_BACKGRA": f"{pftdir}/GRAB_{self.glc_lsm}_CRU3_{self.grid}_v11.nc", # grabfile
            "F_BACKSHR": f"{pftdir}/SHRB_{self.glc_lsm}_CRU3_{self.grid}_v11.nc", # shrbfile
            "F_BACKFOR": f"{pftdir}/FORB_{self.glc_lsm}_CRU3_{self.grid}_v11.nc", # forbfile
            "F_BACKCRO": f"{pftdir}/CROB_{self.glc_lsm}_CRU3_{self.grid}_v11.nc", # crobfile
            "F_BACKURB": f"{pftdir}/URBB_{self.glc_lsm}_CRU3_{self.grid}_v11.nc", # urbbfile
            "F_MCGRATH": f"{mcgdir}/{mcg}_{self.syear}_{self.eyear}_ForestBckgrdMcGrath_{self.grid}.nc", # mcgfile
            "F_MCGRATH": f"{sdir}/irrigation_{self.syear}_{self.eyear}_{self.grid}.nc", # irrfile
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
            "SYEAR":  lutsyear,
            "EYEAR": luteyear,
            "NPFTS" : 16,
            "GRADEF" : 9,
            "CRODEF": 13,
            "SHRDEF": 8,
            "FORDEF": 4,
            "URBDEF": 15,
            "FORPFTS": "1 , 2, 3, 4, 5, 6, 0, 0, 0, 0",
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
        URB=["urban"] # urban
        GRA=["range", "pastr"] # grass
        C3C=["c3ann", "c3per", "c3nfx"] # crops C3
        C4C=["c4ann", "c4per"] # crops C4
        CRO=["c3ann", "c3per", "c3nfx", "c4ann", "c4per"] #crops
        PRI=["primf", "primn"] # primary vegetation
        FOR=["primf", "secdf"] # forest
        NFV=["primn", "secdn"] # non-forest natural vegetation
        NAT=["primn", "secdn", "primf", "secdf"] # natural vegetation
        PAS=["pastr"] # pasture
        RAN=["range"] # rangeland
        IRR=['irrig_c3ann', 'irrig_c4ann', 'irrig_c3per', 'irrig_c4per', 'irrig_c3nfx'] # irrigation fractions
        ICR=['c3ann', 'c4ann', 'c3per', 'c4per', 'c3nfx'] # corrisponding states for irrigation 
        vars_state="primf,primn,secdf,secdn,urban,c3ann,c4ann,c3per,c4per,c3nfx,pastr,range,secmb,secma"
        vars_irrig="irrig_c3ann,irrig_c4ann,irrig_c3per,irrig_c4per,irrig_c3nfx"
        vars_crops="c3ann,c4ann,c3per,c4per,c3nfx"
        vars_trans="primf_to_secdn,primf_to_urban,primf_to_c3ann,primf_to_c4ann,primf_to_c3per,primf_to_c4per,primf_to_c3nfx,primf_to_pastr,primf_to_range,primn_to_secdf,primn_to_urban,primn_to_c3ann,primn_to_c4ann,primn_to_c3per,primn_to_c4per,primn_to_c3nfx,primn_to_pastr,primn_to_range,secdf_to_secdn,secdf_to_urban,secdf_to_c3ann,secdf_to_c4ann,secdf_to_c3per,secdf_to_c4per,secdf_to_c3nfx,secdf_to_pastr,secdf_to_range,secdn_to_secdf,secdn_to_urban,secdn_to_c3ann,secdn_to_c4ann,secdn_to_c3per,secdn_to_c4per,secdn_to_c3nfx,secdn_to_pastr,secdn_to_range,urban_to_secdf,urban_to_secdn,urban_to_c3ann,urban_to_c4ann,urban_to_c3per,urban_to_c4per,urban_to_c3nfx,urban_to_pastr,urban_to_range,c3ann_to_secdf,c3ann_to_secdn,c3ann_to_urban,c3ann_to_c4ann,c3ann_to_c3per,c3ann_to_c4per,c3ann_to_c3nfx,c3ann_to_pastr,c3ann_to_range,c4ann_to_secdf,c4ann_to_secdn,c4ann_to_urban,c4ann_to_c3ann,c4ann_to_c3per,c4ann_to_c4per,c4ann_to_c3nfx,c4ann_to_pastr,c4ann_to_range,c3per_to_secdf,c3per_to_secdn,c3per_to_urban,c3per_to_c3ann,c3per_to_c4ann,c3per_to_c4per,c3per_to_c3nfx,c3per_to_pastr,c3per_to_range,c4per_to_secdf,c4per_to_secdn,c4per_to_urban,c4per_to_c3ann,c4per_to_c4ann,c4per_to_c3per,c4per_to_c3nfx,c4per_to_pastr,c4per_to_range,c3nfx_to_secdf,c3nfx_to_secdn,c3nfx_to_urban,c3nfx_to_c3ann,c3nfx_to_c4ann,c3nfx_to_c3per,c3nfx_to_c4per,c3nfx_to_pastr,c3nfx_to_range,pastr_to_secdf,pastr_to_secdn,pastr_to_urban,pastr_to_c3ann,pastr_to_c4ann,pastr_to_c3per,pastr_to_c4per,pastr_to_c3nfx,pastr_to_range,range_to_secdf,range_to_secdn,range_to_urban,range_to_c3ann,range_to_c4ann,range_to_c3per,range_to_c4per,range_to_c3nfx,range_to_pastr,primf_harv,primn_harv,secmf_harv,secyf_harv,secnf_harv,primf_bioh,primn_bioh,secmf_bioh,secyf_bioh,secnf_bioh"

        scenario_dict = {
            "rcp19": "119",
            "rcp26": "126",
            "rcp34": "434",
            "rcp34OS": "534",
            "rcp45": "245",
            "rcp60": "460",
            "rcp70": "370",
            "rcp85": "585"
        }
        if self.scenario == "historical":
            sdir="historical"
            sfile="states"
            tfile="transitions"
            mfile="management"
        elif self.scenario == "historical_high":
            sdir="historical_high"
            sfile="states"
            tfile="transitions"
            mfile="management"
        elif self.scenario in scenario_dict.keys():
            sdir=f"scenarios/{self.scenario}"
            afile=f"added_tree_cover_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp{scenario_dict[self.scenario]}-2-1-f_gn_2015-2100"
            sfile=f"multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp{scenario_dict[self.scenario]}-2-1-f_gn_2015-2100"
            tfile=f"multiple-transitions_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp{scenario_dict[self.scenario]}-2-1-f_gn_2015-2100"
            mfile=f"multiple-management_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp{scenario_dict[self.scenario]}-2-1-f_gn_2015-2100"
        # here would come the mkdir for the sdir
        Path(luhdir).mkdir(parents=True, exist_ok=True)
        Path(os.path.join(luhdir, sdir)).mkdir(parents=True, exist_ok=True)

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

        # select period and self.region
        if self.region == 'Europe':
            reg = "-56,84,16,79"
        elif self.region == 'Australasia':
            reg = "102,218,-53,4"
        elif self.region == 'NorthAmerica':
            reg = "170,360,0,85"
        elif self.region == 'Germany':
            reg = "6,15.5,46.4,55.5"

        Path(os.path.join(luhdir, sdir, self.region)).mkdir(parents=True, exist_ok=True)
        Path(os.path.join(luhdir, sdir, self.region, self.grid)).mkdir(parents=True, exist_ok=True)
        path_region = os.path.join(luhdir, sdir, self.region)
        path_sdir = os.path.join(luhdir, sdir)

        # still have to fix this part to make it lighter. 
        if self.trans:
            ofile=f"{tfile}_{self.syear}_{self.eyear}_{self.region}.nc"
            cdo.sellonlatbox(reg, input=f"-selyear,{self.syear}/{self.eyear} -selvar,{vars_trans} {datadir}/{tfile}.nc", output=f"{path_sdir}/{self.region}/{ofile}")
        if self.state:
            #cdo.sellonlatbox(reg, input=f"-selyear,{self.syear}/{self.eyear} -selvar,{vars_state} {datadir}/{sfile}.nc", output=f"{path_region}/states_{self.syear}_{self.eyear}_{self.region}.nc")
            if remap_com == "invertlat":
                cdo.invertlat(input=f"{path_region}/states_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/states_{self.syear}_{self.eyear}_{self.grid}.nc")
            elif remap_com == "remapbil":
                cdo.remapbil(f"{scriptsdir}/grid_{self.grid}", input=f"{path_region}/states_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/states_{self.syear}_{self.eyear}_{self.grid}.nc")
            elif remap_com == "remapcon2":
                cdo.remapcon2(f"{scriptsdir}/grid_{self.grid}", input=f"{path_region}/states_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/states_{self.syear}_{self.eyear}_{self.grid}.nc")
        if self.addtree:
            cdo.sellonlatbox(reg, input=f"-selyear,{self.syear}/{self.eyear} -selvar,added_tree_cover {path_sdir}/{afile}.nc", output=f"{path_region}/addtree_{self.syear}_{self.eyear}_{self.region}.nc")
            if remap_com == "invertlat":
                cdo.invertlat(input=f"{path_region}/addtree_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/{self.grid}/addtree_{self.syear}_{self.eyear}_{self.grid}.nc")
            elif remap_com == "remapbil":
                cdo.remapbil(f"{scriptsdir}/grid_{self.grid}", input=f"{path_region}/addtree_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/{self.grid}/addtree_{self.syear}_{self.eyear}_{self.grid}.nc")
            elif remap_com == "remapcon2":
                cdo.remapcon2(f"{scriptsdir}/grid_{self.grid}", input=f"{path_region}/addtree_{self.syear}_{self.eyear}_{self.region}.nc", output=f"{path_region}/{self.grid}/addtree_{self.syear}_{self.eyear}_{self.grid}.nc")
            cdo.copy(input=f'{path_region}/{self.grid}/addtree_{self.syear}_{self.eyear}_{self.grid}.nc', output=f"{path_region}/{self.grid}/addtree_{self.syear}_{self.eyear}_{self.grid}.srv")

        # compute irragtion fraction 
        if self.irri:
            cdo.sellonlatbox(reg, input=f"-selyear,{self.syear}/{self.eyear} -selvar,{vars_irrig} {path_sdir}/{mfile}.nc", output=f"{path_region}/irri_vars_{self.syear}_{self.eyear}_{self.region}.nc")

            if self.scenario in ["historical", "historic_low", "historic_high"]:
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
        
        if self.trans:
            # New classification

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
                #{"varn": "pas2for", "for_1": PAS, "for_2": FOR, "outvar_condition": "primf"},
                {"varn": "for2pas", "for_1": FOR, "for_2": PAS, "outvar_condition": None},
                {"varn": "pas2cro", "for_1": PAS, "for_2": CRO, "outvar_condition": None},
                {"varn": "cro2pas", "for_1": CRO, "for_2": PAS, "outvar_condition": None},
                {"varn": "pas2urb", "for_1": PAS, "for_2": URB, "outvar_condition": None},
                {"varn": "ran2pas", "for_1": RAN, "for_2": PAS, "outvar_condition": None},
                {"varn": "urb2pas", "for_1": URB, "for_2": PAS, "outvar_condition": None}

            ]
            for data in fromto_array:
                self.fromto(data["varn"], data["for_1"], data["for_2"], tfile, ext, cutting, path_region, remap_com, data["outvar_condition"])

    def prepare_mcgrath(self):
        # commands for interpolation to given grid
        Path(glcdir).mkdir(parents=True, exist_ok=True)
        if self.remap == "bilinear":
            ext = "BIL"
            remap_com = f"remapbil, grid_{self.grid}"
            if self.grid == "EUR-011":
                cutting = "-selindexbox,2,434,2,434"
            else:
                cutting = ""
        else:
            remap_com = ""

        # select period and self.region
        if self.region == "Europe":
            reg = "-56,84,16,79"
        elif self.region == "Globe":
            reg = "-180,180,-90,90"
        else:
            reg = None
            # Should not be included the reg from Germany etc?
        ifile = f"{datadir}/{self.lcd}_{self.syear}_{self.mcgrath_eyear}.nc"
        #
        # compute background for LUT classes using zonal mean
        #
        cdo.chname(f"maxvegetfrac,{PFT_TeBrEv}", input=f"-vertsum -sellevel,{TeBrEv} {ifile}", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_TeBrEv.nc")
        cdo.chname(f"maxvegetfrac,{PFT_TeBrDec}", input=f"-vertsum -sellevel,{TeBrDec} {ifile}", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_TeBrDec.nc")
        cdo.chname(f"maxvegetfrac,{PFT_ConEv}", input=f"-vertsum -sellevel,{EvCon} {ifile}", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_ConEv.nc")
        cdo.chname(f"maxvegetfrac,forest", input=f"-vertsum -sellevel,{TeBrEv},{TeBrDec},{EvCon} {ifile}", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_FOR.nc")
        cdo.merge(input=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_TeBrEv.nc {glcdir}/{self.lcd}_{self.syear}_{self.eyear}_TeBrDec.nc {glcdir}/{self.lcd}_{self.syear}_{self.eyear}_ConEv.nc", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_dummy.nc")
        if reg:
            cdo.setmisstoc(-999, input=f"-remapbil,{scriptsdir}/grid_{self.grid} -sellonlatbox,{reg} -div {glcdir}/{self.lcd}_{self.syear}_{self.eyear}_dummy.nc -varssum {glcdir}/{self.lcd}_{self.syear}_{self.eyear}_dummy.nc", output=f"{glcdir}/{self.lcd}_{self.syear}_{self.eyear}_ForestBckgrdMcGrath_{self.grid}.nc")
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
        # combine land-use changes using reclassifcation\
        ifile=f"{tfile}_{self.syear}_{self.eyear}_{self.region}.nc"
        logfile=f"{tfile}_{self.syear}_{self.eyear}_{self.region}.log"
        ofile=f"transitions_{self.syear}_{self.eyear}_{self.region}_{varn}"
        cdo.mulc(0, input=f"-selvar,primf_to_urban {path}/{ifile}", output=f"{path}/dummy.nc")
        for inivar in for_1:
            for outvar in for_2:
                if not outvar_condition:
                    cdo.add(input=f"-selvar,{inivar}_to_{outvar} {path}/{ifile} {path}/dummy.nc", output=f"{path}/{ofile}.nc")
                    cdo.chname(f"{inivar}_to_{outvar},{varn}", input=f"{path}/{ofile}.nc", output=f" {path}/dummy.nc")
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

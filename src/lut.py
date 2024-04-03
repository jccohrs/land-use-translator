import numpy as np
from conf.conf import *

def lucas_lut_forward(years, xsize, ysize, npfts, backgrd, addtree):
    
    nr_grass = 0
    nr_crops = 0
    nr_shrubs = 0
    nr_forest = 0
    nr_urban = 1
    pfts_grass = list(filter(0, grapfts))
    pfts_crops = list(filter(0, cropfts))
    pfts_shrubs = list(filter(0, shrpfts))
    pfts_forest = list(filter(0, forpfts))
    pfts_urban = list(filter(0, urbpfts))

    rcm_lsm = np.zeros((xsize, ysize)) # input: labd-sea mask
    pft_frac = np.zeros((xsize, ysize, npfts)) # input: land cover fractions
    pft_frac_ts = np.zeros((xsize, ysize, npfts, years+1)) # output: time series of land cover fractions
    grass_backgr = np.zeros((xsize, ysize, nr_grass)) # background map for grass 
    crops_backgr = np.zeros((xsize, ysize, nr_crops)) # background map for crops
    shrubs_backgr = np.zeros((xsize, ysize, nr_shrubs)) # background map for shrubs
    forest_backgr = np.zeros((xsize, ysize, nr_forest)) # background map for forest
    urban_backgr = np.zeros((xsize, ysize, nr_urban)) # background map for forest

    verb = False
    pfts_forest_shrubs = [0] * nr_forest
    pfts_forest_shrubs_grass = [0] * nr_forest
    pfts_shrubs_grass = [0] * nr_shrubs
    for ipft in range(nr_forest):
        pfts_forest_shrubs[ipft] = pfts_forest[ipft]
        pfts_forest_shrubs_grass[ipft] = pfts_forest[ipft]
    for ipft in range(nr_shrubs):
        pfts_shrubs_grass[ipft] = pfts_shrubs[ipft]
        pfts_forest_shrubs.append(pfts_shrubs[ipft])
        pfts_forest_shrubs_grass.append(pfts_shrubs[ipft])
    for ipft in range(nr_grass):
        pfts_shrubs_grass.append(pfts_grass[ipft])
        pfts_forest_shrubs_grass.append(pfts_grass[ipft])
    forest_shrubs_backgr = [[[0 for _ in range(nr_forest + nr_shrubs)] for _ in range(ysize)] for _ in range(xsize)]
    shrubs_grass_backgr = [[[0 for _ in range(nr_shrubs + nr_grass)] for _ in range(ysize)] for _ in range(xsize)]
    nat_backgr = [[[0 for _ in range(nr_forest + nr_shrubs + nr_grass)] for _ in range(ysize)] for _ in range(xsize)]
    for i in range(xsize):
        for j in range(ysize):
            for ipft in range(nr_forest):
                forest_shrubs_backgr[i][j][ipft] = forest_backgr[i][j][ipft] / 2.
                nat_backgr[i][j][ipft] = forest_backgr[i][j][ipft] / 3.
            for ipft in range(nr_shrubs):
                forest_shrubs_backgr[i][j][ipft + nr_forest] = shrubs_backgr[i][j][ipft] / 2.
                shrubs_grass_backgr[i][j][ipft] = shrubs_backgr[i][j][ipft] / 2.
                nat_backgr[i][j][ipft + nr_forest] = shrubs_backgr[i][j][ipft] / 3.
            for ipft in range(nr_grass):
                shrubs_grass_backgr[i][j][ipft + nr_shrubs] = grass_backgr[i][j][ipft] / 2.
                nat_backgr[i][j][ipft + nr_forest + nr_shrubs] = grass_backgr[i][j][ipft] / 3.
    pft_frac_ts[:, :, :, 0] = pft_frac[:, :, :]
    pft_help[:, :, :] = pft_frac[:, :, :]
    for i in range(xsize):
        for j in range(ysize):
            for z in range(years):
                if rcm_lsm[i, j] > 0.0:
                    # forest to crops
                    lucas_lut_transrules(for2cro[i, j, z], pft_help[i, j, :], pfts_crops, pfts_forest, pfts_shrubs, pfts_grass, nr_crops, nr_forest, nr_shrubs, nr_grass, pft_crops_default, crops_backgr, npfts, xsize, ysize, false, 3, false, dummy, verb)
                    # non-grass to crops
                    lucas_lut_transrules(nfv2cro[i, j, z], pft_help[i, j, :], pfts_crops, pfts_shrubs, pfts_grass, 0, nr_crops, nr_shrubs, nr_grass, 1, pft_crops_default, crops_backgr, npfts, xsize, ysize, false, 2, false, dummy, verb)
                    # rangeland to crops
                    lucas_lut_transrules(ran2cro[i, j, z], pft_help[i, j, :], pfts_crops, pfts_shrubs, pfts_grass, 0, nr_crops, nr_shrubs, nr_grass, 1, pft_crops_default, crops_backgr, npfts, xsize, ysize, false, 2, false, dummy, verb)
                    # pasture to crops
                    lucas_lut_transrules(pas2cro[i, j, z], pft_help[i, j, :], pfts_crops, pfts_grass, 0, 0, nr_crops, nr_grass, 1, 1, pft_crops_default, crops_backgr, npfts, xsize, ysize, false, 1, false, dummy, verb)
                    # crops to forest
                    lucas_lut_transrules(cro2for[i, j, z], pft_help[i, j, :], pfts_forest, pfts_crops, 0, 0, nr_forest, nr_crops, 1, 1, pft_forest_default, forest_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                    # crops to non-forest
                    lucas_lut_transrules(cro2nfv[i, j, z], pft_help[i, j, :], pfts_shrubs_grass, pfts_crops, 0, 0, nr_shrubs + nr_grass, nr_crops, 1, 1, pft_shrubs_default, shrubs_grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                    # crops to rangeland
                    lucas_lut_transrules(cro2ran[i, j, z], pft_help[i, j, :], pfts_grass, pfts_crops, 0, 0, nr_grass, nr_crops, 1, 1, pft_grass_default, grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                    # crops to rangeland
                    lucas_lut_transrules(cro2pas[i, j, z], pft_help[i, j, :], pfts_grass, pfts_crops, 0, 0, nr_grass, nr_crops, 1, 1, pft_grass_default, grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                    # crops to urban
                    lucas_lut_transrules(cro2urb[i, j, z], pft_help[i, j, :], pfts_urban, pfts_crops, 0, 0, nr_urban, nr_crops, 1, 1, pft_urban_default, urban_backgr, npfts, xsize, ysize, false, 1, false, dummy, verb)
                    # forest to urban
                    lucas_lut_transrules(for2urb[i, j, z], pft_help[i, j, :], pfts_urban, pfts_forest, pfts_shrubs, pfts_grass, nr_urban, nr_forest, nr_shrubs, nr_grass, pft_urban_default, urban_backgr, xsize, ysize, npfts, false, 3, false, dummy, verb)
                    # non-forest to urban
                    lucas_lut_transrules(nfv2urb[i, j, z], pft_help[i, j, :], pfts_urban, pfts_shrubs, pfts_grass, 0, nr_urban, nr_shrubs, nr_grass, 1, pft_urban_default, npfts, xsize, ysize, backgrd, 2, false, dummy, verb)
                    # rangeland to urban
                    lucas_lut_transrules(ran2urb[i, j, z], pft_help[i, j, :], pfts_urban, pfts_shrubs_grass, 0, 0, nr_urban, nr_shrubs + nr_grass, 1, 1, pft_urban_default, urban_backgr[i, j, :], npfts, xsize, ysize, false, 1, false, dummy, verb)
                    # pasture to urban
                    lucas_lut_transrules(pas2urb[i, j, z], pft_help[i, j, :], pfts_urban, pfts_grass, 0, 0, nr_urban, nr_grass, 1, 1, pft_urban_default, urban_backgr[i, j, :], npfts, xsize, ysize, false, 1, false, dummy, verb)
                    # urban to crops
                    lucas_lut_transrules(urb2cro[i, j, z], pft_help[i, j, :], pfts_crops, pfts_urban, 0, 0, nr_crops, nr_urban, 1, 1, pft_crops_default, crops_backgr[i, j, :], npfts, xsize, ysize, false, 1, false, dummy, verb)
                    # urban to not-forest
                    lucas_lut_transrules(urb2nfv[i, j, z], pft_help[i, j, :], pfts_shrubs_grass, pfts_urban, 0, 0, nr_shrubs + nr_grass, nr_urban, 1, 1, pft_shrubs_default, shrubs_grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                    # urban to forest
                    lucas_lut_transrules(urb2for[i, j, z], pft_help[i, j, :], pfts_forest, pfts_urban, 0, 0, nr_forest, nr_urban, 1, 1, pft_forest_default, forest_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                    # urban to rangeland
                    lucas_lut_transrules(urb2ran[i, j, z], pft_help[i, j, :], pfts_shrubs_grass, pfts_urban, 0, 0, nr_shrubs + nr_grass, nr_urban, 1, 1, pft_shrubs_default, shrubs_grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                    # urban to pasture
                    lucas_lut_transrules(urb2pas[i, j, z], pft_help[i, j, :], pfts_grass, pfts_urban, 0, 0, nr_grass, nr_urban, 1, 1, pft_grass_default, grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                    # forest to pasture
                    lucas_lut_transrules(for2pas[i, j, z], pft_help[i, j, :], pfts_grass, pfts_forest, 0, 0, nr_grass, nr_forest, 1, 1, pft_grass_default, grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                    # non-forest to pasture
                    lucas_lut_transrules(nfv2pas[i, j, z], pft_help[i, j, :], pfts_grass, pfts_shrubs, 0, 0, nr_grass, nr_shrubs, 1, 1, pft_grass_default, grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                    # rangeland to pasture
                    lucas_lut_transrules(ran2pas[i, j, z], pft_help[i, j, :], pfts_grass, pfts_shrubs, 0, 0, nr_grass, nr_shrubs, 1, 1, pft_grass_default, grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                    # pasture to forest
                    lucas_lut_transrules(pas2for[i, j, z], pft_help[i, j, :], pfts_forest, pfts_grass, 0, 0, nr_forest, nr_grass, 1, 1, pft_forest_default, forest_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                    # pasture to non-forest
                    lucas_lut_transrules(pas2nfv[i, j, z], pft_help[i, j, :], pfts_shrubs, pfts_grass, 0, 0, nr_shrubs, nr_grass, 1, 1, pft_shrubs_default, shrubs_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                    # forest to rangeland
                    lucas_lut_transrules(for2ran[i, j, z], pft_help[i, j, :], pfts_shrubs_grass, pfts_forest, 0, 0, nr_shrubs + nr_grass, nr_forest, 1, 1, pft_shrubs_default, shrubs_grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                    # non-forest to rangeland
                    lucas_lut_transrules(nfv2ran[i, j, z], pft_help[i, j, :], pfts_grass, pfts_shrubs, 0, 0, nr_grass, nr_shrubs, 1, 1, pft_grass_default, grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                    # rangeland to forest
                    lucas_lut_transrules(ran2for[i, j, z], pft_help[i, j, :], pfts_forest, pfts_shrubs, pfts_grass, 0, nr_forest, nr_shrubs, nr_grass, 1, pft_forest_default, forest_backgr[i, j, :], npfts, xsize, ysize, backgrd, 2, false, dummy, verb)
                    # forest to non-forest
                    lucas_lut_transrules(for2nfv[i, j, z], pft_help[i, j, :], pfts_shrubs_grass, pfts_forest, 0, 0, nr_shrubs + nr_grass, nr_forest, 1, 1, pft_shrubs_default, shrubs_grass_backgr[i, j, :], npfts, xsize, ysize, backgrd, 1, false, dummy, verb)
                    # non-forest to forest
                    lucas_lut_transrules(nfv2for[i, j, z], pft_help[i, j, :], pfts_forest, pfts_shrubs, pfts_grass, 0, nr_forest, nr_shrubs, nr_grass, 1, pft_forest_default, forest_backgr[i, j, :], npfts, xsize, ysize, backgrd, 2, false, dummy, verb)
                # add tree cover (only if addtree = true)
                if addtree:
                    lucas_lut_transrules(nat2for[i, j, z], pft_help[i, j, :], pfts_forest, pfts_shrubs, pfts_grass, 0, nr_forest, nr_shrubs, nr_grass, 1, pft_forest_default, forest_backgr[i, j, :], npfts, xsize, ysize, backgrd, 2, false, dummy, verb)
                pft_frac_ts[i, j, :, z + 1] = pft_help[i, j, :]
    print('land use change finished')
    for i in range(1, xsize + 1):
        for j in range(1, ysize + 1):
            for z in range(1, years + 2):
                if rcm_lsm[i-1][j-1] > 0.0:
                    pft_sum = 0.0
                    for ipft in range(1, npfts + 1):
                        pft_sum += pft_frac_ts[i-1][j-1][ipft-1][z-1]
                    if pft_sum > 0.0:
                        for ipft in range(1, npfts + 1):
                            pft_frac_ts[i-1][j-1][ipft-1][z-1] /= pft_sum
                else:
                    for ipft in range(1, npfts + 1):
                        pft_frac_ts[i-1][j-1][ipft-1][z-1] = -999.
def lucas_lut_transrules(trans, inpfts, pfts_1, pfts_2, pfts_3, pfts_4,
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

    def lucas_lut_irrigation(years, xsize, ysize, npfts, pft_frac_ts, irri_frac, rcm_lsm):
        for i in range(xsize):
            for j in range(ysize):
                if rcm_lsm[i][j] > 0.0:
                    for z in range(years + 1):
                        sum_crops = pft_frac_ts[i][j][13][z] + pft_frac_ts[i][j][14][z]
                        if sum_crops > 0.0:
                            if irri_frac[i][j][z] > 0.0:
                                pft_frac_ts[i][j][13][z] = (1.0 - irri_frac[i][j][z]) * sum_crops
                                pft_frac_ts[i][j][14][z] = irri_frac[i][j][z] * sum_crops
                            else:
                                pft_frac_ts[i][j][13][z] = sum_crops
    
        for i in range(xsize):
            for j in range(ysize):
                for z in range(years + 1):
                    if rcm_lsm[i][j] > 0.0:
                        pft_sum = 0.0
                        for ipft in range(npfts):
                            if pft_frac_ts[i][j][ipft][z] < 0.0:
                                pft_frac_ts[i][j][ipft][z] = 0.0
                            pft_sum += pft_frac_ts[i][j][ipft][z]
                        if pft_sum > 0.0:
                            for ipft in range(npfts):
                                pft_frac_ts[i][j][ipft][z] /= pft_sum

#def lucas_lut_backward():
#    dummy = np.zeros_like(pft_frac)
#
#    pfts_forest_shrubs = np.zeros(nr_forest, dtype=int)
#    pfts_forest_shrubs_grass = np.zeros(nr_forest, dtype=int)
#    pfts_shrubs_grass = np.zeros(nr_shrubs, dtype=int)
#    
#    for ipft in range(nr_forest):
#        pfts_forest_shrubs[ipft] = pfts_forest[ipft]
#        pfts_forest_shrubs_grass[ipft] = pfts_forest[ipft]
#    
#    for ipft in range(nr_shrubs):
#        pfts_shrubs_grass[ipft] = pfts_shrubs[ipft]
#        pfts_forest_shrubs[ipft + nr_forest] = pfts_shrubs[ipft]
#        pfts_forest_shrubs_grass[ipft + nr_forest] = pfts_shrubs[ipft]
#    
#    for ipft in range(nr_grass):
#        pfts_shrubs_grass[ipft + nr_shrubs] = pfts_grass[ipft]
#        pfts_forest_shrubs_grass[ipft + nr_forest + nr_shrubs] = pfts_grass[ipft]
#    
#    print('finished pft groups')
#    
#    forest_shrubs_backgr = np.zeros((xsize, ysize, nr_forest))
#    nat_backgr = np.zeros((xsize, ysize, nr_forest))
#    shrubs_grass_backgr = np.zeros((xsize, ysize, nr_shrubs))
#    
#    for i in range(xsize):
#        for j in range(ysize):
#            for ipft in range(nr_forest):
#                forest_shrubs_backgr[i, j, ipft] = forest_backgr[i, j, ipft] / 2.0
#                nat_backgr[i, j, ipft] = forest_backgr[i, j, ipft] / 3.0
#            
#            for ipft in range(nr_shrubs):
#                forest_shrubs_backgr[i, j, ipft + nr_forest] = shrubs_backgr[i, j, ipft] / 2.0
#                shrubs_grass_backgr[i, j, ipft] = shrubs_backgr[i, j, ipft] / 2.0
#                nat_backgr[i, j, ipft + nr_forest] = shrubs_backgr[i, j, ipft] / 3.0
#    
#    print('finished background map')
#    
#    pft_frac_ts = np.zeros((xsize, ysize, npfts, years + 1))
#    pft_frac_ts[:, :, :, years] = pft_frac[:, :, :]
#    
#    pft_help = np.copy(pft_frac)
#    
#    for i in range(xsize):
#        for j in range(ysize):
#            for z in range(years):
#                zz = (years + 1) - z
#                if rcm_lsm[i, j] > 0.0:
#                    # perform land cover changes
#                    if i == 500 and j == 163:
#                        print(zz)
#                        print(pft_help[i, j, :])
#                        VERB = True
#                    else:
#                        VERB = False
#
#                if i == 500 and j == 163:
#                    print(nfv2cro[i, j, zz], PFT_help[i, j, 13], PFTS_SHRUBS, PFTS_CROPS, NR_SHRUBS, NR_CROPS)
#
#                    LUCAS_LUT_TRANSRULES(nfv2cro[i, j, zz], PFT_help[i, j], PFTS_SHRUBS, PFTS_CROPS, 0, 0, NR_SHRUBS, NR_CROPS, 1, 1, PFT_SHRUBS_DEFAULT, SHRUBS_BACKGR[i, j], NPFTS, XSIZE, YSIZE, BACKGRD, 1, False, DUMMY, VERB)
#
#                    if i == 500 and j == 163:
#                        print(PFT_help[i, j])
#
#                    # Non-forest to crops
#                    if i == 500 and j == 163:
#                        print(cro2nfv[i, j, zz], PFT_help[i, j, 13])
#
#                    LUCAS_LUT_TRANSRULES(cro2nfv[i, j, zz], PFT_help[i, j], PFTS_CROPS, PFTS_SHRUBS, PFTS_GRASS, 0, NR_CROPS, NR_SHRUBS, NR_GRASS, 1, PFT_CROPS_DEFAULT, CROPS_BACKGR[i, j], NPFTS, XSIZE, YSIZE, BACKGRD, 2, False, DUMMY, VERB)
#
#                    # Non-forest to crops
#                    if i == 500 and j == 163:
#                        print(for2cro[i, j, zz], PFT_help[i, j, 13])
#
#                    LUCAS_LUT_TRANSRULES(for2cro[i, j, zz], PFT_help[i, j], PFTS_FOREST, PFTS_CROPS, 0, 0, NR_FOREST, NR_CROPS, 1, 1, PFT_FOREST_DEFAULT, FOREST_BACKGR[i, j], NPFTS, XSIZE, YSIZE, BACKGRD, 1, MCGRATH, MCGRATH_FRAC[i, j, :, zz], VERB)
#
#                    # Crops to non-forest
#                    if i == 500 and j == 163:
#                        print(cro2for[i, j, zz], PFT_help[i, j, 13])
#
#                    LUCAS_LUT_TRANSRULES(cro2for[i, j, zz], PFT_help[i, j], PFTS_CROPS, PFTS_FOREST, PFTS_SHRUBS, 0, NR_CROPS, NR_FOREST, NR_SHRUBS, 1, PFT_CROPS_DEFAULT, CROPS_BACKGR[i, j], NPFTS, XSIZE, YSIZE, False, 2, False, DUMMY, VERB)
#
#                    # Crops to rangeland
#                    if i == 500 and j == 163:
#                        print(ran2cro[i, j, zz], PFT_help[i, j, 13])
#
#                    LUCAS_LUT_TRANSRULES(ran2cro[i, j, zz], PFT_help[i, j], PFTS_GRASS, PFTS_CROPS, 0, 0, NR_GRASS, NR_CROPS, 1, 1, PFT_GRASS_DEFAULT, GRASS_BACKGR[i, j], NPFTS, XSIZE, YSIZE, BACKGRD, 1, False, DUMMY, VERB)
#
#                    # Crops to rangeland
#                    if i == 500 and j == 163:
#                        print(cro2ran[i, j, zz], PFT_help[i, j, 13])
#
#                    LUCAS_LUT_TRANSRULES(cro2ran[i, j, zz], PFT_help[i, j], PFTS_CROPS, PFTS_GRASS, 0, 0, NR_CROPS, NR_GRASS, 1, 1, PFT_CROPS_DEFAULT, CROPS_BACKGR[i, j], NPFTS, XSIZE, YSIZE, False, 1, False, DUMMY, VERB)
#

# Function to read data from files
def read_data(file_name, years):
    data = {}
    for i in range(1, years+1):
        with open(file_name, "rb") as file:
            header = np.fromfile(file, dtype=np.float64, count=8)
            data[i] = np.fromfile(file, dtype=np.float64)  # Assuming data is stored as float64
    return data

def lucas_lut_input(years, npfts, f_rcm_lsm_in, f_lc_in, f_for2cro, f_cro2for, f_for2ran, f_ran2for, f_for2pas,
                    f_pas2for, f_for2urb, f_nfv2cro, f_cro2nfv, f_nfv2ran, f_ran2nfv, f_nfv2urb, f_ran2cro, f_cro2ran,
                    f_ran2urb, f_pas2cro, f_cro2pas, f_pas2urb, f_cro2urb, f_nfv2pas, f_pas2nfv, f_ran2pas, f_for2nfv,
                    f_nfv2for, f_urb2for, f_urb2nfv, f_urb2cro, f_urb2pas, f_urb2ran, f_lc_out, f_backcro, f_backshr,
                    f_backfor, f_backgra, f_backurb, f_irri_in, f_mcgrath, f_addtree, backgrd, irri, mcgrath, addtree):
    """
    -------------------------------------
        READ IN RCM LAND SEA MASK
    -------------------------------------
    """
    header = np.fromfile(file, dtype=np.float64, count=8)
    rcm_lsm = np.fromfile(file, dtype=np.float64)

    # RCM PFTS
    with open(f_lc_in, 'rb') as file:
        for i in range(npfts):
            header_8 = np.fromfile(file, dtype=np.dtype('int32'), count=8)
            pft_frac[:,:,i] = np.fromfile(file, dtype=np.dtype('float32'))
    file.close()

    # grass vegetation 
    
    # to crop
    for2cro = read_data(f_for2cro, years)

    # to urb
    for2urb = read_data(f_for2urb, years)

    # to grass
    for2ran = read_data(f_for2ran, years)
    
    
     non-grass vegetation
    cro2for = read_data(f_cro2for, years)

    # Crops to grass
    for2ran = read_data(f_for2ran, years)

    # Non-grass vegetation to grass
    ran2for = read_data(f_ran2for, years)

    # Grass to non-grass vegetation
    for2pas = read_data(f_for2pas, years)

    # Grass to crops
    pas2for = read_data(f_pas2for, years)

    # Crops to urban
    for2urb = read_data(f_for2urb, years)

    # Non-grass vegetation to urban
    nfv2urb = read_data(f_nfv2urb, years)


    # Non-grass vegetation to urban
    for2urb = read_data(f_for2urb, years)

    # Grass to non-grass vegetation
    nfv2cro = read_data(f_nfv2cro, years)

    # Grass to urban
    cro2nfv = read_data(f_cro2nfv, years)

    # Remaining conversions follow the same pattern

    # Read background maps if BACKGRD is true

    if backgrd:
        # Read grass background
        grass_backgr = read_data(f_backgra, nr_grass)  # NR_GRASS needs to be defined

        # Read crops background
        crops_backgr = read_data(f_backcro, nr_crops)  # NR_CROPS needs to be defined

        # Read shrubs background
        shrubs_backgr = read_data(f_backshr, nr_shrubs)  # NR_SHRUBS needs to be defined

        # Read forest background
        forest_backgr = read_data(f_backfor, nr_forest)  # NR_FOREST needs to be defined

        # Read urban background
        urban_backgr = read_data(f_backurb, nr_urban)  # NR_URBAN needs to be defined

    # Read irrigation fraction if IRRI is true

    if irri:
        irri_frac = read_data(f_irri_in, years+1)

    # Read MCGRATH relative fractions of forest types if MCGRATH is true

    if mcgrath:
        mcgrath_frac = {}
        for i in range(1, years+1):
            mcgrath_frac[i] = {}
            for j in range(1, 4):  # Assuming 3 forest types
                with open(f_mcgrath, "rb") as file:
                    header = np.fromfile(file, dtype=np.float64, count=8)
                    mcgrath_frac[i][j] = np.fromfile(file, dtype=np.float64)

    # Read additional tree data if ADDTREE is true

    if addtree:
        nat2for = read_data(f_addtree, years)


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

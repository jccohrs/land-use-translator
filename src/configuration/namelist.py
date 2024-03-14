from input_params import *
import json

def generate_namelist_dict():
    if res == 250:
        ext = "NINT"
    else:
        if remap == "bilinear":
            ext = "BIL"
        else:
            ext = ""
    
    # Select period and region
    if region == "Europe":
        reg = [-56, 84, 16, 79]
        if res == 100:
            xsize = 1400
            ysize = 630
        elif res == 250:
            xsize = 560
            ysize = 252
    elif region == "Global":
        reg = [-180, 180, -90, 90]
        xsize = 3600
        ysize = 1800
    elif region == "Australasia":
        reg = [102, 218, -53, 4]
        xsize = 1160
        ysize = 570
    elif region == "NorthAmerica":
        reg = [175, 360, 0, 85]
        xsize = 1900
        ysize = 850
    elif region == "GAR011":
        xsize = 145
        ysize = 129
    elif region == "Germany":
        if res ==  25:
            reg = [6, 15, 5, 46.4, 55, 5]
            xsize = 371
            ysize = 351
        elif res == 100:
            xsize = 95
            ysize = 89
        elif res == 250:
            xsize = 38
            ysize = 36
        elif res == 500:
            xsize = 19
            ysize = 18
    
    lmcg = False
    if scenario == "historical":
        sdir = f"{luhdir}/historic/{region}/{grid}"
        if mcgback:
            lmcg = True
    elif scenario == "historical_high":
        sdir = f"{luhdir}/historic_high/{region}/{grid}"
        if mcgback:
            lmcg = True
    elif scenario == "historical_low":
        sdir = f"{luhdir}/historic_low/{region}/{grid}"
        if mcgback:
            lmcg = True
    else:
        sdir = f"{luhdir}/scenarios/{scenario}/{region}/{grid}"

    if irri:
        ofile = f"{oname}_{scenario}_{syear}_{eyear}_{grid}_irri"
        lirri = True
    else:
        ofile = f"{oname}_{scenario}_{syear}_{eyear}_{grid}"
        lirri = False
    #PH need to include other options when using other LC dataset
    if mcgback and scenario in ["historical", "historical_high", "historical_low"]:
        ofile = f"{ofile}_mcg2"
    
    lforward = forward 
    lutsyear = syear
    luteyear = eyear
    lbackgrd = backgrd
    laddtr = addtr

    if addtr and scenario not in ["historical", "historical_high", "historical_low"]:
        ofile = f"{ofile}_addtr"

    namelist_dict = {
        # FILES
        "F_RCM_LSM_IN": f"{lsmdir}/{glc_lsm}_LSM_{grid}.srv", # lsmfile
        "F_LC_IN": f"{pftdir}/PFTS_{glc}_CRU3_{grid}_v11.srv", # pftfile
        "F_BACKGRA": f"{pftdir}/GRAB_{glc_lsm}_CRU3_{grid}_v11.srv", # grabfile
        "F_BACKSHR": f"{pftdir}/SHRB_{glc_lsm}_CRU3_{grid}_v11.srv", # shrbfile
        "F_BACKFOR": f"{pftdir}/FORB_{glc_lsm}_CRU3_{grid}_v11.srv", # forbfile
        "F_BACKCRO": f"{pftdir}/CROB_{glc_lsm}_CRU3_{grid}_v11.srv", # crobfile
        "F_BACKURB": f"{pftdir}/URBB_{glc_lsm}_CRU3_{grid}_v11.srv", # urbbfile
        "F_MCGRATH": f"{mcgdir}/{mcg}_{syear}_{eyear}_ForestBckgrdMcGrath_{grid}.srv", # mcgfile
        "F_MCGRATH": f"{sdir}/irrigation_{syear}_{eyear}_{grid}.srv", # irrfile
        "F_LC_OUT": f"{odir}/{ofile}.srv", # outfile
        "F_FOR2CRO": f"{sdir}/transitions_{syear}_{eyear}_{region}_for2cro_{ext}_{grid}.srv", # for2cro
        "F_CRO2FOR": f"{sdir}/transitions_{syear}_{eyear}_{region}_cro2for_{ext}_{grid}.srv", # cro2for
        "F_FOR2RAN": f"{sdir}/transitions_{syear}_{eyear}_{region}_for2ran_{ext}_{grid}.srv", # for2ran
        "F_RAN2FOR": f"{sdir}/transitions_{syear}_{eyear}_{region}_ran2for_{ext}_{grid}.srv", # ran2for
        "F_FOR2PAS": f"{sdir}/transitions_{syear}_{eyear}_{region}_for2pas_{ext}_{grid}.srv", # for2pas
        "F_PAS2FOR": f"{sdir}/transitions_{syear}_{eyear}_{region}_pas2for_{ext}_{grid}.srv", # pas2for
        "F_FOR2URB": f"{sdir}/transitions_{syear}_{eyear}_{region}_for2urb_{ext}_{grid}.srv", # for2urb
        "F_NFV2CRO": f"{sdir}/transitions_{syear}_{eyear}_{region}_nfv2cro_{ext}_{grid}.srv", # nfv2cro
        "F_CRO2NFV": f"{sdir}/transitions_{syear}_{eyear}_{region}_cro2nfv_{ext}_{grid}.srv", # cro2nfv
        "F_NFV2RAN": f"{sdir}/transitions_{syear}_{eyear}_{region}_nfv2ran_{ext}_{grid}.srv", # nfv2ran
        "F_RAN2NFV": f"{sdir}/transitions_{syear}_{eyear}_{region}_ran2nfv_{ext}_{grid}.srv", # ran2nfv
        "F_NFV2URB": f"{sdir}/transitions_{syear}_{eyear}_{region}_nfv2urb_{ext}_{grid}.srv", # nfv2urb
        "F_RAN2CRO": f"{sdir}/transitions_{syear}_{eyear}_{region}_ran2cro_{ext}_{grid}.srv", # ran2cro
        "F_CRO2RAN": f"{sdir}/transitions_{syear}_{eyear}_{region}_cro2ran_{ext}_{grid}.srv", # cro2ran
        "F_RAN2URB": f"{sdir}/transitions_{syear}_{eyear}_{region}_ran2urb_{ext}_{grid}.srv", # ran2urb
        "F_PAS2CRO": f"{sdir}/transitions_{syear}_{eyear}_{region}_pas2cro_{ext}_{grid}.srv", # pas2cro
        "F_CRO2PAS": f"{sdir}/transitions_{syear}_{eyear}_{region}_cro2pas_{ext}_{grid}.srv", # cro2pas
        "F_PAS2URB": f"{sdir}/transitions_{syear}_{eyear}_{region}_pas2urb_{ext}_{grid}.srv", # pas2urb
        "F_CRO2URB": f"{sdir}/transitions_{syear}_{eyear}_{region}_cro2urb_{ext}_{grid}.srv", # cro2urb
        "F_NFV2PAS": f"{sdir}/transitions_{syear}_{eyear}_{region}_nfv2pas_{ext}_{grid}.srv", # nfv2pas
        "F_PAS2NFV": f"{sdir}/transitions_{syear}_{eyear}_{region}_pas2nfv_{ext}_{grid}.srv", # pas2nfv
        "F_RAN2PAS": f"{sdir}/transitions_{syear}_{eyear}_{region}_ran2pas_{ext}_{grid}.srv", # ran2pas
        "F_FOR2NFV": f"{sdir}/transitions_{syear}_{eyear}_{region}_for2nfv_{ext}_{grid}.srv", # for2nfv
        "F_NFV2FOR": f"{sdir}/transitions_{syear}_{eyear}_{region}_nfv2for_{ext}_{grid}.srv", # nfv2for
        "F_URB2FOR": f"{sdir}/transitions_{syear}_{eyear}_{region}_urb2for_{ext}_{grid}.srv", # urb2for
        "F_URB2NFV": f"{sdir}/transitions_{syear}_{eyear}_{region}_urb2nfv_{ext}_{grid}.srv", # urb2nfv
        "F_URB2CRO": f"{sdir}/transitions_{syear}_{eyear}_{region}_urb2cro_{ext}_{grid}.srv", # urb2cro
        "F_URB2PAS": f"{sdir}/transitions_{syear}_{eyear}_{region}_urb2pas_{ext}_{grid}.srv", # urb2pas
        "F_URB2RAN": f"{sdir}/transitions_{syear}_{eyear}_{region}_urb2ran_{ext}_{grid}.srv", # urb2ran
        "F_ADDTREE": f"{sdir}/addtree_{syear}_{eyear}_{grid}.srv", # nat2for

        # CONTROL
        "XSIZE": xsize,
        "YSIZE": ysize,
        "MCGRATH": lmcg,
        "FORWARD": lforward,
        "BACKGRD": lbackgrd,
        "ADDTREE": laddtr,
        "IRRI": lirri,
        "SYEAR":  syear,
        "EYEAR": eyear,
        "NPFTS" : 16,
        "GRADEF" : 9,
        "CRODEF": 13,
        "SHRDEF": 8,
        "FORDEF": 4,
        "URBDEF": 15,
        "FORPFTS": [1 , 2, 3, 4, 5, 6, 0, 0, 0, 0],
        "SHRPFTS": [7, 8, 0, 0, 0, 0, 0, 0, 0, 0],
        "GRAPFTS": [9, 10, 11, 0, 0, 0, 0, 0, 0, 0],
        "CROPFTS": [13, 14, 0, 0, 0, 0, 0, 0, 0, 0],
        "URBPFTS": [15, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "CONPFTS": [12, 16, 0, 0, 0, 0, 0, 0, 0, 0]
    }
    return namelist_dict
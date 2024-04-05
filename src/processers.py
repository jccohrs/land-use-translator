from pathlib import Path
import shutil
from cdo import *
from src.config import *
import os

cdo = Cdo()

def generate_namelist(config):
    region = config.region
    res = config.res
    remap = config.remap
    scenario = config.scenario
    grid = f"reg01_{region}"
    irri = config.irri
    syear = config.syear
    eyear = config.eyear
    mcgback = config.mcgback
    forward = config.forward
    backgrd = config.backgrd
    addtr = config.addtree
    glc = f"{config.lcd}-{config.esayear}-{config.vers}"
    glc_lsm = f"{config.lcd}-2015-{config.vers}"
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
    else:
        ofile = f"{oname}_{scenario}_{syear}_{eyear}_{grid}"
    #PH need to include other options when using other LC dataset
    if mcgback and scenario in ["historical", "historical_high", "historical_low"]:
        ofile = f"{ofile}_mcg2"
    if forward:
        lutsyear = syear
        luteyear = eyear
    else:
        lutsyear = eyear
        luteyear = syear

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
        "FORWARD": forward,
        "BACKGRD": backgrd,
        "ADDTREE": addtr,
        "IRRI": irri,
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
    cdo.setgrid(os.path.join(scriptsdir, f"grid_{grid}"), input=f"{odir}/{ofile}.srv", output=f"{odir}/{ofile}.nc", options=f"-f nc -settime,12:00:00 -setctomiss,-999") 
    return namelist_dict

def prepare_luh2_data():
    """
    V2.0
    """
    URB="urban" # urban
    GRA="range pastr" # grass
    C3C="c3ann c3per c3nfx" # crops C3
    C4C="c4ann c4per" # crops C4
    CRO="c3ann c3per c3nfx c4ann c4per" #crops
    PRI="primf primn" # primary vegetation
    FOR="primf secdf" # forest
    NFV="primn secdn" # non-forest natural vegetation
    NAT="primn secdn primf secdf" # natural vegetation
    PAS="pastr" # pasture
    RAN="range" # rangeland
    IRR=('irrig_c3ann' 'irrig_c4ann' 'irrig_c3per' 'irrig_c4per' 'irrig_c3nfx') # irrigation fractions
    ICR=('c3ann' 'c4ann' 'c3per' 'c4per' 'c3nfx') # corrisponding states for irrigation 
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
    if scenario == "historic":
        sdir="historic"
        sfile="states"
        tfile="transitions"
        mfile="management"
    elif scenario == "historic_high":
        sdir="historic_high"
        sfile="states"
        tfile="transitions"
        mfile="management"
    elif scenario in scenario_dict.keys():
        sdir=f"scenarios/{scenario}"
        afile=f"added_tree_cover_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp{scenario_dict[scenario]}-2-1-f_gn_2015-2100"
        sfile=f"multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp{scenario_dict[scenario]}-2-1-f_gn_2015-2100"
        tfile=f"multiple-transitions_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp{scenario_dict[scenario]}-2-1-f_gn_2015-2100"
        mfile=f"multiple-management_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp{scenario_dict[scenario]}-2-1-f_gn_2015-2100"

    # here would come the mkdir for the sdir
    path = Path(luh2dir).mkdir(parents=True, exist_ok=True)
    Path(path / sdir).mkdir(parents=True, exist_ok=True)

    # interpolation
    if grid == "reg025_Europe":
        ext="NINT"
        remap_com="invertlat"
        cutting=''
    else:
        if remap == "bilinear":
            ext="BIL"
            remap_com=f"remapbil"
            if grid == "EUR-011":
                cutting = "-selindexbox,2,434,2,434"
            else:
                cutting = ""
        elif remap == "con2":
            ext = "CON2"
            remap_com = f"remapcon2"
            if grid == 'EUR-011':
                cutting = "-selindexbox,2,434,2,434"
            else:
                cutting = ''

    # select period and region
    if region == 'Europe':
        reg = "-56, 84, 16, 79"
    elif region == 'Australasia':
        reg = "102, 218, -53, 4"
    elif region == 'NorthAmerica':
        reg = "170, 360, 0, 85"
    elif region == 'Germany':
        reg = "6, 15.5, 46.4, 55.5"

    Path(f"{sdir}/{region}").mkdir(parents=True, exist_ok=True)
    shutil.copy(f"grid_{grid}", f"{sdir}/{region}")
    Path(f"{sdir}/{region}/{grid}").mkdir(parents=True, exist_ok=True)

    if trans:
        ofile=f"{tfile}_{syear}_{eyear}_{region}.nc"
        cdo.sellonlatbox(reg, input=f"-selyear,{syear}/{eyear} -selvar,{vars_trans} {sdir}/{tfile}.nc", output=f"{sdir}/{region}/{ofile}")
        # cdo sellonlatbox,$reg -selyear,$syear/$eyear -selvar,$vars_trans $sdir/$tfile.nc $sdir/$region/$ofile
    if state:
        #cdo sellonlatbox,$reg -selyear,$syear/$eyear -selvar,$vars_state $sdir/$sfile.nc $sdir/$region/states_${syear}_${eyear}_$region.nc
        cdo.sellonlatbox(reg, input=f"-selyear,{syear}/{eyear} -selvar,{vars_state} {sdir}/{sfile}.nc", output=f"{sdir}/{region}/states_{syear}_{eyear}_{region}.nc")
        #cdo $remap_com $sdir/$region/states_${syear}_${eyear}_$region.nc $sdir/$region/states_${syear}_${eyear}_$grid.nc
        if remap_com == "invertlat":
            cdo.invertlat(input=f"{sdir}/{region}/states_{syear}_{eyear}_{region}.nc", output=f"{sdir}/{region}/states_{syear}_{eyear}_{grid}.nc")
        elif remap_com == "remapbil":
            cdo.remapbil(f"grid_{grid}", input=f"{sdir}/{region}/states_{syear}_{eyear}_{region}.nc", output=f"{sdir}/{region}/states_{syear}_{eyear}_{grid}.nc")
        elif remap_com == "remapcon2":
            cdo.remapcon2(f"grid_{grid}", input=f"{sdir}/{region}/states_{syear}_{eyear}_{region}.nc", output=f"{sdir}/{region}/states_{syear}_{eyear}_{grid}.nc")
    if addtr:
        #cdo sellonlatbox,$reg -selyear,$syear/$eyear -selvar,added_tree_cover $sdir/$afile.nc $sdir/$region/addtree_${syear}_${eyear}_$region.nc
        #cdo $remap_com $sdir/$region/addtree_${syear}_${eyear}_$region.nc $sdir/$region/$grid/addtree_${syear}_${eyear}_$grid.nc
        #cdo -f srv copy $sdir/$region/$grid/addtree_${syear}_${eyear}_$grid.nc $sdir/$region/$grid/addtree_${syear}_${eyear}_$grid.srv
        cdo.sellonlatbox(reg, input=f"-selyear,{syear}/{eyear} -selvar,added_tree_cover {sdir}/{afile}.nc", output=f"{sdir}/{region}/addtree_{syear}_{eyear}_{region}.nc")
        if remap_com == "invertlat":
            cdo.invertlat(input=f"{sdir}/{region}/addtree_{syear}_{eyear}_{region}.nc", output=f"{sdir}/{region}/{grid}/addtree_{syear}_{eyear}_{grid}.nc")
        elif remap_com == "remapbil":
            cdo.remapbil(f"grid_{grid}", input=f"{sdir}/{region}/addtree_{syear}_{eyear}_{region}.nc", output=f"{sdir}/{region}/{grid}/addtree_{syear}_{eyear}_{grid}.nc")
        elif remap_com == "remapcon2":
            cdo.remapcon2(f"grid_{grid}", input=f"{sdir}/{region}/addtree_{syear}_{eyear}_{region}.nc", output=f"{sdir}/{region}/{grid}/addtree_{syear}_{eyear}_{grid}.nc")
        cdo.copy(input=f'{sdir}/{region}/{grid}/addtree_{syear}_{eyear}_{grid}.nc', output=f"{sdir}/{region}/{grid}/addtree_{syear}_{eyear}_{grid}.srv")
    
    # compute irragtion fraction 
    if irrig:
        cdo.sellonlatbox(reg, input=f"-selyear,{syear}/{eyear} -selvar,{vars_irrig} {sdir}/{mfile}.nc", output=f"{sdir}/{region}/irri_vars_{syear}_{eyear}_{region}.nc")
    
        if scenario in ["historic", "historic_low", "historic_high"]:
            if remap_com == "invertlat":
                cdo.invertlat(input=f"-chname,{IRR[0]},irrig_frac -selvar,{IRR[0]} {sdir}/{region}/irri_vars_{syear}_{eyear}_{region}.nc", output=f"{sdir}/{region}/{grid}/irrigation_{syear}_{eyear}_{grid}.nc")
            elif remap_com == "remapbil":
                cdo.remapbil(f"grid_{grid}", input=f"-chname,{IRR[0]},irrig_frac -selvar,{IRR[0]} {sdir}/{region}/irri_vars_{syear}_{eyear}_{region}.nc", output=f"{sdir}/{region}/{grid}/irrigation_{syear}_{eyear}_{grid}.nc")
            elif remap_com == "remapcon2":
                cdo.remapcon2(f"grid_{grid}", input=f"-chname,{IRR[0]},irrig_frac -selvar,{IRR[0]} {sdir}/{region}/irri_vars_{syear}_{eyear}_{region}.nc", output=f"{sdir}/{region}/{grid}/irrigation_{syear}_{eyear}_{grid}.nc")
            cdo.copy(options="-setmisstoc,-999", input=f'{sdir}/{region}/{grid}/irrigation_{syear}_{eyear}_{grid}.nc', output=f"{sdir}/{region}/{grid}/irrigation_{syear}_{eyear}_$grid.srv")
        else:
            # cdo setmisstoc,-999 -$remap_com -varssum -selvar,$vars_crops $sdir/$region/states_${syear}_${eyear}_$grid.nc $sdir/$region/sum_crop_frac.nc
            if remap_com in ["remapbil", "remapcon2"]:
                cdo.setmisstoc(input=f"-999 -{remap_com}, grid_{grid} -varssum -selvar,{vars_crops} {sdir}/{region}/states_{syear}_{eyear}_{grid}.nc", output=f"{sdir}/{region}/sum_crop_frac.nc")
            else:
                cdo.setmisstoc(input=f"-999 -{remap_com} -varssum -selvar,{vars_crops} {sdir}/{region}/states_{syear}_{eyear}_{grid}.nc", output=f"{sdir}/{region}/sum_crop_frac.nc")

            # Change variables names

            for n in range(5):
                if remap_com in ["remapbil", "remapcon2"]:
                    cdo.mul(input=f"-{remap_com}, grid_{grid} -selvar,{IRR[n]} {sdir}/{region}/irri_vars_{syear}_{eyear}_{region}.nc -{remap_com}, grid_{grid} -selvar,{ICR[n]} {sdir}/{region}/states_{syear}_{eyear}_{grid}.nc", output=f"{sdir}/{region}/dummy.nc")
                else:
                    cdo.mul(input=f"-{remap_com} -selvar,{IRR[n]} {sdir}/{region}/irri_vars_{syear}_{eyear}_{region}.nc -{remap_com} -selvar,{ICR[n]} {sdir}/{region}/states_{syear}_{eyear}_{grid}.nc", output=f"{sdir}/{region}/dummy.nc")
                if n == 0:
                    cdo.chname(input=f"{IRR[0]},irrig_frac -selvar,{IRR[0]} {sdir}/{region}/dummy.nc", output=f"{sdir}/{region}/sum_irri_frac_{syear}_{eyear}_{grid}.nc")
                else:
                    cdo.add(input=f"-chname,{IRR[n]},irrig_frac -selvar,{IRR[n]} {sdir}/{region}/dummy.nc {sdir}/{region}/sum_irri_frac_{syear}_{eyear}_{grid}.nc")   
                    shutil.move(f"{sdir}/{region}/dummy2.nc", f"{sdir}/{region}/sum_irri_frac_{syear}_{eyear}_{grid}.nc")
                os.remove(f"{sdir}/{region}/dummy.nc")
            cdo.div(input=f"{sdir}/{region}/sum_irri_frac_{syear}_{eyear}_{grid}.nc {sdir}/{region}/sum_crop_frac.nc", output=f"{sdir}/{region}/{grid}/irrigation_{syear}_{eyear}_{grid}.nc")
            cdo.copy(input=f"-setmisstoc,-999 {sdir}/{region}/{grid}/irrigation_{syear}_{eyear}_{grid}.nc", output=f"{sdir}/{region}/{grid}/irrigation_{syear}_{eyear}_{grid}.srv")
    if trans:

        # New classification

        fromto_array = [
            {"varn": for2urb, "for_1": FOR, "for_2": URB, "outvar_condition": None},
            {"varn": urb2for, "for_1": URB, "for_2": FOR, "outvar_condition": "primf"},
            {"varn": for2nfv, "for_1": FOR, "for_2": NFV, "outvar_condition": "primn"},
            {"varn": for2cro, "for_1": FOR, "for_2": CRO, "outvar_condition": None},
            {"varn": nfv2for, "for_1": NFV, "for_2": FOR, "outvar_condition": "primf"},
            {"varn": cro2for, "for_1": CRO, "for_2": FOR, "outvar_condition": "primf"},
            {"varn": cro2urb, "for_1": CRO, "for_2": URB, "outvar_condition": None},
            {"varn": urb2cro, "for_1": URB, "for_2": CRO, "outvar_condition": None},
            {"varn": cro2nfv, "for_1": CRO, "for_2": NFV, "outvar_condition": "primn"},
            {"varn": nfv2cro, "for_1": NFV, "for_2": CRO, "outvar_condition": None},
            {"varn": nfv2urb, "for_1": NFV, "for_2": URB, "outvar_condition": None},
            {"varn": urb2nfv, "for_1": URB, "for_2": NFV, "outvar_condition": None},
            {"varn": ran2nfv, "for_1": RAN, "for_2": NFV, "outvar_condition": "primn"},
            {"varn": nfv2ran, "for_1": NFV, "for_2": RAN, "outvar_condition": None},
            {"varn": ran2for, "for_1": RAN, "for_2": FOR, "outvar_condition": "primf"},
            {"varn": for2ran, "for_1": FOR, "for_2": RAN, "outvar_condition": None},
            {"varn": ran2cro, "for_1": RAN, "for_2": CRO, "outvar_condition": None},
            {"varn": cro2ran, "for_1": CRO, "for_2": RAN, "outvar_condition": None},
            {"varn": ran2urb, "for_1": RAN, "for_2": URB, "outvar_condition": None},
            {"varn": urb2ran, "for_1": URB, "for_2": RAN, "outvar_condition": None},
            {"varn": pas2nfv, "for_1": PAS, "for_2": NFV, "outvar_condition": "primn"},
            {"varn": nfv2pas, "for_1": NFV, "for_2": PAS, "outvar_condition": None},
            {"varn": pas2for, "for_1": PAS, "for_2": FOR, "outvar_condition": None},
            {"varn": for2pas, "for_1": FOR, "for_2": PAS, "outvar_condition": None},
            {"varn": pas2cro, "for_1": PAS, "for_2": CRO, "outvar_condition": None},
            {"varn": cro2pas, "for_1": CRO, "for_2": PAS, "outvar_condition": None},
            {"varn": pas2urb, "for_1": PAS, "for_2": URB, "outvar_condition": None},
            {"varn": ran2pas, "for_1": RAN, "for_2": PAS, "outvar_condition": None},
            {"varn": urb2pas, "for_1": URB, "for_2": PAS, "outvar_condition": None}

        ]

        for data in fromto_array:
            fromto(data["varn"], data["for_1"], data["for_2"], tfile, ext, cutting, data["outvar_condition"])

def prepare_mcgrath():
    # commands for interpolation to given grid
    if remap == "bilinear":
        ext = "BIL"
        remap_com = f"remapbil, grid_{grid}"
        if grid == "EUR-011":
            cutting = "-selindexbox,2,434,2,434"
        else:
            cutting = ""
    else:
        remap_com = ""
    
    # select period and region
    if region == "Europe":
        reg = "-56,84,16,79"
    elif region == "Globe":
        reg = "-180,180,-90,90"
    ifile = f"{glcdir}/{lcd}_{syear}_{eyear}.nc"
    #
    # compute background for LUT classes using zonal mean
    #
    cdo.chname(input=f"-vertsum -sellevel,{TeBrEv} {ifile}", output=f"{glcdir}/{lcd}_{syear}_{eyear}_TeBrEv.nc")
    cdo.maxvegetfrac(input=f"-vertsum -sellevel,{TeBrEv} {ifile}", output=f"{glcdir}/{lcd}_{syear}_{eyear}_TeBrEv.nc")
    cdo.var803(input=f"-vertsum -sellevel,{TeBrEv} {ifile}", output=f"{glcdir}/{lcd}_{syear}_{eyear}_TeBrEv.nc")
    cdo.chname(input=f"-vertsum -sellevel,{TeBrDec} {ifile}", output=f"{glcdir}/{lcd}_{syear}_{eyear}TeBrDec.nc")
    cdo.maxvegetfrac(input=f"-vertsum -sellevel,{TeBrDec} {ifile}", output=f"{glcdir}/{lcd}_{syear}_{eyear}TeBrDec.nc")
    cdo.var804(input=f"-vertsum -sellevel,{TeBrDec} {ifile}", output=f"{glcdir}/{lcd}_{syear}_{eyear}TeBrDec.nc")
    cdo.chname(input=f"-vertsum -sellevel,{EvCon} {ifile}", output=f"{glcdir}/{lcd}_{syear}_{eyear}_ConEv.nc")
    cdo.maxvegetfrac(input=f"-vertsum -sellevel,{EvCon} {ifile}", output=f"{glcdir}/{lcd}_{syear}_{eyear}_ConEv.nc")
    cdo.var805(input=f"-vertsum -sellevel,{EvCon} {ifile}", output=f"{glcdir}/{lcd}_{syear}_{eyear}_ConEv.nc")
    cdo.chname(input=f"-vertsum -sellevel,{TeBrEv},{TeBrDec},{EvCon} {ifile}", output=f"{glcdir}/{lcd}_{syear}_{eyear}_FOR.nc")
    cdo.maxvegetfrac(input=f"-vertsum -sellevel,{TeBrEv},{TeBrDec},{EvCon} {ifile}", output=f"{glcdir}/{lcd}_{syear}_{eyear}_FOR.nc")
    cdo.forest(input=f"-vertsum -sellevel,{TeBrEv},{TeBrDec},{EvCon} {ifile}", output=f"{glcdir}/{lcd}_{syear}_{eyear}_FOR.nc")
    cdo.merge(input=f"{glcdir}/{lcd}_{syear}_{eyear}_TeBrEv.nc {glcdir}/{lcd}_{syear}_{eyear}_TeBrDec.nc {glcdir}/{lcd}_{syear}_{eyear}_ConEv.nc", output=f"{glcdir}/{lcd}_{syear}_{eyear}_dummy.nc")
    cdo.setmisstoc(input=f"-999 -remapbil,grid_{grid} -sellonlatbox,{reg} -div {glcdir}/{lcd}_{syear}_{eyear}_dummy.nc -varssum {glcdir}/{lcd}_{syear}_{eyear}_dummy.nc", output=f"{glcdir}/{lcd}_{syear}_{eyear}_ForestBckgrdMcGrath_{grid}.nc")

    for year in [2011, 2012, 2013, 2014, 2015]:
        cdo.setdate(input=f"{year}-06-15 -selyear,2010 {glcdir}/{lcd}_{syear}_{eyear}_ForestBckgrdMcGrath_{grid}.nc", output=f"{glcdir}/dummy_{year}.nc")
    cdo.mergetime(input=f"{glcdir}/dummy_????.nc", output=f"{glcdir}/{lcd}_2011_2015_ForestBckgrdMcGrath_{grid}.nc")
    cdo.mergetime(input=f"{glcdir}/{lcd}_{syear}_{eyear}_ForestBckgrdMcGrath_{grid}.nc {glcdir}/{lcd}_2011_2015_ForestBckgrdMcGrath_{grid}.nc", output=f"{glcdir}/{lcd}_{syear}_2015_ForestBckgrdMcGrath_{grid}.nc")
    cdo.mergetime(input=f"{glcdir}/dummy_????.nc", output=f"{glcdir}/{lcd}_2011_2015_ForestBckgrdMcGrath_{grid}.nc")
    cdo.copy(input=f"{glcdir}/{lcd}_{syear}_2015_ForestBckgrdMcGrath_{grid}.nc", output=f"{glcdir}/{lcd}_{syear}_2015_ForestBckgrdMcGrath_{grid}.srv")

def fromto(varn, for_1, for_2, tfile, ext, cutting, outvar_condition=None):
    odir = grid
    # combine land-use changes using reclassifcation\
    ifile=f"{tfile}_{syear}_{eyear}_{region}.nc"
    logfile=f"{tfile}_{syear}_{eyear}_{region}.log"
    ofile=f"transitions_{syear}_{eyear}_{region}_{varn}"
    cdo.mulc(input=f"0 -selvar,primf_to_urban {ifile}", output=f"dummy.nc")
    for inivar in for_1:
        for outvar in for_2:
            if not outvar_condition:
                cdo.add(input=f"-selvar,{inivar}_to_{outvar} {ifile} dummy.nc", output=f"{ofile}.nc")
                cdo.chname(input=f"{inivar}_to_{outvar},{varn} {ofile}.nc", output=f"dummy.nc")
            if outvar != outvar_condition:
                print(f"{inivar}_to_{outvar}")
                cdo.add(input=f"-selvar,{inivar}_to_{outvar} {ifile} dummy.nc", output=f"{ofile}.nc")
                cdo.chname(input=f"{inivar}_to_{outvar},{varn} {ofile}.nc", output=f"dummy.nc")
    shutil.move(f"{sdir}/{region}/dummy.nc", f"{ofile}.nc")
    if remap_com == "invertlat":
        cdo.invertlat(input=f"{ofile}.nc", output=f"grid/{ofile}_{ext}_{grid}.nc")
    elif remap_com == "remapbil":
        cdo.remapbil(input=f"grid_{grid} {ofile}.nc", output=f"grid/{ofile}_{ext}_{grid}.nc")
    elif remap_com == "remapcon2":
        cdo.remapcon2(input=f"grid_{grid} {ofile}.nc", output=f"grid/{ofile}_{ext}_{grid}.nc")
    cdo.copy(options=f"{cutting} -setmisstoc,-999.", input=f'{grid}/{ofile}_{ext}_{grid}.nc', output=f"{grid}/{ofile}_{ext}_{grid}.srv")

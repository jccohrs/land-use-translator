from pathlib import Path
import shutil
from cdo import *
from src.config import *
import os

cdo = Cdo()

class ProcesserClass:

    def __init__(self, config):
        self.region = config.region
        self.res = config.res
        self.remap = config.remap
        self.scenario = config.scenario
        self.grid = f"reg01_{self.region}"
        self.irri = config.irri
        self.syear = config.syear
        self.eyear = config.eyear
        self.mcgback = config.mcgback
        self.forward = config.forward
        self.backgrd = config.backgrd
        self.addtree = config.addtree
        self.glc = f"{config.lcd}-{config.esayear}-{config.vers}"
        self.glc_lsm = f"{config.lcd}-2015-{config.vers}"
        self.lcd = config.lcd
        self.trans = config.trans
        self.state = config.state
        self.mcgrath_eyear = config.mcgrath_eyear

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

        namelist_dict = {
            # FILES
            "F_RCM_LSM_IN": f"{lsmdir}/{self.glc_lsm}_LSM_{self.grid}.srv", # lsmfile
            "F_LC_IN": f"{pftdir}/PFTS_{self.glc}_CRU3_{self.grid}_v11.srv", # pftfile
            "F_BACKGRA": f"{pftdir}/GRAB_{self.glc_lsm}_CRU3_{self.grid}_v11.srv", # grabfile
            "F_BACKSHR": f"{pftdir}/SHRB_{self.glc_lsm}_CRU3_{self.grid}_v11.srv", # shrbfile
            "F_BACKFOR": f"{pftdir}/FORB_{self.glc_lsm}_CRU3_{self.grid}_v11.srv", # forbfile
            "F_BACKCRO": f"{pftdir}/CROB_{self.glc_lsm}_CRU3_{self.grid}_v11.srv", # crobfile
            "F_BACKURB": f"{pftdir}/URBB_{self.glc_lsm}_CRU3_{self.grid}_v11.srv", # urbbfile
            "F_MCGRATH": f"{mcgdir}/{mcg}_{self.syear}_{self.eyear}_ForestBckgrdMcGrath_{self.grid}.srv", # mcgfile
            "F_MCGRATH": f"{sdir}/irrigation_{self.syear}_{self.eyear}_{self.grid}.srv", # irrfile
            "F_LC_OUT": f"{odir}/{ofile}.srv", # outfile
            "F_FOR2CRO": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_for2cro_{ext}_{self.grid}.srv", # for2cro
            "F_CRO2FOR": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_cro2for_{ext}_{self.grid}.srv", # cro2for
            "F_FOR2RAN": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_for2ran_{ext}_{self.grid}.srv", # for2ran
            "F_RAN2FOR": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_ran2for_{ext}_{self.grid}.srv", # ran2for
            "F_FOR2PAS": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_for2pas_{ext}_{self.grid}.srv", # for2pas
            "F_PAS2FOR": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_pas2for_{ext}_{self.grid}.srv", # pas2for
            "F_FOR2URB": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_for2urb_{ext}_{self.grid}.srv", # for2urb
            "F_NFV2CRO": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_nfv2cro_{ext}_{self.grid}.srv", # nfv2cro
            "F_CRO2NFV": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_cro2nfv_{ext}_{self.grid}.srv", # cro2nfv
            "F_NFV2RAN": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_nfv2ran_{ext}_{self.grid}.srv", # nfv2ran
            "F_RAN2NFV": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_ran2nfv_{ext}_{self.grid}.srv", # ran2nfv
            "F_NFV2URB": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_nfv2urb_{ext}_{self.grid}.srv", # nfv2urb
            "F_RAN2CRO": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_ran2cro_{ext}_{self.grid}.srv", # ran2cro
            "F_CRO2RAN": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_cro2ran_{ext}_{self.grid}.srv", # cro2ran
            "F_RAN2URB": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_ran2urb_{ext}_{self.grid}.srv", # ran2urb
            "F_PAS2CRO": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_pas2cro_{ext}_{self.grid}.srv", # pas2cro
            "F_CRO2PAS": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_cro2pas_{ext}_{self.grid}.srv", # cro2pas
            "F_PAS2URB": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_pas2urb_{ext}_{self.grid}.srv", # pas2urb
            "F_CRO2URB": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_cro2urb_{ext}_{self.grid}.srv", # cro2urb
            "F_NFV2PAS": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_nfv2pas_{ext}_{self.grid}.srv", # nfv2pas
            "F_PAS2NFV": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_pas2nfv_{ext}_{self.grid}.srv", # pas2nfv
            "F_RAN2PAS": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_ran2pas_{ext}_{self.grid}.srv", # ran2pas
            "F_FOR2NFV": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_for2nfv_{ext}_{self.grid}.srv", # for2nfv
            "F_NFV2FOR": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_nfv2for_{ext}_{self.grid}.srv", # nfv2for
            "F_URB2FOR": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_urb2for_{ext}_{self.grid}.srv", # urb2for
            "F_URB2NFV": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_urb2nfv_{ext}_{self.grid}.srv", # urb2nfv
            "F_URB2CRO": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_urb2cro_{ext}_{self.grid}.srv", # urb2cro
            "F_URB2PAS": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_urb2pas_{ext}_{self.grid}.srv", # urb2pas
            "F_URB2RAN": f"{sdir}/transitions_{self.syear}_{self.eyear}_{self.region}_urb2ran_{ext}_{self.grid}.srv", # urb2ran
            "F_ADDTREE": f"{sdir}/addtree_{self.syear}_{self.eyear}_{self.grid}.srv", # nat2for

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

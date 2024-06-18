pftdir = "data/ESA-CCI/PFT_FILES" # directory of ESA-CCI LC file
odir = "data/LUCAS_LUC" # directory that stores orginal GLOBCOVER files
scriptsdir = "scripts"
luhdir = "data/LUH_V2" # directory that stores output files from SAGA
mcgdir = "data/MCGRATH" # directory of McGrath forest fraction file
lut = "data/LUCAS_LUC/imove-preprocessing/lut/build/lucas_lut_levante.exe"
luh2dir = "/work/ch0636/g300089/LUCAS/DATA/LUH_V2"
datadir = "data"
plotdir = "plots"

URB = ["urban"] # urban
GRA = ["range", "pastr"] # grass
C3C = ["c3ann", "c3per", "c3nfx"] # crops C3
C4C = ["c4ann", "c4per"] # crops C4
CRO = ["c3ann", "c3per", "c3nfx", "c4ann", "c4per"] #crops
PRI = ["primf", "primn"] # primary vegetation
FOR = ["primf", "secdf"] # forest
NFV = ["primn", "secdn"] # non-forest natural vegetation
NAT = ["primn", "secdn", "primf", "secdf"] # natural vegetation
PAS = ["pastr"] # pasture
RAN = ["range"] # rangeland
IRR = ('irrig_c3ann' 'irrig_c4ann' 'irrig_c3per' 'irrig_c4per' 'irrig_c3nfx') # irrigation fractions
ICR = ('c3ann' 'c4ann' 'c3per' 'c4per' 'c3nfx') # corrisponding states for irrigation


glcdir = "data/MCGRATH" # directory of GLOBCOVER file
orgdir_mcgrath=f"{glcdir}/ORG_FILES" # directory that stores orginal GLOBCOVER files
sagadir_mcgrath=f"{glcdir}/SAGA_FILES" # directory that stores output files from SAGA
ncdir_mcgrath=f"{glcdir}/NC_FILES" # directory that stores netcdf output files
srvdir_mcgrath=f"{glcdir}/SRV_FILES" # directory that stores service binary files
pftdir_mcgrath=f"{glcdir}/PFT_FILES"
tdir=f"{glcdir}/tmp"

TeBrEv = 10
TeBrDec = "12,13,14,16,21"
EvCon = "5,6,7,8,18,19"
EvDec = 22
PFT_TeBrEv = "var803"
PFT_TeBrDec = "var804"
PFT_ConEv = "var805"
PFT_ConDec = "var806"


oname = "LUCAS_LUC7_ESACCI_LUH2"
mcg = "combined_species_mtc"

FORPFTS = [1, 2, 3, 4, 5, 6, 0, 0, 0, 0]
SHRPFTS = [7, 8, 0, 0, 0, 0, 0, 0, 0, 0]
GRAPFTS = [9, 10, 11, 0, 0, 0, 0, 0, 0, 0]
CROPFTS = [13, 14, 0, 0, 0, 0, 0, 0, 0, 0]
URBPFTS = [15, 0, 0, 0, 0, 0, 0, 0, 0, 0]
CONPFTS = [12, 16, 0, 0, 0, 0, 0, 0, 0, 0] 

nr_grass = 0
nr_crops = 0
nr_shrubs = 0
nr_forest = 0
nr_urban = 1

for i in range(0, 10):
    if CROPFTS[i] > 0:
        nr_crops += 1
    if FORPFTS[i] > 0:
        nr_forest += 1
    if GRAPFTS[i] > 0:
        nr_grass += 1
    if SHRPFTS[i] > 0:
        nr_shrubs += 1

# prepare LUH2 configuration
vars_state="primf,primn,secdf,secdn,urban,c3ann,c4ann,c3per,c4per,c3nfx,pastr,range,secmb,secma"
vars_irrig="irrig_c3ann,irrig_c4ann,irrig_c3per,irrig_c4per,irrig_c3nfx"
vars_crops="c3ann,c4ann,c3per,c4per,c3nfx"
vars_trans="primf_to_secdn,primf_to_urban,primf_to_c3ann,primf_to_c4ann,primf_to_c3per,primf_to_c4per,primf_to_c3nfx,primf_to_pastr,primf_to_range,primn_to_secdf,primn_to_urban,primn_to_c3ann,primn_to_c4ann,primn_to_c3per,primn_to_c4per,primn_to_c3nfx,primn_to_pastr,primn_to_range,secdf_to_secdn,secdf_to_urban,secdf_to_c3ann,secdf_to_c4ann,secdf_to_c3per,secdf_to_c4per,secdf_to_c3nfx,secdf_to_pastr,secdf_to_range,secdn_to_secdf,secdn_to_urban,secdn_to_c3ann,secdn_to_c4ann,secdn_to_c3per,secdn_to_c4per,secdn_to_c3nfx,secdn_to_pastr,secdn_to_range,urban_to_secdf,urban_to_secdn,urban_to_c3ann,urban_to_c4ann,urban_to_c3per,urban_to_c4per,urban_to_c3nfx,urban_to_pastr,urban_to_range,c3ann_to_secdf,c3ann_to_secdn,c3ann_to_urban,c3ann_to_c4ann,c3ann_to_c3per,c3ann_to_c4per,c3ann_to_c3nfx,c3ann_to_pastr,c3ann_to_range,c4ann_to_secdf,c4ann_to_secdn,c4ann_to_urban,c4ann_to_c3ann,c4ann_to_c3per,c4ann_to_c4per,c4ann_to_c3nfx,c4ann_to_pastr,c4ann_to_range,c3per_to_secdf,c3per_to_secdn,c3per_to_urban,c3per_to_c3ann,c3per_to_c4ann,c3per_to_c4per,c3per_to_c3nfx,c3per_to_pastr,c3per_to_range,c4per_to_secdf,c4per_to_secdn,c4per_to_urban,c4per_to_c3ann,c4per_to_c4ann,c4per_to_c3per,c4per_to_c3nfx,c4per_to_pastr,c4per_to_range,c3nfx_to_secdf,c3nfx_to_secdn,c3nfx_to_urban,c3nfx_to_c3ann,c3nfx_to_c4ann,c3nfx_to_c3per,c3nfx_to_c4per,c3nfx_to_pastr,c3nfx_to_range,pastr_to_secdf,pastr_to_secdn,pastr_to_urban,pastr_to_c3ann,pastr_to_c4ann,pastr_to_c3per,pastr_to_c4per,pastr_to_c3nfx,pastr_to_range,range_to_secdf,range_to_secdn,range_to_urban,range_to_c3ann,range_to_c4ann,range_to_c3per,range_to_c4per,range_to_c3nfx,range_to_pastr,primf_harv,primn_harv,secmf_harv,secyf_harv,secnf_harv,primf_bioh,primn_bioh,secmf_bioh,secyf_bioh,secnf_bioh"

scenario_dict = {
    "rcp19": "IMAGE-ssp119",
    "rcp26": "IMAGE-ssp126",
    "rcp34": "GCAM-ssp434",
    "rcp34OS": "MAGPIE-ssp534",
    "rcp45": "MESSAGE-ssp245",
    "rcp60": "GCAM-ssp460",
    "rcp70": "AIM-ssp370",
    "rcp85": "MAGPIE-ssp585"
}
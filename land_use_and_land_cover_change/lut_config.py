from utils import create_backgr_vars

pftdir = "data/ESA_CCI/PFT_FILES" # directory of ESA_CCI LC file
odir = "data/LUCAS_LUC" # directory that stores orginal GLOBCOVER files
scriptsdir = "scripts"
luhdir = "data/LUH_V2" # directory that stores output files from SAGA
mcgdir = "data/MCGRATH" # directory of McGrath forest fraction file
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
IRR = ['irrig_c3ann', 'irrig_c4ann', 'irrig_c3per', 'irrig_c4per', 'irrig_c3nfx'] # irrigation fractions
ICR = ['c3ann', 'c4ann', 'c3per', 'c4per', 'c3nfx'] # corrisponding states for irrigation


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
lcd = "LSM" # version of ESA-file
mcg = "combined_species_mtc" # name of McGrath file
vers = "v2.0.7" # version of the model

FORPFTS = [1, 2, 3, 4, 5, 6, 0, 0, 0, 0]
SHRPFTS = [7, 8, 0, 0, 0, 0, 0, 0, 0, 0]
GRAPFTS = [9, 10, 11, 0, 0, 0, 0, 0, 0, 0]
CROPFTS = [13, 14, 0, 0, 0, 0, 0, 0, 0, 0]
URBPFTS = [15, 0, 0, 0, 0, 0, 0, 0, 0, 0]
CONPFTS = [12, 16, 0, 0, 0, 0, 0, 0, 0, 0] 

GRADEF = 9
CRODEF = 13
SHRDEF = 8
FORDEF = 1
URBDEF = 15

nr_grass = 0
nr_crops = 0
nr_shrubs = 0
nr_forest = 0
nr_urban = 0

for i in range(10):
    if CROPFTS[i] > 0:
        nr_crops += 1
    if FORPFTS[i] > 0:
        nr_forest += 1
    if GRAPFTS[i] > 0:
        nr_grass += 1
    if SHRPFTS[i] > 0:
        nr_shrubs += 1
    if URBPFTS[i] > 0:
        nr_urban += 1

# prepare LUH2 configuration
vars_state="primf,primn,secdf,secdn,urban,c3ann,c4ann,c3per,c4per,c3nfx,pastr,range,secmb,secma"
vars_irrig="irrig_c3ann,irrig_c4ann,irrig_c3per,irrig_c4per,irrig_c3nfx"
vars_crops="c3ann,c4ann,c3per,c4per,c3nfx"
vars_trans="primf_to_secdn,primf_to_urban,primf_to_c3ann,primf_to_c4ann,primf_to_c3per,primf_to_c4per,primf_to_c3nfx,primf_to_pastr,primf_to_range,primn_to_secdf,primn_to_urban,primn_to_c3ann,primn_to_c4ann,primn_to_c3per,primn_to_c4per,primn_to_c3nfx,primn_to_pastr,primn_to_range,secdf_to_secdn,secdf_to_urban,secdf_to_c3ann,secdf_to_c4ann,secdf_to_c3per,secdf_to_c4per,secdf_to_c3nfx,secdf_to_pastr,secdf_to_range,secdn_to_secdf,secdn_to_urban,secdn_to_c3ann,secdn_to_c4ann,secdn_to_c3per,secdn_to_c4per,secdn_to_c3nfx,secdn_to_pastr,secdn_to_range,urban_to_secdf,urban_to_secdn,urban_to_c3ann,urban_to_c4ann,urban_to_c3per,urban_to_c4per,urban_to_c3nfx,urban_to_pastr,urban_to_range,c3ann_to_secdf,c3ann_to_secdn,c3ann_to_urban,c3ann_to_c4ann,c3ann_to_c3per,c3ann_to_c4per,c3ann_to_c3nfx,c3ann_to_pastr,c3ann_to_range,c4ann_to_secdf,c4ann_to_secdn,c4ann_to_urban,c4ann_to_c3ann,c4ann_to_c3per,c4ann_to_c4per,c4ann_to_c3nfx,c4ann_to_pastr,c4ann_to_range,c3per_to_secdf,c3per_to_secdn,c3per_to_urban,c3per_to_c3ann,c3per_to_c4ann,c3per_to_c4per,c3per_to_c3nfx,c3per_to_pastr,c3per_to_range,c4per_to_secdf,c4per_to_secdn,c4per_to_urban,c4per_to_c3ann,c4per_to_c4ann,c4per_to_c3per,c4per_to_c3nfx,c4per_to_pastr,c4per_to_range,c3nfx_to_secdf,c3nfx_to_secdn,c3nfx_to_urban,c3nfx_to_c3ann,c3nfx_to_c4ann,c3nfx_to_c3per,c3nfx_to_c4per,c3nfx_to_pastr,c3nfx_to_range,pastr_to_secdf,pastr_to_secdn,pastr_to_urban,pastr_to_c3ann,pastr_to_c4ann,pastr_to_c3per,pastr_to_c4per,pastr_to_c3nfx,pastr_to_range,range_to_secdf,range_to_secdn,range_to_urban,range_to_c3ann,range_to_c4ann,range_to_c3per,range_to_c4per,range_to_c3nfx,range_to_pastr,primf_harv,primn_harv,secmf_harv,secyf_harv,secnf_harv,primf_bioh,primn_bioh,secmf_bioh,secyf_bioh,secnf_bioh"
vars_cro = create_backgr_vars(nr_crops, CROPFTS[0])
vars_for = create_backgr_vars(nr_forest, FORPFTS[0])
vars_shr = create_backgr_vars(nr_shrubs, SHRPFTS[0])
vars_gra = create_backgr_vars(nr_grass, GRAPFTS[0])
vars_urb = create_backgr_vars(nr_urban, URBPFTS[0])
vars_pfts = ""
for i in range(1, 17):
    vars_pfts += "var"+str(800+i)+"," if i < 16 else "var"+str(800+i)

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

coords = {
    "Germany": "6.,15.5,46.5,55.4",
    "Europe": "-55.95,84.,16.1,79.",
    "WestAfrica": "-26.64,20.88,-1.52,28.18",
    "NorthAmerica": "170.,360.,0.,85",
    "AustralAsia":  "102.,218.,-53.,4."
}

th_file_syear = 850 # transitions historic file starts from 850
tf_file_syear = 2015 # transitions future file starts from 2015

def get_output_file_title(forward, syear, eyear):
    """
    Returns the title of the LUCAS LUC dataset, dynamically adjusting for historical or future datasets.
    """
    dataset_type = "historical" if not forward else "future"
    return f"Land Use and Land Cover (LULC) {dataset_type} ({syear}-{eyear}) dataset"

output_file_comment = """
Land Use and Land Cover (LULC) Dataset. 
Plant functional types and non-vegetated classes: 
1  - Tropical broadleaf evergreen trees 
2  - Tropical deciduous trees 
3  - Temperate broadleaf evergreen trees 
4  - Temperate deciduous trees 
5  - Evergreen coniferous trees 
6  - Deciduous coniferous trees 
7  - Coniferous shrubs 
8  - Deciduous shrubs 
9  - C3 grass 
10 - C4 grass 
11 - Tundra 
12 - Swamp 
13 - Non-irrigated crops 
14 - Irrigated crops 
15 - Urban 
16 - Bare
""".strip()  # Replace newlines with a space


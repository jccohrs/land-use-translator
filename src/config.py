pftdir = "data/ESA-CCI/PFT_FILES" # directory of ESA-CCI LC file
odir = "data/LUCAS_LUC" # directory that stores orginal GLOBCOVER files
scriptsdir = "scripts"
luhdir = "data/LUH_V2" # directory that stores output files from SAGA
mcgdir = "data/MCGRATH" # directory of McGrath forest fraction file
lsmdir = "data/ESA-CCI/SRV_FILES" # land-sea mask file
lut = "data/LUCAS_LUC/imove-preprocessing/lut/build/lucas_lut_levante.exe"
luh2dir = "/work/ch0636/g300089/LUCAS/DATA/LUH_V2"
datadir = "data"

URB = "urban" # urban
GRA = "range pastr" # grass
C3C = "c3ann c3per c3nfx" # crops C3
C4C = "c4ann c4per" # crops C4
CRO = "c3ann c3per c3nfx c4ann c4per" #crops
PRI = "primf primn" # primary vegetation
FOR = "primf secdf" # forest
NFV = "primn secdn" # non-forest natural vegetation
NAT = "primn secdn primf secdf" # natural vegetation
PAS = "pastr" # pasture
RAN = "range" # rangeland
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
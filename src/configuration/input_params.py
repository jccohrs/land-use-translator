# Input Parameter Generate

region = "NorthAmerica"
forward = True
mcgback = False
backgrd = True
irri = True
addtr = False

syear = 2015
eyear = 2100
esayear = 2015
oname = "LUCAS_LUC7_ESACCI_LUH2"
lcd = "ESACCI-LC-L4-LCCS-Map-300m-P1Y" # version of ESA-file
mcg = "combined_species_mtc" # name of McGrath file
vers = "v2.0.7"
lut = "../DATA/LUCAS_LUC/imove-preprocessing/lut/build/lucas_lut_levante.exe"
grid = f"reg01_{region}"
res = 100 # reso^lution in 10th of a degree
remap = "bilinear"
scenario = "rcp85"
pftdir = "../DATA/ESA-CCI/PFT_FILES" # directory of ESA-CCI LC file
odir = "../DATA/LUCAS_LUC" # directory that stores orginal GLOBCOVER files
luhdir = "../DATA/LUH_V2" # directory that stores output files from SAGA
mcgdir = "../DATA/MCGRATH" # directory of McGrath forest fraction file
lsmdir = "../DATA/ESA-CCI/SRV_FILES" # land-sea mask file


# Input Parameter Prepare

#primf: forested primary land - Forest
#primn: non-forested primary land - Natural
#secdf: potentially forested secondary land - secondary forest
#secdn: potentially non-forested secondary land - secondary natural
#pastr: managed pasture - Grass
#range: rangeland - Grass or Shrups
#urban: urban land - Urban
#c3ann: C3 annual crops - C3 crops
#c3per: C3 perennial crops - C3 crops
#c4ann: C4 annual crops - C4 crops
#c4per: C4 perennial crops - C4 crops
#c3nfx: C3 nitrogen-fixing crops - C3 crops

trans = True
state = False
luh2dir = "/work/ch0636/g300089/LUCAS/DATA/LUH_V2"

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

glc=f"{lcd}-{esayear}-{vers}"
glc_lsm=f"{lcd}-2015-{vers}"

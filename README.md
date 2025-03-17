 
# LUCAS Land Use Translator 


This Python program generates high-resolution land use and land cover (LULC) datasets for regional climate modeling across both historical and future periods using land use change information from the LUH2. The original project was developed in Fortran by Peter Hoffmann (see: [ESSD, 2023](https://essd.copernicus.org/articles/15/3819/2023/)) as part of HICSS project LANDMATE (https://hicss-hamburg.de/projects/landmate/index.php.en) and the CORDEX Flagship Pilot Study LUCAS (https://ms.hereon.de/cordex_fps_lucas/). By integrating this LULC data, climate models can more accurately simulate the impacts of land use changes on regional climate dynamics.

## Installation Guide

1. **Install pipenv**:

   ```bash
   pip install pipenv

2. **Change to the project directory**:

   ```bash
   cd project_name

3. **Install the required libraries**:

   ```bash
   pipenv install 

4. **Activate your pipenv environment**:
   ```bash
   pipenv shell

## Quick start (Example historical scenario 1950-2015 for Germany)

To get started with this project, follow these steps:

1. **Change to the projects data directory**:

   ```bash
   cd land_use_and_land_cover_change/data/
   ```
2. **Download required data files using `wget`:**

   ```bash
   wget --continue --progress=bar --no-check-certificate \
    "https://luh.umd.edu/LUH2/LUH2_v2h/transitions.nc" \
    "https://luh.umd.edu/LUH2/LUH2_v2h/states.nc" \
    "https://luh.umd.edu/LUH2/LUH2_v2h/management.nc" \
    "https://zenodo.org/records/14981619/files/CROB_reg01_Global.nc?download=1" \
    "https://zenodo.org/records/14981619/files/FORB_reg01_Global.nc?download=1" \
    "https://zenodo.org/records/14981619/files/GRAB_reg01_Global.nc?download=1" \
    "https://zenodo.org/records/14981619/files/SHRB_reg01_Global.nc?download=1" \
    "https://zenodo.org/records/14981619/files/URBB_reg01_Global.nc?download=1" \
    "https://zenodo.org/records/14981619/files/PFTS_reg01.nc?download=1"
   ```
These files will be downloaded and save in `data` directory.

3. **Run the program from the project directory:**

```bash
python3 main.py
```

5. **Output file:** The generated output file will be located in the  `data/LUCAS_LUC/` directory. This output file contains the Plant Functional Type (PFT) fraction for the 14 PFTs and two non-vegetated land cover classes across the selected region, scenario, and timeline.

## Dataset Requirements

The following datasets are required to run the program:

- **Land Use (LU) Transitions provided by LUH2**: LU changes for the selected time period (will be filtered in case of providing larger datasets). The following files can be downloaded from the [LUH Data Portal](https://luh.umd.edu/data.shtml) for different scenarios (historical, historical_high, rcp19 etc.):
   - transitions.nc: The land-use transitions are the annual changes between land-use states. 
   - states.nc: The land-use states are the fractions of each grid-cell occupied by various land-uses in a given year.
   - management.nc: Contains irrigation data; this file is only required if irrigation is enabled (i.e., if `irri` is set to True).
   - add_tree_cover.nc: Tree cover data; this file is only required if if addtree is enabled (i.e., if `addtree` is set to True).

These files should be then moved to `land_use_and_land_cover_change/data/`.
- **Landmate PFTs Maps**: Contains data for 14 Plant Functional Types (PFTs) and 2 non-vegetated land cover classes (urban and bare ground), providing detailed vegetation characterization. The data for Europe can be downloaded from [WDC Climate](https://www.wdc-climate.de/ui/entry?acronym=LM_PFT_EUR_v1.1_afts). **From [WDC Climate](https://www.wdc-climate.de/ui/entry?acronym=LM_PFT_EUR_v1.1_afts).** you can download the .nc file for Europe (you will have to log in first). This file should be then moved to `land_use_and_land_cover_change/data/` and renamed as `PFTS_reg01.nc`. 

- **Mcgrath Data (Optional)**: For the backward extension of historical forest type distribution, additional information on the relative distribution of broadleaf and needleleaf forests, derived from the McGrath dataset, can be utilized. For more information about obtaining this dataset, please contact the maintainers of this project.
- **Land-sea Mask (Optional)**: 
   - By default, the land-sea mask will be calculated from the Landmate PFT maps based on land classification.
   - If you want to use a custom land-sea mask, you can provide the path to the file via the `path_file_lsm` parameter in the configuration file.
- **Background Data (Optional)**: In certain cases, where a certain vegetation type is not present within a grid cell but should be increased according to the LUH2 and the rules provided by the LUT, a background map of potential vegetation is needed.
   - The project already provides global background data. This data can be used to enhance the simulation (i.e. set `backgrd` to True). If you prefer to use your own regional or global background data, you can specify the path to the new data files via `path_file_back*` and `path_file_back*_global` parameters.

By default, the datasets should be located in a designated `data/` directory. If different storage locations are prefered, alternative paths for each dataset can be specified in a configuration file. This will be further detailed in the Usage section bellow.

 

## Configuration
The main configuration file is located at `config/main.yaml`. Modify the parameters in this file to customize the program:

#### LUT configuration
   - (`region`): Choose from pre-configured regions ("Germany", "Europe", "WestAfrica", "NorthAmerica", "Australasia"), or add a new region by providing the necessary grid files (located into `scripts/`) and coordinates (`coords`parameter).
   - (`forward`): **True** for computation of future scenarios or **False** for historical reconstruction.
   - (`backgrd`): True/False. Optionally include background data.
   - (`mcgrath`): True/False for using mcgrath data in the LUT. 
   - (`addtree`): True/False for using addtree data in the LUT. 
   - (`irri`): True/False. Enable or disable irrigation data, if the irrigation dataset is available.
   - (`syear`)/(`eyear`): Specify the time period for LU calculations by setting the starting year (`syear`) and ending year (`eyear`).
   - (`mcgrath_eyear`): end year of mcgrath file (in case that its different from eyear).
   - (`npfts`): number of npfts used in the LUT. Currently, 16 is required but might be changed if different land cover dataset is used.
   - (`xsize`): xsize of the region.
   - (`ysize`): ysize of the region. 


   ### LUH2 prepare data configuration
   - `prepare_luh2_data`: True/False. Files preparation for the LUT by extracting required data from the given transitions and landmate-pft-maps files.
   - `prepare_mcgrath`: True/False. Preparation of Mcgrath data (Optional) for the LUT by extracting required data from given mcgrath file.
   - `remap`: remapping method used for extracting data in the preparation data section (**bilinear** or **con2**).
   - `scenario`: choose scenario name ("historical", "historical_high", "historical_low", "rcp19", "rcp26", "rcp34", "rcp45", "rcp60", "rcp70", "rcp85").
   - `grid`: choose the resolution in degrees.
   - `coords`(Optional): add custom coordinates for the selected region.
   ### Land-seas Mask configuration
   By default, the land-sea mask will be calculated from the Landmate PFT maps based on land classification. However, if you want to use a custom land-sea mask, you can specify its file path:
   - (`path_file_lsm`): Specify the path to the land-sea mask file if you prefer to use a custom one.
   - (`rcm_lsm_var`): variable name in the RCM land-sea mask file in case of having specified another file in `path_file_lsm`.
   ### Background Data configuration
   The project already includes global background data by default, which will be used unless you specify an alternative:
   - (`path_file_back*`): path to a custom background dataset.
   ### Optional file paths
   - (`path_file_*`): Specify file path for the used files. In case of not specifying default locations will be checked by the programm. 


 

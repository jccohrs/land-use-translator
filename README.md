 
# High-Resolution Land Use and Land Cover (LULC) Dataset Calculator for Regional Climate Modeling 


This Python program generates high-resolution land use and land cover (LULC) datasets for regional climate modeling across historical and future periods. By integrating this LULC data, climate models can more accurately simulate the impacts of land use changes on regional climate dynamics.

## Dataset Requirements

The following datasets are required to run the program:

- **Land Use (LU) Transitions**: LU changes for the selected time period, with an option to filter larger datasets to match the specified configuration. These files can be downloaded from the [LUH Data Portal](https://luh.umd.edu/data.shtml) for various scenarios. for various scenarios.
   - transitions.nc: Main land use transitions data.
   - states.nc: States of land use, for each grid cell and year.
   - management.nc: Contains forest management data; this   file is only required if McGrath Forest Types is enabled (i.e., if mcgrath is set to True).
   - add_tree_cover.nc: Tree cover data, used if AddTree Data is active (i.e., if addtree is set to True).
- **Landmate PFTs Maps**: Contains data for 16 Plant Functional Types (PFTs), providing detailed vegetation characterization. The data can be downloaded from [WDC Climate](https://www.wdc-climate.de/ui/entry?acronym=LM_PFT_EUR_v1.1).
- **Sea-Land Mask (ESA-CCI)**: A sea-land mask file specific to the intended region.
- **Irrigation Data** *(Optional)*: Required if irrigation data is included in the calculations.
- **Background Data** *(Optional)*: Additional background data for enhanced simulation detail.

By default, all required and optional datasets should be located in a designated data/ directory. If different storage locations are prefered, alternative paths for each dataset can be specified in a configuration file. This will be further detailed in the Usage section bellow.

## Usage Instructions

1. **Configuration File**: The main configuration file is located at `config/main.yaml`. Modify the options in this file to customize the program:

   ### LUT configuration
   - (`region`): Choose from pre-configured regions ("Germany", "Europe", "WestAfrica"), or add a new region by providing the necessary grid files.
   - (`forward`): **True** for future simulation or **False** for historical simulation.
   - (`backgrd`): True/False. Optionally include background data if available.
   - (`mcgrath`): True/False for using mcgrath data in the LUT. 
   - (`addtree`): True/False for using addtree data in the LUT. 
   - (`irri`): True/False. Enable or disable irrigation data, if the irrigation dataset is available.
   - (`syear`)/(`eyear`): Specify the time period for LU calculations by setting the starting year (`syear`) and ending year (`eyear`).
   - (`mcgrath_eyear`): end year of mcgrath file (in case that its different from eyear).
   - (`npfts`): number of npfts used in the LUT.
   - (`xsize`): xsize of the region.
   - (`ysize`): ysize of the region. 
   - (`rcm_lsm_var`): variable name in the RCM LSM file.


   ### LUH2 prepare data configuration
   - (`prepare_luh2_data`): True/False. Files preparation for the LUT by extracting required data from the given transitions and landmate-pft-maps files (required if first time using this program).
   - (`prepare_mcgrath`): True/False. Preparation of Mcgrath data (Optional) for the LUT by extracting required data from given mcgrath file.
   - (`remap`): remapping method used for extracting data in the preparation data section (**bilinear** or **con2**).
   - (`scenario`): choose scenario name ("historical", "historical_high", "historical_low", "rcp19", "rcp26", "rcp34", "rcp45", "rcp60", "rcp70", "rcp85").
   - (`grid`): choose the resolution in degrees.

   ### Optional file paths
   - (`path_file_*`): Specify file path for the used files. In case of not specifying default locations will be checked by the programm. 

2. **Execution**: To run it just go to the project directory and


   ```bash
   $ python3 main.py

## Installation Guide

1. **Install pipenv**:

   ```bash
   $ pip install pipenv

2. **Change to the project directory**:

   ```bash
   $ cd project_name

3. **Install the required libraries**:

   ```bash
   $ pipenv install 

4. **Activate your pipenv environment**:
   ```bash
   $ pipenv shell 

## References

 * https://essd.copernicus.org/articles/14/1735/2022/
 * https://essd.copernicus.org/articles/15/3819/2023/

 

 
# High-Resolution Land Use and Land Cover (LULC) Dataset Calculator for Regional Climate Modeling 


This Python program generates high-resolution land use and land cover (LULC) datasets for regional climate modeling across historical and future periods. By integrating this LULC data, climate models can more accurately simulate the impacts of land use changes on regional climate dynamics.

## Dataset Requirements

The following datasets are required to run the program:

- **Land Use (LU) Transitions**: LU changes for the selected time period, with an option to filter larger datasets to match the specified configuration.
- **Landmate PFTs Maps**: Contains data for 16 Plant Functional Types (PFTs), providing detailed vegetation characterization.
- **Sea-Land Mask (ESA-CCI)**: A sea-land mask file specific to the intended region.
- **McGrath Forest Types** *(Optional)*: Optional forest cover data for the specified period.
- **Irrigation Data** *(Optional)*: Required if irrigation data is included in the calculations.
- **Background Data** *(Optional)*: Additional background data for enhanced simulation detail.

## Regional and Temporal Scope

The program is pre-configured to calculate LULC for the following regions:

- **Africa**
- **West Africa**
- **Europe**
- **Germany**

To add other regions, include the appropriate grid files in the configuration.

## Usage Instructions

1. **Configuration File**: The main configuration file is located at `config/main.yaml`. Modify the options in this file to customize the program:

   - **Select Region**: Choose from pre-configured regions, or add a new region by providing the necessary grid files.
   - **Define Time Period**: Specify the time period for LU calculations by setting the starting year (`syear`) and ending year (`eyear`).
   - **Optional Dataset Inclusion**:
     - **McGrath Forest Types**: Choose to include or exclude this dataset.
     - **Irrigation Data**: Enable or disable irrigation data, if the irrigation dataset is available.
     - **Background Data**: Optionally include background data if available.

2. **Execution**: Once configured, the program will process the data, filtering and aligning datasets according to the selected options. To run it just go to the project directory and


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

 

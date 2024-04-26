import yaml
from src.lut import LUT
from src.utils import dotdict
from config.validation import validate_config
import sys


def main():
    with open("config/main.yaml") as stream:
        try:
            config = yaml.safe_load(stream)
            validate_config(config)
            return config
        except yaml.YAMLError as exc:
            print(exc)

if __name__ == '__main__':
    config = dotdict(main())

    LUT = LUT(config)
    print("Generating Namelist")
    print("...")
    namelist = LUT.generate_namelist()
    print("Succesfully generated")
    print("______________________")
    print("Preparing LUH2 data")
    print("...")
    LUT.prepare_luh2_data()
    print("Succesfully generated")
    print("______________________")
    print("Preparing MCGRATH data")
    print("...")
    LUT.prepare_mcgrath()
    print("Succesfully generated")
    print("______________________")
    ## perform changes in landcover changes
    if config.forward:
        print('CALCULATE LAND USE CHANGES FORWARD IN TIME')
        print("...")
        LUT.lucas_lut_forward()
    else:
        print('CALCULATE LAND USE CHANGES BACKWARD IN TIME')
        print("...")
        LUT.lucas_lut_backward()
    print('FINSHED LAND USE CHANGES')
    print("______________________")
    # split crops into irrigated and non-irrigated crops
    if config.irri:
        print('SPLIT CROPS USING IRRIGATION FRACTION')
        print("...")
        LUT.lucas_lut_irrigation()

    # split crops into irrigated and non-irrigated crops
    if config.mcgrath:
        print('ADJUST FOREST RELATIVE FOREST FRACTIONS USING MCGRATH DATA')
        print("...")
        LUT.lucas_lut_mcgrath()
    
    # Write out the data
    print('WRITE OUT DATA')
    LUT.lucas_lut_output()
    #print(f"RESULTS WRITTEN TO {namelist["F_LC_OUT"]}")
    print('LUCAS LUT SUCCESSFULLY FINISHED')
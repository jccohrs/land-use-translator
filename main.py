import yaml
from src.lut import LUT
from src.utils import dotdict
from config.validation import validate_config, validate_namelist_files

def load_configuration():
    with open("config/main.yaml") as stream:
        try:
            config = yaml.safe_load(stream)
            return dotdict(config)
        except yaml.YAMLError as exc:
            print(exc)

def prepare_lut(config):
    validate_config(config)
    lut = LUT(config)
    namelist = lut.generate_namelist()
    validate_namelist_files(namelist)
    return lut

def print_section_heading(heading):
    print("______________________")
    print(heading)
    print("...")

def main():
    config = load_configuration()
    lut = prepare_lut(config)
    
    print_section_heading("Preparing LUH2 data")
    lut.prepare_luh2_data()
    print("Successfully generated")

    print_section_heading("Preparing MCGRATH data")
    lut.prepare_mcgrath()
    print("Successfully generated")

    print_section_heading("Calculating land use changes")
    if config.forward:
        print('CALCULATE LAND USE CHANGES FORWARD IN TIME')
        lut.lucas_lut_forward()
    else:
        print('CALCULATE LAND USE CHANGES BACKWARD IN TIME')
        lut.lucas_lut_backward()
    print('FINISHED LAND USE CHANGES')

    if config.irri:
        print_section_heading("Splitting crops using irrigation fraction")
        lut.lucas_lut_irrigation()

    if config.mcgrath:
        print_section_heading("Adjusting forest relative forest fractions using McGrath data")
        lut.lucas_lut_mcgrath()

    print_section_heading("Writing out data")
    lut.lucas_lut_output()
    print('LUCAS LUT SUCCESSFULLY FINISHED')

if __name__ == '__main__':
    main()

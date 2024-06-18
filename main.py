import yaml
from src.lut import LUT
from src.utils import dotdict, print_section_heading
from config.validation import validate_config, validate_prepared_files

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
    return lut

def main():
    config = load_configuration()
    lut = prepare_lut(config)
    
    if config.prepare_luh2_data:
        print_section_heading("Preparing LUH2 data")
        lut.prepare_luh2_data()
        print("Successfully generated")
    if config.prepare_mcgrath:
        print_section_heading("Preparing MCGRATH data")
        lut.prepare_mcgrath()
        print("Successfully generated")
    namelist = lut.generate_namelist()
    validate_prepared_files(namelist)
    print_section_heading("Calculating land use changes")
    if config.forward:
        print('CALCULATE LAND USE CHANGES FORWARD IN TIME')
        lut.lucas_lut_forward()
    else:
        print('CALCULATE LAND USE CHANGES BACKWARD IN TIME')
        lut.lucas_lut_backward()
    print('FINISHED LAND USE CHANGES')

    print_section_heading("Writing out data")
    lut.lucas_lut_output()
    print('LUCAS LUT SUCCESSFULLY FINISHED')

if __name__ == '__main__':
    main()

import yaml
from src.lut import LUT
from src.utils import dotdict, print_section_heading
from config.validation import validate_config, validate_prepared_files, validate_main_files

def load_configuration():
    with open("config/main.yaml") as stream:
        try:
            config = yaml.safe_load(stream)
            return dotdict(config)
        except yaml.YAMLError as exc:
            print(exc)

def main():
    # Load configuration
    config = load_configuration()

    # Validate configuration
    validate_config(config)

    # Initialize LUT class
    lut = LUT(config)

    # Generate namelist
    namelist = lut.generate_namelist()

    # Validate main files
    validate_main_files(namelist, config)

    # Preparing the data for the lut calculation
    if config.prepare_luh2_data:
        print_section_heading("Preparing LUH2 data")
        lut.func_prepare_luh2_data()
    if config.prepare_mcgrath and not config.forward:
        print_section_heading("Preparing MCGRATH data")
        lut.func_prepare_mcgrath()
    
    # validating the prepared files
    validate_prepared_files(namelist, config)

    # Running the LUT calculation
    print_section_heading("Calculating land use changes")
    if config.forward:
        print('FORWARD IN TIME')
        lut.lucas_lut_forward()
    else:
        print('BACKWARD IN TIME')
        lut.lucas_lut_backward()
    
    # Writing out the data
    print_section_heading("Writing out data")
    lut.lucas_lut_output()
    print('LUCAS LUT SUCCESSFULLY FINISHED')

if __name__ == '__main__':
    main()

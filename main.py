import yaml
from src.processers import ProcesserClass
from src.utils import dotdict
from config.validation import validate_config
import sys

sys.tracebacklimit = 0

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

    PRC = ProcesserClass(config)
    print("Generating Namelist dicttionary")
    print("...")
    namelist = PRC.generate_namelist()
    print("Succesfully generated")
    print("Preparing LUH2 data")
    print("...")
    PRC.prepare_luh2_data()
    print("Succesfully generated")
    print("Preparing MCGRATH data")
    print("...")
    PRC.prepare_mcgrath()
    print("Succesfully generated")
    ## Next Steps 
    # LUT
    # OUTPUT
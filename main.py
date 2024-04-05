import yaml
from src.processers import ProcesserClass
from src.utils import dotdict

def main():
    with open("config/main.yaml") as stream:
        try:
            config = yaml.safe_load(stream)
            return config
        except yaml.YAMLError as exc:
            print(exc)

if __name__ == '__main__':
    config = dotdict(main())

    # Generate Namelist
    PRC = ProcesserClass(config)
    namelist = PRC.generate_namelist()
    PRC.prepare_luh2_data()
    print(namelist)
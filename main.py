import yaml
from src.processers import generate_namelist
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
    generate_namelist(config)

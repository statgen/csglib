import yaml

def read_yaml(filepath):
  with open(filepath,"rt") as fp:
    data = yaml.safe_load(fp)

  return data


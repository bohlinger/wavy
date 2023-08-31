from dataclasses import dataclass
from wavy.wconfig import load_or_default
from typing import Any, Dict, List

# initialize defaults from config files

@dataclass
class config_class:
    nID: str = None
    name: dict = None
    download: dict = None
    wavy_input: dict = None
    wavy_output: dict = None
    reader: str = None
    collector: str = None
    vardef: dict = None
    coords: dict = None
    misc: dict = None
    tags: list = None

def parse_config_file(obs_type: str, nID: str) -> dict:
    config_file_str = obs_type + '_cfg.yaml'
    parsed_file = load_or_default(config_file_str)
    return parsed_file

def dict_to_class(parsed_file: dict[Any, Any]) -> config_class:
    return config_class(**parsed_file)

def init_class(obs_type=None, nID=None) -> config_class:
    parsed_file = parse_config_file(obs_type, nID)
    cfg = dict_to_class(parsed_file[nID])
    cfg.misc['obs_type'] = obs_type
    dc = dict_to_class(parsed_file[nID])
    dc.nID = nID
    return dc

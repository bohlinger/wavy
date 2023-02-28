from dataclasses import dataclass
from wavy.wconfig import load_or_default
from typing import Any, Dict, List

# initialize defaults from config files

@dataclass
class config:
    nID: str = None
    names: dict = None
    retrieval: dict = None
    output: dict = None
    reader: str = None
    variables: dict = None
    misc: dict = None

    #def __post_init(self):
    #    self.


def parse_config_file(obs_platform_type: str, nID: str) -> dict:
    config_file_str = obs_platform_type + '_cfg.yaml'
    parsed_file = load_or_default(config_file_str)
    return parsed_file


def dict_to_class(parsed_file: Dict[Any, Any]) -> config:
    return config(**parsed_file)


def init_class_defaults(obs_platform_type,nID) -> config:
    parsed_file = parse_config_file(obs_platform_type,nID)
    cfg = dict_to_class(parsed_file[nID])
    cfg.misc['obs_platform_type'] = obs_platform_type
    return dict_to_class(parsed_file[nID])

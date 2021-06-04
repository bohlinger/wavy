import os
import yaml
import logging
logger = logging.getLogger(__name__)

# get envs for wavy
from dotenv import load_dotenv
load_dotenv()
confdir = os.getenv('WAVY_CONFIG')

def load_or_default(name):
    """
    Load the YAML configuration file from the current directory or use default.

    e.g.:

    .. code::

        c = load_or_default('/model_specs.yaml')

    """
    logging.debug('attempting to load: %s..' % name)
    try:
        #with open(name, 'r') as s:
        if confdir is not None:
            filestr = confdir + name
        else: filestr = name
        with open(filestr, 'r') as s:
            return yaml.safe_load(s)

    except FileNotFoundError:
        logging.debug('could not load from local directory, using default.')
        from pkg_resources import resource_stream
        #return yaml.safe_load(resource_stream(__name__, name + '.default'))
        return yaml.safe_load(resource_stream(__name__, 'config/' \
                            +  name + '.default'))

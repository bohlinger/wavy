import os
import yaml
import logging
logger = logging.getLogger(__name__)


def load_or_default(name):
    """
    Load the YAML configuration file from the current directory or use default. The approach is threefold and shown in order:
    1. check if env 'WAVY_CONFIG' is set
    2. check if .env exists
    3. check if a config folder exists using xdg
    4. fall back to default files within the package

    e.g.:

    .. code::

        c = load_or_default('/model_specs.yaml')

    """
    logging.debug('attempting to load: %s..' % name)
    # get envs for wavy
    # 1. check if already in envs
    confdir = os.getenv('WAVY_CONFIG')
    # 2. if not check if envs are in .env
    if confdir is None:
        from dotenv import load_dotenv
        load_dotenv()
        confdir = os.getenv('WAVY_CONFIG')
    # 3. user wide, xdg
    import xdg
    confdir = os.path.join(xdg.xdg_config_home(),'wavy')

    try:
        #with open(name, 'r') as s:
        if confdir is not None:
            filestr = os.path.join(confdir,name)
        else: filestr = name
        with open(filestr, 'r') as s:
            return yaml.safe_load(s)

    except FileNotFoundError:
        logging.debug('could not load from local directory, using default.')
        from pkg_resources import resource_stream
        #return yaml.safe_load(resource_stream(__name__, name + '.default'))
        return yaml.safe_load(resource_stream(__name__, 'config/' \
                            +  name + '.default'))

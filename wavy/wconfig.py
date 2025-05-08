import os
import yaml
import logging

logger = logging.getLogger(__name__)
import dotenv
import xdg



def __get_confdir__():
    """
    Configuration directory is specified through:

    .env
    WAVY_CONFIG environment variable
    XDG configuration directory (wavy)
    """
    dotenv.load_dotenv(dotenv_path=dotenv.find_dotenv(usecwd=True))
    c = os.getenv('WAVY_CONFIG', None)
    if c is None:
        c = os.path.join(xdg.xdg_config_home(), 'wavy')

        if not os.path.exists(c):
            c = None

    logger.debug('config directory: %s' % c)
    return c


def load_or_default(name):
    """
    Load the YAML configuration file from the current directory or use default. The approach is threefold and shown in order:

    1. check if env 'WAVY_CONFIG' is set or specified in .env
    2. check if a config folder exists using xdg
    3. fall back on default files within the package

    e.g.:

    .. code::

        c = load_or_default('model_cfg.yaml')

    """
    logging.debug('attempting to load: %s..' % name)

    confdir = __get_confdir__()

    try:
        if confdir is not None:
            filestr = os.path.join(confdir, name)
        else:
            raise FileNotFoundError()
        with open(filestr, 'r') as s:
            return yaml.safe_load(s)

    except FileNotFoundError:
        logging.debug('could not load from local directory, using default.')
        from pkg_resources import resource_stream
        return yaml.safe_load(
            resource_stream(__name__,
                            os.path.join('config', name + '.default')))

def load_minimal(name):
    logging.debug('attempting to load: %s..' % name)

    from pkg_resources import resource_stream
    return yaml.safe_load(resource_stream(__name__,
                          os.path.join('config', name + '.minimal')))

def load_dir(name):
    from pkg_resources import resource_stream
    #return resource_stream('wavy', name + '.py')
    return resource_stream(__name__, name + '.py')

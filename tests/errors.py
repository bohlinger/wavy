"""
Central exception hierarchy for wavy.

Guidelines:
- Every wavy-raised exception inherits from WavyError, so callers can
  do `except WavyError` to catch "anything wavy went wrong" without
  swallowing unrelated Python exceptions (KeyError, ValueError, etc.
  from bugs or bad input should NOT be silently caught this way).
- Group by domain (config, model, satellite, station, io) so callers
  can be as broad or as specific as they need.
- Keep messages actionable: say what was wrong AND what to check.
"""


class WavyError(Exception):
    """Base class for all wavy-specific exceptions."""

    pass


# --- configuration -------------------------------------------------- #


class ConfigError(WavyError):
    """Something is wrong with a config file (model_cfg.yaml, etc.)."""

    pass


class MissingConfigKeyError(ConfigError):
    """A required key is missing from a config entry."""

    pass


# --- model data ------------------------------------------------------- #


class ModelError(WavyError):
    """Base class for model_class-related errors."""

    pass


class ModelFileSearchError(ModelError):
    """
    Raised when the 'best guess' search for an accessible model file
    exceeds the configured number of attempts (max_iter). Almost
    always means src_tmplt/fl_tmplt in the config are wrong.
    """

    pass


class ModelFileNotFoundError(ModelError):
    """No model file could be found for the requested period at all."""

    pass


# --- satellite data --------------------------------------------------- #


class SatelliteError(WavyError):
    """Base class for satellite_module-related errors."""

    pass


# --- station data ------------------------------------------------------ #


class StationError(WavyError):
    """Base class for station_module-related errors."""

    pass

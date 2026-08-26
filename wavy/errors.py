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


<<<<<<< HEAD
class ModelPathTemplateError(ConfigError):
    """
    Raised when building a filesystem path from a model's
    src_tmplt/strsub configuration fails. Indicates a config problem,
    not a per-date data-availability problem.
    """

    pass


class RegionNotDefinedError(ConfigError):
    """
    Raised when a requested region is not defined in any known region
    source (region_cfg.yaml's rect/poly/geojson sections, or as a
    model grid in model_cfg.yaml).
    """

    pass


=======
>>>>>>> 13f3198 (fixup! feat(error module): error classes wavy specific - each class targets one aspect of wavy to filter specifically the errors from one module.)
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


class ModelReadError(ModelError):
    """Raised when a model's reader function fails to read data."""

    pass


class ModelProcessingError(ModelError):
    """
    Raised when post-read processing (renaming, CF compliance,
    convention enforcement, longitude formatting) fails after the
    reader itself succeeded.
    """

    pass


class ModelFileNotFoundError(ModelError):
    """No accessible model files were found for the requested period."""

    pass


class GridRetrievalError(ModelError):
    """
    Raised when model grid coordinates needed for region matching
    could not be retrieved, even after falling back to the model's
    default/configured grid_date.
    """

    pass


# --- satellite data --------------------------------------------------- #


class SatelliteError(WavyError):
    """Base class for satellite_module-related errors."""

    pass


<<<<<<< HEAD
class SatellitePathTemplateError(ConfigError):
    """
    Raised when building a filesystem path from a satellite's
    src_tmplt/strsub configuration fails. Indicates a config problem,
    not a per-date data-availability problem.
    """

    pass


class SatelliteFileNotFoundError(SatelliteError):
    """No accessible satellite files were found for the requested period."""

    pass


class SatelliteReadError(SatelliteError):
    """
    Raised when none of the candidate satellite files could be read
    successfully.
    """

    pass


class SatelliteProcessingError(SatelliteError):
    """
    Raised when post-read processing (renaming, CF compliance,
    convention enforcement, longitude formatting) fails after at
    least one file was read successfully.
    """

    pass


class SatelliteVariableError(SatelliteError):
    """
    Raised when a required coordinate/variable name could not be
    resolved or renamed in the satellite data.
    """

    pass


# --- station data ------------------------------------------------------ #


class StationError(WavyError):
    """Base class for station_module-related errors."""

    pass


# --- collocation ---------------------------------------------------------- #


class CollocationError(WavyError):
    """Base class for collocation_class-related errors."""

    pass


class CollocationInputError(CollocationError):
    """
    Raised when collocation preconditions are not met - e.g. no
    observation data (oco), no model specified, or an unrecognized
    collocation method.
    """

    pass


class CollocationBuildError(CollocationError):
    """
    Raised when the collocated xarray Dataset could not be assembled
    from the intermediate collocation results (missing/mismatched
    keys between results_dict and variable_def.yaml).
    """

    pass


class CollocationRunError(CollocationError):
    """
    Raised when the collocation process itself fails - either an
    unexpected error during populate()/collocate(), or every
    candidate forecast date failed to collocate.
    """

    pass

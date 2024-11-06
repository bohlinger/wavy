Setting up a project and **wavy** config files
##############################################
**wavy** can be run with your personalized config-files. For basically any data wavy is using or writing you need the apropriate config file customized to your needs. This is straight forward, easy to achieve and will be demonstrated in this section.

Default **wavy** config files are delivered with the code. These can be copied to your project folder and adjusted to your needs. The following default config files exist in wavy/wavy/config:


.. code-block:: bash

   $ ls
   model_cfg.yaml.default          region_cfg.yaml.default
   d22_var_dicts.yaml.default      satellite_cfg.yaml.default
   insitu_cfg.yaml.default         validation_metrics.yaml.default
   quicklook_cfg.yaml.default      variable_def.yaml.default
   variable_frost.yaml.default



The naming of the files is descriptive but here some brief description:
        * insitu_cfg* -> for insitu module, when using insitu data
        * satellite_cfg* -> for satellite module, when using satellite data
        * model_cfg* -> for model module, when using model output data
        * d22_var_dicts* -> extra config for specifying content of .d22 files
        * region_cfg* -> for specifying your regions of interest
        * validation_metrics* -> for validation paths
        * variable_def* -> specifying standard names, variable abbreviations, variable attributes
        * quicklook_cfg* -> for customizing quicklook figures
        * variable_frost* -> for variables as defined in FROST API

In order to create your own customized project, copy the needed config default files to your directory of choice (aka project directory), remove the ".default" extension, and ammend them. E.g. when only using satellite products within your regions of interest you only need (if you do not want to change default variable settings):

        * satellite_cfg* -> for satellite module, when using satellite data
        * region_cfg* -> for specifying your regions of interest

**wavy** browses the directory structure as follows:

    * check if env 'WAVY_CONFIG' is set or specified in .env
    * fall back on default files within the package


This means, to use your own version of the config files you can create an .env file in your **wavy** directory whith an environmental variable pointing to where your config files are located:

.. code-block:: bash

        WAVY_CONFIG=path/to/your/config/files

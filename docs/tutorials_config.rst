Setting up a project and **wavy** config files
##############################################
**wavy** can be run for various projects with different settings. This is organized by your personalized config-files. For basically any data wavy is using or writing you need the apropriate config file customized to your needs. This is straight forward, easy to achieve and will be demonstrated in this section.

Default **wavy** config files are delivered with the code. These can be copied to your project folder and adjusted to your needs. The following default config files exist in wavy/wavy/config:


.. code-block:: bash

   $ ls
   collocation_specs.yaml.default  region_specs.yaml.default
   d22_var_dicts.yaml.default      satellite_specs.yaml.default
   insitu_specs.yaml.default       validation_specs.yaml.default
   model_specs.yaml.default        variable_info.yaml.default
   quicklook_specs.yaml.default


The naming of the files is descriptive but here some brief description:
        * insitu_specs* -> for insitu module, when using insitu data
        * satellite_specs* -> for satellite module, when using satellite data
        * model_specs* -> for model module, when using model output data
        * d22_var_dicts* -> extra config for specifying content of .d22 files
        * region_specs* -> for specifying your regions of interest
        * collocation_specs* -> for the collocation paths
        * validation_specs* -> for validation paths
        * variable_info* -> specifying standard names, variable abbreviations, variable attributes
        * quicklook_specs* -> for customizing quicklook figures


In order to create your own customized project, copy the needed config default files to your directory of choice (aka project directory), remove the ".default" extension, and ammend them. E.g. when only using satellite products within your regions of interest you only need (if you do not want to change default variable settings):

        * satellite_specs* -> for satellite module, when using satellite data
        * region_specs* -> for specifying your regions of interest

**wavy** browses the directory structure as follows:

    * check if env 'WAVY_CONFIG' is set or specified in .env
    * check if a config folder exists using xdg
    * fall back on default files within the package


This means you can create an .env file in your project directory. This .env file needs an environmental variable pointing to where your project files are located:

.. code-block:: bash

        WAVY_CONFIG=path/to/your/config/files

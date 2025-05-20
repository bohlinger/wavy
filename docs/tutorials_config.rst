Setting up a project and **wavy** config files
##############################################

Setup
-----
For most use cases the following workflow may be advantageous:

#. create a project directory e.g.

   .. code-block:: bash

       $ mkdir ~/my_wavy_project

#. create a config directory therein

   .. code-block:: bash

       $ cd ~/my_wavy_project
       $ mkdir config

#. create an .env file that will contain the path to your config directory
   
   .. code-block:: bash

       $ cd ~/my_wavy_project
       $ touch .env

In the .env-file you write for this case WAVY_CONFIG=/home/{USER}/my_wavy_project/config/, where you substitute {USER} with your username. Once these steps are concluded you are good to go!

On the config files
-------------------
**wavy** can be run with your personalized config-files. For basically any data **wavy** is using you can customize the according config file to your needs. This is straight forward, easy to achieve and will be demonstrated in this section.

Default **wavy** config files are delivered with the code. These can be copied to your project folder and adjusted to your needs. The following default config files are obtainable:


.. code-block:: bash

   $ ls
   model_cfg.yaml          region_cfg.yaml
   d22_var_dicts.yaml      satellite_cfg.yaml
   insitu_cfg.yaml         validation_metrics.yaml
   quicklook_cfg.yaml      variable_def.yaml
   variable_frost.yaml

The naming of the files is descriptive but here some brief description:
        * insitu_cfg* -> for insitu module, when using insitu data
        * satellite_cfg* -> for satellite module, when using satellite data
        * model_cfg* -> for model module, when using model output data
        * d22_var_dicts* -> extra config for specifying content of .d22 files
        * region_cfg* -> for specifying your regions of interest
        * validation_metrics* -> for defining names of validation metrics
        * variable_def* -> specifying standard names, variable abbreviations, variable attributes
        * quicklook_cfg* -> for customizing quicklook figures
        * variable_frost* -> for variables as defined in FROST API

In all config files there are some default settings which you usually have to customize. You can obtain the default config files by using **wavy**'s function:

.. code-block:: bash

   $ mamba activate wavyopen
   $ wavyCFG --help
   $ wavyCFG --path ~/my_wavy_project/config/. --f satellite_cfg.yaml

"wavyCFG --help" will give you instruction on how to proceed but in general the above line is how wavyCFG can be executed. The satellite, insitu, and model config files have minimal version of the config files that are easier to ammend and this can be evoked by adding the --t flag for overwriting the default:

.. code-block:: bash

   $ wavyCFG --path ~/my_wavy_project/config/. --f satellite_cfg.yaml --t minimal

E.g. when only using satellite products within your regions of interest you only need (if you do not want to change default variable settings):

        * satellite_cfg* -> for satellite module, when using satellite data
        * region_cfg* -> for specifying your regions of interest

**wavy** browses the directory structure as follows:

    * check if env 'WAVY_CONFIG' is set or specified in .env
    * if nothing is found, fall back on default files within the package

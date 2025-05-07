Tutorials
==========
The following tutorials describe how to setup the wavy config files, to retrieve satellite and insitu data, to perform collocation and validation tasks. **wavy** is intended to be fleksibel such that customization can be achieved with minimal changes in the code. The **wavy** config files build the fundament for this approach and serve a similar purpose as namelist files often used for runtime changes in numerical modelling.

In general, executable files usually have help function which can be read using e.g.:

.. code-block:: bash

   $ ./{Executable}.py --help

e.g.:

.. code-block:: bash

   $ cd ~/wavy/wavy/apps
   $ ./wavyDownload.py --help

.. toctree::
   :maxdepth: 2

   tutorials_config
   tutorials_sat
   tutorials_modelmod
   tutorials_insitu
   tutorials_stormtrack
   tutorials_consolidate
   tutorials_collocmod
   tutorials_validate
   tutorials_stats
   tutorials_triple_collocation
   tutorials_gridder
   tutorials_ais

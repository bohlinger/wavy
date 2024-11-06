Installation and setup
======================

Installing wavy can be done via conda. The steps are as follows:

#. clone the github repo like:

   .. code-block:: bash

      $ cd ~
      $ git clone https://github.com/bohlinger/wavy.git

      or for one single branch try:
      $ git clone --single-branch --branch master https://github.com/bohlinger/wavy.git


#. install wavy:

   .. code-block:: bash

      $ cd ~/wavy
      $ conda env create -f environment.yml
      $ conda activate wavy

A much faster installation method would be using mamba if you have that installed.

   .. code-block:: bash

      $ cd ~/wavy
      $ mamba env create -f environment.yml
      $ conda activate wavy


Now, append wavy root directory to $PYTHONPATH, for instance add the following to your .bashrc:

   .. code-block:: bash

      export PYTHONPATH=$PYTHONPATH:/path/to/your/wavy
      
.. note::

   /path/to/your/wavy/ should be replace with the full path of your wavy folder. It will be the case throughout all this documentation.

Create an .env file in your wavy directory and point to the location where you plan to store your custom config files. Your .env-file could look like:

   .. code-block:: bash

      WAVY_CONFIG=/path/to/your/config/

How to start your own project and how to manage the config files is explained in a separate tutorial.

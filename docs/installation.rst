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

      export PYTHONPATH=$PYTHONPATH:/home/patrikb/wavy

.. note::

   My user *patrikb* was used in this example. This will be the case in more examples and needs to be adjusted for other users.

Create an .env file in your wavy directory and point to the location where you plan to store your custom config files, as well as your wavy directory. Your .env-file could look like:

   .. code-block:: bash

      WAVY_CONFIG=path/to/your/config/files
      WAVY_DIR=/path/to/your/wavy/

How to start your own project and how to manage the config files is explained in a separate tutorial.

Installation and setup
======================

Installation for regular use
---------------------------
For regular use **wavy** wavy can simply be installed via pip, conda, or mamba by installing the wavyopen package. So please choose between the following options:


.. code-block:: bash

   $ pip install wavyopen


.. code-block:: bash

   $ mamba create --name wavyopen --channel=conda-forge
   $ mamba activate wavyopen
   $ mamba install wavyopen


.. code-block:: bash

   $ conda create --name wavyopen --channel=conda-forge
   $ conda activate wavyopen
   $ conda install wavyopen -c conda-forge

When activated, you can use wavyopen in any directory of your computer. In case you want to renew or update your environment, it is often best to remove it completely and install again as described above. The removal can be done with:

.. code-block:: bash

   $ conda remove -n wavyopen --all

.. note::

   All code examples and tutorials will assume that you have installed wavy using conda or mamba. In case you still have, from older installations or the development installation described below, environmental variables in e.g. your .bashrc that attribute paths for wavy, please remove them before you follow the installation procedure above.

Installation for development
----------------------------
Installing **wavy** can be done via conda. The steps are as follows:

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

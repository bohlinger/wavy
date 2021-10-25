Installing
==========
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

Now, append wavy root directory to $PYTHONPATH, for instance add the following to your .bashrc:

   .. code-block:: bash

      export PYTHONPATH=$PYTHONPATH:/home/patrikb/wavy

Either create a xdg project path or an .env file in the wavy root directory and point to the location wher you plan to store your custom config files. This could look like:

   .. code-block:: bash

      WAVY_DIR=/home/patrikb/wavy/
      WAVY_CONFIG=/home/patrikb/wavy/config/

#. adjust config files --> done in tutorials

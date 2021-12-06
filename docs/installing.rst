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

.. note::

   My user *patrikb* was used in this example. This will be the case in more examples and needs to be adjusted for other users.

Either create a xdg project path or an .env file in the wavy root directory and point to the location wher you plan to store your custom config files. This could look like:

   .. code-block:: bash

      WAVY_DIR=/home/patrikb/wavy/
      WAVY_CONFIG=/home/patrikb/wavy/config/


#. prepare config files

   Copy config files to WAVY_CONFIG and remove the .default suffix. Now edit the config files you will be using --> edits done in tutorial

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
      $ pip install -e .

#. adjust config files --> done in tutorials

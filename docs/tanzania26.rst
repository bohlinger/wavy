Tanzania Validation Workshop 2026
================================

The following examples are tailored to the **wavy** Tanzania Validation Workshop 2026.

0. Installation of wavyopen and preparation
###########################################
**wavyopen** is wavy's pip and conda package which can easily be installed with pip, conda, or mamba. In this workshop we will install the conda package in a fresh conda environment like:

.. code::

   conda create --name wavyopen --channel=conda-forge
   conda activate wavyopen
   conda install wavyopen -c conda-forge

Assuming you have mamba installed you can, instead of conda, use mamba as an inplace replacement like:

.. code::

   mamba create --name wavyopen --channel=conda-forge
   mamba activate wavyopen
   mamba install wavyopen

When activated, you can use wavyopen in any directory of your computer. In case you want to renew or update your environment, it is often best to remove it completely and install again as described above. The removal can be done with:

.. code::

   conda remove -n wavyopen --all
   # or
   mamba remove -n wavyopen --all

Since you probably have multiple, independent projects at the same time it may make sense to follow a general workflow of creating a project directory, specifying your environment variables, preparing config files, creating scripts. Let's go through each of these steps together:

#. Create a project directory:

   .. code::
   
      mkdir ~/wavyopen_tanzania26
      mkdir -p ~/wavyopen_tanzania26/config
      mkdir -p ~/wavyopen_tanzania26/scripts
      mkdir -p ~/wavyopen_tanzania26/data


#. Specify enfironmental variables

   .. code::

      cd ~/wavyopen_tanzania26
      touch .env

   Specify now where you will have your **wavy** config files in the .env-file.

   .. code::

      WAVY_CONFIG=/home/{USER}/wavyopen_tanzania26/config/

   {USER} has to be substituted with your username. Also, if you would like to download data you need to have you copernicus marine credentials stored in the environment. This can be done anywhere (.bashrc, .profilerc, any shell wrapper script, manually, ...) and it could also be done in the .env file which we will do for this workshop.

   .. code::

      COPERNICUSMARINE_SERVICE_USERNAME=
      COPERNICUSMARINE_SERVICE_PASSWORD=

#. Preparing your config files:
   The config files can be established following some prepared **wavy** examples which you can obtain using the wavyCFG script like:

   .. code::

      cd ~/wavyopen_tanzania26/config
      wavyCFG --help
      wavyCFG --path ./. --f satellite_cfg.yaml --t minimal
      wavyCFG --path ./. --f model_cfg.yaml --t minimal
      wavyCFG --path ./. --f region_cfg.yaml --t default

#. Adjustments and further preparations are done in the subsequent sections.


1. Downloading data satellite data
##################################
First we need to adjust the satellite_cfg.yaml file for your purpose. Please open the satellite_cfg.yaml file and ammend it. Specifically, the target template for the downloads (trgt_tmplt) and the source template for files that wavy should use (src_tmplt) need to be defined. I am choosing here:

.. code::

   trgt_tmplt: /home/{USER}/wavyopen_tanzania26/data/name/%Y/%m/
   src_tmplt: /home/{USER}/wavyopen_tanzania26/data/name/%Y/%m/

Again, pleaser substitute {USER} with your username. 

In the same file, the time increment for copernicus can be changed to days, such that downloading the satellite files for several days is faster.

.. code::

   time_incr: d

Now you can use the wavyDownload script to download data from the copernicus marine service.

.. code::

   wavyDownload --help
   wavyDownload --sd 2026-08-17 --ed 2026-08-21 --nID cmems_L3_NRT --name s3a

You can repeat this for all the other satellites as well (s3a, c2, j3, h2b, al, cfo, s6a, swon). If you like to retrieve all satellite missions in the list **name** then you can replace the name of the satellite with **all** like:

.. code::

   wavyDownload --sd 2026-08-17 --ed 2026-08-21 --nID cmems_L3_NRT --name all


2. Process satellite data
#########################


Now, you can start preparing python scripts reading, processing, and plotting your data. This may look like:

.. code-block:: python3

   from wavy import sc, gc, ms
   from wavy.grid_stats import apply_metric

   # satellite data from directory
   sco = sc(nID='cmems_L3_NRT',
            name='s3a',
            sd='2026-08-17', ed='2026-08-21',
            region="NorthSea").populate()

   # plot results
   # map
   sco.quicklook(m=True)
   # time series
   sco.quicklook(ts=True)
   # all
   sco.quicklook(a=True)

   # multiple variables at once (significant wave height and 10m wind speed)
   sco = sc(nID='cmems_L3_NRT',
            name='s3a',
            sd='2026-08-17', ed='2026-08-21',
            region="NorthSea",
            varalias=['Hs','U']).populate()

   # plot results
   sco.quicklook(a=True, varalias='Hs')
   sco.quicklook(a=True, varalias='U')


   # satellite data from multiple sources
   mso = ms(nID=['cmems_L3_NRT'],
            name=['s3a', 's3b', 'c2', 'cfo', 'h2b', 'j3', 'al', 's6a', 'swon'],
            sd='2026-08-17', ed='2026-08-21',
            region='NorthSea')

   # plot results
   mso.quicklook(a=True, mode='indiv')

   # grid satellite data
   bb = (-5, 12, 50, 62)  # lonmin,lonmax,latmin,latmax
   res = (1, 1)  # lon/lat
   gco = gc(oco=mso, bb=bb, res=res)

   # compute metrics
   gridvar, lon_grid, lat_grid = apply_metric(gco=gco)

   # plot results
   gco.quicklook(val_grid=gridvar, lon_grid=lon_grid, lat_grid=lat_grid,
                 title="", metric='mor', land_mask_resolution='i')


Now, introduce your custom region in region_cfg.yaml and rerun the script by replacing "NorthSea" with what you defined. Here for example we define a bounding box corresponding to the domain of the ww3 model file. 

.. code-block:: yaml

   rect:
       global: {llcrnrlon: -180.0, llcrnrlat: -90.0, urcrnrlon: 180.0, urcrnrlat: 90.0}
       tanzania: {llcrnrlon: 30.0, llcrnrlat: -34, urcrnrlon: 80., urcrnrlat: 15.0}



3. Add custom model to wavy
###########################
Add the Tanzania relevant model output files to the model_cfg.yaml file. For instance you can add your ww3 model like:

.. code-block:: yaml

   ww3_tanzania:
       name: null
       download: null
       wavy_input: {src_tmplt: '/home/{USER}/wavyopen_tanzania26/data/model/ww3/',
                    fl_tmplt: ww3.nc}
       reader: read_era5 #read_ww3_4km
       collector: null
       vardef: {Hs: hs, time: time, lons: lon, lats: lat}
       coords: null
       misc:
           init_times: [0, 12]
           init_step: 12
           grid_date: null
           convention: meteorological
           date_incr_unit: h
           date_incr: 3

Where {USER} should be replaced with your own user name. 

Check if your model data is readable by wavy with:

.. code-block:: python3

    from wavy import mc

    mco = mc(nID='ww3_tanzania', sd='2026-08-20').populate()
    mco.quicklook(m=True)



4. Collocate satellite with model
#################################

Access to model and observations enables you to validate the model against the observations. This can be done using the collocation module like:

.. code-block:: python3

    from wavy import cc, ms

    mso = ms(nID=['cmems_L3_NRT'],
             name=['s3a', 's3b', 'c2', 'cfo', 'h2b', 'j3', 'al', 's6a', 'swon'],
             sd='2026-08-17', ed='2026-08-21',
             region='tanzania')

    mso = mso.filter_landMask().filter_distance_to_coast(llim=50000)

    cco = cc(model='ww3_tanzania', oco=mso, leadtime='best').populate()

    cco.quicklook(ts=True, m=True, sc=True, hist=True,
                  std_regression_line=True,
                  std_regression_col='b',
                  std_regression_lw=1)


5. Validate with model against satellite observations
#####################################################

Validation is quick and easy. Using the collocation class object **cco** you do:

.. code-block:: python3

   cco.validate_collocated_values()

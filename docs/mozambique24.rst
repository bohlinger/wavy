Mozambique Validation Workshop 2024
===================================


The following examples are tailored to the **wavy** Mozambique Validation Workshop in fall 2024. This workshop will focus on some simple examples that can be used as python code snippets in your workflow.

0. **wavy** installation
##############################

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


1. **wavy** config files
########################

Create a new project directory, in this example we will call it Moz_ws24_wavy. Within this folder, create another folder config, where you will store the configuration files used by wavy for your project. This could look like:

.. code-block:: bash

    :~$ mkdir Moz_ws24_wavy
    :~$ cd Moz_ws24_wavy
    :~/Moz_ws24_wavy$ mkdir config  
        
In the config folder located in the path/to/your/wavy/wavy/ directory, you will find some default config files. For this workshop you will need the following:

.. code-block:: bash

   satellite_cfg.yaml.default  region_cfg.yaml.default
   quicklook_cfg.yaml.default    validation_metrics.yaml.default
   model_cfg.yaml.default        variable_def.yaml.default     
        
Now copy these files listed above into the config folder you just created in your project directory (/Moz_ws24_wavy/config/) and remove the suffix *.default*. It should now look like:
 
 .. code-block:: bash

    :~/Moz_ws24_wavy/config$ ls
    satellite_cfg.yaml  region_cfg.yaml
    quicklook_cfg.yaml    validation_metrics.yaml
    model_cfg.yaml        variable_def.yaml  
 
At the root of your project directory, establish an .env file such that **wavy** knows where to find the config files it should use. This could look like:
 
.. code-block:: bash        
    :~/Moz_ws24_wavy$ touch .env

And the .env file should contain the following line:

.. code-block:: bash

     WAVY_CONFIG=/home/USER/Moz_ws24_wavy/config

Replace USER with you username that you get when typing

.. code-block:: bash

     echo ${USER}

In order to download satellite data, you also need to copy the wavyDownload.py file from /wavy/apps/standalone/ into your project directory. 
        
This should be the structure of your project directory:

.. code-block:: bash

        :~/Moz_ws24_wavy$ ls -la
        total 16
        drwxrwxr-x  3 user user 4096 Nov 14 09:04 .
        drwx------ 79 user user 4096 Nov 14 09:06 ..
        drwxrwxr-x  2 user user 4096 Nov 14 09:10 config
        -rwxr-xr-x  1 user user   44 Nov 14 09:04 .env
        -rwxrwxr-x  1 user user 3814 Nov  5 09:12 wavyDownload.py


2. Download L3 satellite altimetry data
#######################################

L3 satellite data is obtained from Copernicus with the product identifier WAVE_GLO_WAV_L3_SWH_NRT_OBSERVATIONS_014_001. User credentials are required for this task. So before you can start you have to get a Copernicus account (free of costs). Prepare access to Copernicus products. Your should add the following lines to your .bashrc, adapted with your username and password from Copernicus. 

.. code::

    export COPERNICUSMARINE_SERVICE_USERNAME=YOUR_COPERNICUS_USERNAME
    export COPERNICUSMARINE_SERVICE_PASSWORD=YOUR_COPERNICUS_PASSWORD


Adjust the satellite config file called *satellite_cfg.yaml*. Remember, this is the file you copied to *~/Moz_ws24_wavy/config*. In this file you should adapt the default paths with the ones from your project. It should include the following section and could look like:

.. code-block:: yaml

   --- # specifications for satellite missions

   cmems_L3_NRT:
        # mandatory
        name:
            s3a: s3a
            s3b: s3b
            c2: c2
            j3: j3
            h2b: h2b
            al: al
            cfo: cfo
            s6a: s6a
            swon: swon
        # mandatory when downloading
        # where to store downloaded data
        download:
            ftp: # downloading method
                src_tmplt: "/Core/\
                            WAVE_GLO_PHY_SWH_L3_NRT_014_001/\
                            cmems_obs-wave_glo_phy-swh_nrt_name-l3_PT1S/\
                            %Y/%m/"
                trgt_tmplt: /path/to/Moz_ws24_wavy/altimeter_data/L3/name/%Y/%m
                path_date_incr_unit: 'm'
                path_date_incr: 1
                search_str: '%Y%m%dT'
                strsub: ['name']
                server: "nrt.cmems-du.eu"
           copernicus:
                dataset_id: cmems_obs-wave_glo_phy-swh_nrt_name-l3_PT1S
                trgt_tmplt: /path/to/Moz_ws24_wavy/altimeter_data/L3/name/%Y/%m
                path_date_incr_unit: 'm'
                path_date_incr: 1
                strsub: ['name']
                server: "nrt.cmems-du.eu"
                time_incr: 'd' # 'h', 'd', 'm'
        # optional: where to read from
        #           can be defined directly when calling wavy
        wavy_input:
            src_tmplt: /path/to/Moz_ws24_wavy/altimeter_data/L3/name/%Y/%m
            fl_tmplt: "varalias_name_region_\
                            %Y%m%d%H%M%S_%Y%m%d%H%M%S.nc"
            strsub: ['name']
            path_date_incr_unit: 'm'
            path_date_incr: 1
        # optional: where to write to
        #           can be defined directly when calling wavy
        wavy_output:
            trgt_tmplt: /path/to/Moz_ws24_wavy/altimeter_data/L3/name/%Y/%m
            fl_tmplt: "varalias_name_region_\
                            %Y%m%d%H%M%S_%Y%m%d%H%M%S.nc"
            strsub: ['varalias','name','region']
            file_date_incr: m
        # optional, if not defined the class default is used
        reader: read_local_ncfiles
        collector: get_remote_files_copernicusmarine
        # optional, needs to be defined if not cf and in variable_info.yaml
        vardef:
            Hs: VAVH
            U: WIND_SPEED
        coords:
        # optional, info that can be used by class functions
        misc:
            processing_level:
            provider:
            obs_type:
        # optional, to ease grouping
        tags:


You can proceed now and download L3 data using the wavyDownload.py script you copied in your project folder. You can get help with:

.. code-block:: bash

   $ ./wavyDownload.py -h

And then download some satellite altimeter data:

.. code-block:: bash

   $ ./wavyDownload.py --name s3a --sd 20241017T07 --ed 20241017T08 --nID cmems_L3_NRT

If you need to download satellite data from Copernicus for more than a day or month, you can change the time increment in time_incr. 'h' will download 3-hours files at a time, 'd' will download all available files for a day at a time and 'm' all available files for a month at a time. Make sure to change this parameter if you need to download long periods of data as this will considerably shorten the time it takes to do so.  

.. code-block:: yaml
           copernicus:
                dataset_id: cmems_obs-wave_glo_phy-swh_nrt_name-l3_PT1S
                trgt_tmplt: /path/to/Moz_ws24_wavy/altimeter_data/L3/name/%Y/%m
                path_date_incr_unit: 'm'
                path_date_incr: 1
                strsub: ['name']
                server: "nrt.cmems-du.eu"
                time_incr: 'm' # 'h', 'd', 'm'

You can also download the data directly with python as follows:

.. code-block:: python3

    >>> from wavy.satellite_module import satellite_class as sc
    >>> nID = "cmems_L3_NRT"
    >>> name = "s3a"
    >>> sd = "2024-10-17 07"
    >>> ed = "2024-10-17 09"
    >>> sco = sc(sd=sd, nID=nID, name=name, ed=ed).download()

3. Read satellite data
######################
Once the satellite data is downloaded one can access and read the data for further use with **wavy**. Let's have a look at some examples in a python script.

In python L3-data can be read by importing the satellite_class, choosing a region of interest, the variable of interest (Hs or U), the satellite mission, which product should be used, and whether a time window should be used as well as a start and possibly an end date. This could look like:

.. code-block:: python3

    >>> from wavy.satellite_module import satellite_class as sc
    >>> nID = "cmems_L3_NRT"
    >>> name = "s3a"
    >>> sd = "2024-10-17 07"
    >>> ed = "2024-10-17 09"
    >>> sco = sc(sd=sd, nID=nID, name=name, ed=ed).populate()
    
This would result in a satellite_class object and a similar output message as::

    # ----- 
     ### Initializing satellite_class object ###
 
     Given kwargs:
    {'sd': '2024-10-17 07', 'nID': 'cmems_L3_NRT', 'name': 's3a', 'ed': '2024-10-17 09'}
 
     ### satellite_class object initialized ###
    # ----- 
     ### Read files and populate satellite_class object
     ## Find and list files ...
    path is None -> checking config file
    Object is iterable
    9 valid files found
    source template: /path/to/Moz_ws24_wavy/altimeter_data/L3/name/%Y/%m


    ....


     ## Summary:
    5238 footprints retrieved.
    Time used for retrieving data:
    0.3 seconds
 
     ### satellite_class object populated ###

    # ----- 

Investigating the satellite_object you will find something like::

    >>> sco.
    sco.apply_limits(                             sco.filter_main(
    sco.cfg                                       sco.filter_NIGP(
    sco.cleaner_blockQ(                           sco.filter_runmean(
    sco.cleaner_blockStd(                         sco.get_item_child(
    sco.compute_pulse_limited_footprint_radius()  sco.get_item_parent(
    sco.coords                                    sco.list_input_files(
    sco.crop_to_period(                           sco.meta
    sco.crop_to_poi(                              sco.name
    sco.crop_to_region(                           sco.nID
    sco.despike_blockQ(                           sco.pathlst
    sco.despike_blockStd(                         sco.poi
    sco.despike_GP(                               sco.populate(
    sco.despike_linearGAM(                        sco.quick_anim(
    sco.despike_NIGP(                             sco.quicklook(
    sco.distlim                                   sco.reader(
    sco.download(                                 sco.region
    sco.ed                                        sco.sd
    sco.filter                                    sco.slider_chunks(
    sco.filter_blockMean(                         sco.stdvarname
    sco.filter_distance_to_coast(                 sco.time_gap_chunks(
    sco.filter_footprint_land_interaction(        sco.twin
    sco.filter_footprint_radius(                  sco.units
    sco.filter_GP(                                sco.varalias
    sco.filter_lanczos(                           sco.varname
    sco.filter_landMask(                          sco.vars
    sco.filter_linearGAM(                         sco.write_to_nc(

With the retrieved data in sco.vars::

   >>> sco.vars
    Dimensions:  (time: 5238)
    Coordinates:
      * time     (time) datetime64[ns] 42kB 2024-10-17T06:30:00 ... 2024-10-17T09...
    Data variables:
        Hs       (time) float64 42kB 2.068 2.065 2.063 2.063 ... 2.406 2.386 2.374
        lons     (time) float64 42kB -138.7 -138.7 -138.7 ... -168.3 -168.3 -168.4
        lats     (time) float64 42kB 52.22 52.27 52.33 52.39 ... -25.42 -25.36 -25.3
    Attributes:
        title:    wavy dataset

Using the quicklook function you can quickly visualize the data you have retrieved::

   >>> sco.quicklook(ts=True) # for time series
   >>> sco.quicklook(m=True) # for a map
   >>> sco.quicklook(a=True) # for all


4. Define your own region
#########################

In wavy you can define your own region over which you want to gather satellite data. The region has to be defined in the region_cfg.yaml file. It can either be defined as a rectangular region, a polynom, a geojson format, or a model. If region is a model defined in model_specs.yaml, this will automatically be noticed and a model file will be loaded to cross-check the model domain with the satellite footprints. Let's define Mozambique as a new region:

.. code-block:: yaml

    Moz:
        llcrnrlon: 28.3
        llcrnrlat: -27.8
        urcrnrlon: 46
        urcrnrlat: -10

Now, we use this region to retrieve only data over this region.

.. code-block:: python3

    >>> from wavy.satellite_module import satellite_class as sc
    >>> nID = "cmems_L3_NRT"
    >>> name = "s3a"
    >>> sd = "2024-10-17 07"
    >>> ed = "2024-10-17 09"
    >>> sco = sc(sd=sd, nID=nID, name=name, ed=ed, region="Moz").populate()
    >>> sco.quicklook(m=True)
    
Another option is to define the region directly in the script:

.. code-block:: python3

    >>> from wavy.satellite_module import satellite_class as sc
    >>> region_dict = {'name': 'Moz',
    >>>                'region': {
    >>>                    'llcrnrlon': 28.3,
    >>>                    'llcrnrlat': -27.8,
    >>>                    'urcrnrlon': 46,
    >>>                    'urcrnrlat': -10}}
    >>> nID = "cmems_L3_NRT"
    >>> name = "s3a"
    >>> sd = "2024-10-17 07"
    >>> ed = "2024-10-17 09"
    >>> sco = sc(sd=sd, nID=nID, name=name, ed=ed).populate(region=region_dict)

You can adapt the window for the map as well as follows:

.. code-block:: python3

    >>> sco.quicklook(m=True, map_extent_llon=28.3, map_extent_ulon=46,
                      map_extent_llat=-27.8, map_extent_ulat=-10)
    

5. access/read model data
#########################

Model output can be accessed and read using the model_module module. The model_module config file model_cfg.yaml needs adjustments if you want to include a model that is not present as default. Given that the model output file you would like to read follows the cf-conventions and standard_names are unique, the minimum information you have to provide are usually:

.. code-block:: yaml

   modelname:
       vardef:
           Hs: 
           time: 
           lons: 
           lats: 
       wavy_input:
           fl_tmplt:
       reader: 
       misc:
        init_times: 
        init_step: 
        grid_date: 
        date_incr_unit: 
        date_incr: 


The variable aliases (left hand side below vardef) need to be specified in the variable_def.yaml. Basic variables are already defined. Adding your model output files to wavy means to add something like:

.. code-block:: yaml

    era5Moz:
        name:
        download:
        vardef:
            Hs: swh
            time: valid_time
            lons: longitude
            lats: latitude
        coords:
        wavy_input:
            src_tmplt: /path/to/Moz_ws24_wavy/data/
            fl_tmplt: era5_reanalysis_241016_241024_moz.nc
        reader: read_era
        collector:
        misc:
            init_times: [0]
            init_step: 24
            grid_date: 2024-10-23 00:00:00
            convention: meteorological
            date_incr_unit: h
            date_incr: 1
            proj4: "+proj=longlat +a=6367470 +e=0 +no_defs"

Now you can proceed to load your model in wavy. Start python and type:

.. code-block:: python3

    >>> from wavy.model_module import model_class as mc
    >>> nID = 'era5Moz' # default
    >>> varalias = 'Hs' # default
    >>> sd = "2024-10-17 07"
    >>> mco = mc(nID=nID,sd=sd).populate() # one time slice

Whenever the keyword "leadtime" is None, a best estimate is assumed and retrieved. In this case you are using reanalysis data, meaning that there is no leadtime to take into account. The output will be something like::

   >>> mco.
   mco.cfg                mco.filter             mco.meta               mco.quick_anim(        mco.stdvarname         
   mco.coords             mco.get_item_child(    mco.model              mco.quicklook(         mco.units              
   mco.crop_to_period(    mco.get_item_parent(   mco.nID                mco.reader(            mco.varalias           
   mco.distlim            mco.leadtime           mco.pathlst            mco.region             mco.varname            
   mco.ed                 mco.list_input_files(  mco.populate(          mco.sd                 mco.vars     


   >>> mco.vars
   
       Dimensions:  (time: 1, lats: 43, lons: 41)
    Coordinates:
      * lons     (lons) float64 328B 30.0 30.5 31.0 31.5 ... 48.5 49.0 49.5 50.0
      * lats     (lats) float64 344B -9.0 -9.5 -10.0 -10.5 ... -29.0 -29.5 -30.0
        number   int64 8B ...
      * time     (time) datetime64[ns] 8B 2024-10-16T23:00:00
        expver   <U4 16B ...
    Data variables:
        Hs       (time, lats, lons) float32 7kB nan nan nan ... 5.322 5.505 5.725
    Attributes:
        GRIB_centre:             ecmf
        GRIB_centreDescription:  European Centre for Medium-Range Weather Forecasts
        GRIB_subCentre:          0
        Conventions:             CF-1.7
        institution:             European Centre for Medium-Range Weather Forecasts
        history:                 2024-10-30T14:48 GRIB to CDM+CF via cfgrib-0.9.1...

For the model_class objects a quicklook function exists to depict a certain time step of what you loaded::

   >>> mco.quicklook(m=True) # for a map
   >>> mco.quicklook(a=True) # for a map

6. Collocating model and observations
#####################################
One main focus of **wavy** is to ease the collocation of observations and numerical wave models for the purpose of model validation. If you have available the necessary satellite data and model data you can proceed with collocation:

Collocation of satellite and wave model
****************************************

.. code-block:: python3

    >>> from wavy.satellite_module import satellite_class as sc
    >>> from wavy.collocation_module import collocation_class as cc

    >>> # retrieve the satellite data for the region
    >>> nID = "cmems_L3_NRT"
    >>> name = "s3a"
    >>> sd = "2024-10-17 07"
    >>> ed = "2024-10-17 09"
    >>> sco = sc(sd=sd, nID=nID, name=name, ed=ed, region='Moz').populate()
    >>> # collocate the model
    >>> model = 'era5Moz'
    >>> cco = cc(oco=sco, model=model, leadtime='best', distlim=6, twin=180)

*distlim* is the distance limit for collocation in *km* and date_incr is the time step increase in hours. One can also add a keyword for the collocation time window. The default is +-30min which is equivalent to adding *twin=30*. In this case ERA only had 6h time steps which makes it a bit more unlikely that satellite crossings and model time steps coincide. Increasing *twin* helps, however, it means we assume quasi-stationarity for this time period.

Using the quicklook function again (*cco.quicklook(a=True)*) will enable three plots this time, a time series plot (*ts=True*), a map plot (*m=True*), and a scatter plot (*sc=True*)

.. code-block:: python3
    
    cco.quicklook(a=True)

7. Validate the collocated time series
#######################################
Having collocated a quick validation can be performed using the validationmod. validation_specs.yaml can be adjusted.

.. code-block:: python3

   >>> val_dict = cco.validate_collocated_values()

    # ---
    Validation stats
    # ---
    Correlation Coefficient: 0.60
    Mean Absolute Difference: 0.54
    Root Mean Squared Difference: 0.63
    Normalized Root Mean Squared Difference: 0.20
    Debiased Root Mean Squared Difference: 0.59
    Bias: -0.22
    Normalized Bias: -0.07
    Scatter Index: 18.82
    Model Activity Ratio: 1.31
    Mean of Model: 2.91
    Mean of Observations: 3.13
    Number of Collocated Values: 6

The entire validation dictionary will then be in val_dict.


Download Satellite data
#######################

Preparations
------------

One strength of **wavy** is to ease obtaining and using satellite altimetry data. To get a flavor of how this can be done illustrated with cmems level 3 altimeter data please execute **wavy**'s command line download script:

.. code-block:: bash

   $ conda activate wavyopen
   $ wavyDownload --help


The help-message displayed would give you, among other information, the following options:


.. code-block:: bash

    cmems_L3_NRT:            
     s3a - Sentinel-3A            
     s3b - Sentinel-3B            
     j3 - Jason-3 (deprecated reference mission)
     c2 - Cryosat-2            
     al - SARAL/AltiKa            
     cfo - CFOSAT            
     h2b - HaiYang-2B            
     s6a - Sentinel-6A Michael Freilich (reference mission)
     swon - SWOT nadir altimeter
                
This means that for product cmems_L3_NRT you can choose among 9 satellite missions. In the world of **wavy**, cmems_L3_NRT is called the *nID* (name ID) and the individual missions are *names*. These are arbitrary names. Although they could be choosen freely by each user, it makes sense to keep them descriptive. As for cmems data, unfortunatly, most of other openly available satellite altimeter data is not accessible via thredds or similar options but needs to be downloaded from e.g. an FTP server. To do that you would need the credentials for CEDA, or for the AVISO cataloque as these are the main other sources that **wavy** currently exploits. Both CEDA and AVISO are possible to use as source in wavy but are not included in the *help* function for wavyDownload due to the large amount of choices.

**wavy** relies on the coprenicusmarine toolbox for CMEMS products. In case of remote access via FTP **wavy** needs you to store the respective usernames and passwords in your local .netrc file. This could look like:

.. code::

   machine ftp.ceda.ac.uk    login {USER}  password {PASSWORD}
   machine ftp-access.aviso.altimetry.fr    login {USER}  password {PASSWORD}

In case of using the copernicusmarine toolbox the user needs to make sure that the CMEMS credentials are available by storing them in the .env file, e.g. located in the project directory:

.. code::

   COPERNICUSMARINE_SERVICE_USERNAME=YOUR_COPERNICUS_USERNAME
   COPERNICUSMARINE_SERVICE_PASSWORD=YOUR_COPERNICUS_PASSWORD

Download using config files
===========================

Ammending config files
----------------------
In a validation context, especially operational, download operation needs to be performed many times and it makes sense to adjust the satellite_cfg.yaml file to your needs. Assuming you established a project directory and therein config directory you first establish a minimal satellite_cfg.yaml file:

.. code-block:: bash

   $ cd ~/my_wavy_project/config
   $ conda activate wavyopen
   $ wavyCFG --path ~/my_wavy_project/config/. --f satellite_cfg.yaml --t minimal

It looks like this:

.. code-block:: yaml

   cmems_L3_NRT:
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
    download:
        copernicus:
            dataset_id: cmems_obs-wave_glo_phy-swh_nrt_name-l3_PT1S
            trgt_tmplt:
            path_date_incr_unit: 'm'
            path_date_incr: 1
            strsub: ["name"]
            server: "nrt.cmems-du.eu"
            time_incr: 'h'
    wavy_input:
        src_tmplt:
        fl_tmplt:
        strsub: ["name"]
        path_date_incr_unit: 'm'
        path_date_incr: 1
    reader: read_local_ncfiles
    collector: get_remote_files_copernicusmarine
    vardef:
        Hs: VAVH
        U: WIND_SPEED
    coords:
    misc:
        processing_level:
        provider:
        obs_type:


Now, you can ammend it to your needs. Here is an explanation of the most important variables.

* cmems_L3_NRT - this is the name ID (*nID*) which often refers to a product which has multiple subproducts that are called *name*

* name - attributes names used in wavy (left hand side) to names used in the product (right hand side)

* download - consists of two download types: copernicus marine toolbox or FTP. Both need specifications such as the target path (trgt_tmplt) and for dynamic paths a substituting list (strbsub). In the example above this list contains only 'name' which means that the string name in the dataset_id *cmems_obs-wave_glo_phy-swh_nrt_name-l3_PT1S* will be replaced by any of the names specified above. The same is valid for the trgt_tmplt.

*  wavy_input - specifies, among other things, the source path (src_tmplt) from where wavy should find the data. Typically this is the same as trgt_tmplt described above, but it can also be something different. The strbsub under, works in the same way as explained above.

* vardef - here it is important to specify the exact netcdf name. The variable names on the left hand side are the names **wavy** is using internally and the right hand side are the names used in the files. Any variable name on the left hand side must be described in the variable_def.yaml file.

Since most things are already defined for this product, the only thing we need to change is the path to where we should download. An example is given below:

.. code-block:: yaml

   cmems_L3_NRT:
       download:
           copernicus:
               dataset_id: cmems_obs-wave_glo_phy-swh_nrt_name-l3_PT1S
               trgt_tmplt: /chosen/path/to/satellite/data/L3/name/%Y/%m


The str "name" in your path_template will be replaced by the satellite mission that you download because it was defined in the strsub list. So for Sentinel-3a the final path for your downloaded files will be automatically /chosen/path/to/satellite/data/L3/s3a with subfolders on year and month.

You can now proceed and download like:

.. code-block:: bash

   $ wavyDownload --nID cmems_L3_NRT --name s3a --sd 2025010100 --ed 2025011000

You can find the downloaded files in your chosen download directory.

You can also download altimeter data directly in the python script with the following lines. 

.. code-block:: python3

   >>> from wavy.satellite_module import satellite_class as sc
   >>> nID = 'cmems_L3_NRT'
   >>> name = 's3a'
   >>> sd = '2023-11-10 00'
   >>> ed = '2023-11-10 10'
   >>> # Initialize sc object
   >>> sco = sc(sd=sd,ed=ed,nID=nID,name=name)
   >>> # Download the data to a chosen directory
   >>> path = '/chosen/path/to/satellite/data/L3/s3a'
   >>> sco.download(path=path)

In case of ftp downloads the config setup is similar but you have to make the adjustments under the ftp section:

.. code-block:: yaml

   download:
       ftp: # downloading method
           src_tmplt: "/path/to/remote/dir/%Y/%m"
           trgt_tmplt: /chosen/path/to/satellite/data/L3/name/%Y/%m
           strsub: ['name']

Also the collector needs to be a suitable one like e.g. the one used for CCI files:

.. code-block:: yaml

   collector: get_remote_files_ftp

With ftp, parallel python can be used with a keyword specifying the number of processes, e.g.:

.. code-block:: python3

   >>> sco.download(nproc=4, path=path)

Download using without preparing config files
=============================================

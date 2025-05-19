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
                
This means that for product cmems_L3_NRT you can choose among 9 satellite missions. In the world of **wavy**, cmems_L3_NRT is called the *nID* (name ID) and the individual missions are *names*. These are arbitrary names. Although they could be choosen freely by each user, it makes sense to keep them descriptive. As for cmems data, unfortunatly, most of other openly available satellite altimeter data is not accessible via thredds or similar options but needs to be downloaded from e.g. a FTP server. To do that you would need the credentials for CEDA, or for the AVISO cataloque as these are the main other sources that **wavy** currently exploits.

**wavy** relies on the coprenicusmarine toolbox for CMEMS products. In case of remote access via FTP **wavy** needs you to store the respective usernames and passwords in your local .netrc file. This could look like:

.. code::

   machine ftp.ceda.ac.uk    login {USER}  password {PASSWORD}
   machine ftp-access.aviso.altimetry.fr    login {USER}  password {PASSWORD}

In case of using the copernicusmarine toolbox the user needs to make sure that the CMEMS credentials are available by storing them in your .env file in your e.g. project directory:

.. code::

   COPERNICUSMARINE_SERVICE_USERNAME=YOUR_COPERNICUS_USERNAME
   COPERNICUSMARINE_SERVICE_PASSWORD=YOUR_COPERNICUS_PASSWORD

Ammending config files
----------------------
In a validation context, especially operational, download operation needs to be performed many times and it makes sense to adjust the satellite_cfg.yaml file to your needs. Assuming you established a project directory and therein config directory you first establish a minimal satellite_cfg.yaml file:

.. code-block:: bash

   $ cd ~/my_wavy_project/config
   $ conda activate wavyopen
   $ wavyCFG --path ~/my_wavy_project/config/. --f satellite_cfg.yaml --t minimal

Now, you can ammend it to your needs.

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



Now, prepare your **wavy** environment with providing the directories for satellite data and model data. Add your path for satellite data here demonstrated for CMEMS using the copercnicusmarine toolbox, indicating the path of your choice where you want your data to be stored:

.. code-block:: yaml

   cmems_L3_NRT:
       download:
           copernicus:
               dataset_id: cmems_obs-wave_glo_phy-swh_nrt_name-l3_PT1S
               trgt_tmplt: /chosen/path/to/satellite/data/L3/name/%Y/%m


There exists also something called strsub which defines strings that are o substituted. In this case some are predefined as:

.. code-block:: yaml

   strsub: ['name']

The str "name" in your path_template will be replaced by the satellite mission that you download. So for Sentinel-3a the final path for your downloaded files will be automatically /chosen/path/to/satellite/data/L3/s3a with subfolders on year and month.

You can proceed now and download CMEMS NRT L3 data using the wavyDownload.py script:

.. code-block:: bash

   $ cd ~/wavy/wavy/apps

To get help check ...

.. code-block:: bash

   $ ./wavyDownload.py -h

... or download some satellite altimeter data:

.. code-block:: bash

   $ ./wavyDownload.py --name s3a --sd 2020110100 --ed 2020111000 --nID cmems_L3_NRT

You can find the downloaded files in your chosen download directory.

Similarily one can download L2P and L3 multi-mission altimetry data from the CEDA Climate Change Initiative. This spans a long time period from 1991 to 2018 and enables climate related research and wave model hindcast validation.

.. code-block:: bash

   $ ./wavyDownload.py -sat multi -sd 2017112000 -ed 2017112100 -product cci_L3
   
You can also download altimeter data directly from python with the following lines. 

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

In case of ftp downloads the config setup is similar but you have to make adjustments under the ftp section:

.. code-block:: yaml

   download:
       ftp: # downloading method
           src_tmplt: "/path/to/remote/dir/%Y/%m"
           trgt_tmplt: /chosen/path/to/satellite/data/L3/name/%Y/%m
           strsub: ['name']

With ftp, parallel python can be used with a keyword specifying the number of processes, e.g.:

.. code-block:: python3

   >>> sco.download(nproc=4, path=path)




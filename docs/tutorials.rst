Tutorials
==========
The following tutorials describe how to setup the wavy config files, to perform collocation and validation tasks. **wavy** is intended to be fleksibel such that customization can be achieved with minimal changes in the code. The **wavy** config files build the fundament for this approach and serve a similar purpose as namelist files often used for runtime changes in numerical modelling.

In general, executable files usually have help function which can be read using e.g.:

.. code-block:: bash

   $ ./{Executable}.py -h

e.g.:

.. code-block:: bash

   $ cd ~/wavy/apps/standalone
   $ ./wavyDownload.py -h

1. **wavy** config files, a brief overview.
###########################################
**wavy** obtains custom information from config files. There are default versions which can be adjusted to user needs. The following config files exist as for now:


.. code-block:: bash

   $ ls
   collocation_specs.yaml.default  region_specs.yaml.default
   d22_var_dicts.yaml.default      satellite_specs.yaml.default
   insitu_specs.yaml.default       validation_specs.yaml.default
   model_specs.yaml.default        variable_info.yaml.default
   quicklook_specs.yaml.default

**wavy** browses the directory structure as follows:

    * check if env 'WAVY_CONFIG' is set or specified in .env
    * check if a config folder exists using xdg
    * fall back on default files within the package

It is probably easiest to start off copying the default config files you would like to edit to a directory of your choice. Then remove the ".default" extension. For simplicity, let's assume that you have your custom config files in *~/wavy/config/*.

2. download L3 satellite altimetry data
#######################################

L3 satellite data is obtained from Copernicus with the product identifier WAVE_GLO_WAV_L3_SWH_NRT_OBSERVATIONS_014_001. User credentials are required for this task. So before you can start you have to get a Copernicus account (free of costs).
Prepare access to Copernicus products. Enter your account credentials into the .netrc-file. Your .netrc should look something like:

.. code::

   machine nrt.cmems-du.eu    login {USER}  password {PASSWORD}

Prepare your **wavy** environment with providing the directories for satellite data and model data. There are multiple config files but we only need to worry about a few for now. Explore the config file for satellites like this:

.. code-block:: bash

   $ cd ~/wavy/config
   $ vim satellite_specs.yaml

Add your path for satellite data here under cmems

.. code-block:: yaml

   cmems_L3:
      dst:
         path_template: /home/patrikb/tmp_altimeter/L3/mission

You can proceed now and download L3 data using the wavyDownload.py script:

.. code-block:: bash

   $ cd ~/wavy/apps/standalone

To get help check ...

.. code-block:: bash

   $ ./wavyDownload.py -h

... then download some satellite altimeter data:

.. code-block:: bash

   $ ./wavyDownload.py -sat s3a -sd 2020110100 -ed 2020111000 -product cmems_L3

You can find the downloaded files in your chosen download directory.

3. download L2P and L3 CCI multi-mission satellite altimetry data
#################################################################
Similarily one can download L2P and L3 multi-mission altimetry data from the CEDA Climate Change Initiative. This spans a long time period from 1991 to 2018 and enables climate related research and wave model hindcast validation.

For instance for Jason-3:

.. code-block:: bash

   $ ./wavyDownload.py -sat j3 -sd 2017112000 -ed 2017112100 -product cci_L2P


Or for instance for a multi-mission file:

.. code-block:: bash

   $ ./wavyDownload.py -sat multi -sd 2017112000 -ed 2017112100 -product cci_L3


4. download L2 stallite altimetry data
######################################

.. note::

   There are currently problems with L2 from eumetsat/colhub which
   will be fixed again hopefully soon.

L2 satellite data are obtained from eumetsat and colhub using the SentinelAPI. This requires user credentials for eumetsat and colhub, which are free of costs as well.
Enter your account credentials into the .netrc-file as you did for the L3 data. Your .netrc should have included the following:

.. code::

   machine https://colhub.met.no/ login {USER} password {PASSWORD}
   machine https://coda.eumetsat.int/search login {USER} password {PASSWORD}

Ammend the satellite config file for L2 data and add the download directory of your choice like:

.. code-block:: yaml

   eumetsat_L2:
      L2:
         dst:
             path_template: /home/patrikb/tmp_altimeter/L2/mission

As you can see, this is customized to my username patrikb. Adjust this and continue with downloading some satellite altimeter data:

.. code-block:: bash

   $ ./wavyDownload.py -sat s3a -sd 2020110100 -ed 2020111000 -product eumetsat_L2

4. read satellite data
######################
Once the satellite data is downloaded one can access and read the data for further use with **wavy** or other software.

L3 data from cmems
******************

In python L3 data can be read by importing the satellite_class, choosing a region of interest, the variable of interest (Hs or U), the satellite mission, which product should be used, and whether a time window should be used as well as a start and possibly an end date. This could look like:

.. code-block:: python3

   >>> from wavy.satmod import satellite_class as sc
   >>> region = 'NorwegianSea'
   >>> varalias = 'Hs' # default
   >>> mission = 's3a' # default
   >>> product = 'cmems_L3' # default
   >>> twin = 30 # default
   >>> sd = "2020-11-1" # can also be datetime object
   >>> ed = "2020-11-2" # not necessary if twin is specified
   >>> sco = sc(sdate=sd,edate=ed,region=region)

This would result in a satellite_class object and the following output message::

   >>> sco = sc(sdate=sd,edate=ed,region=region)
   # -----
   ### Initializing satellite_class object ###

   Parsing date
   Translate to datetime
   Parsing date
   Translate to datetime
   Requested time frame: 2020-11-01 00:00:00 - 2020-11-02 00:00:00
   Chosen time window is: 30 min
   No download initialized, checking local files

    ## Find files ...
   path_local is None -> checking config file
   /home/patrikb/tmp_altimeter/L3/s3a/2020/10
   /home/patrikb/tmp_altimeter/L3/s3a/2020/11
   26 valid files found

    ## Read files ...
   Get filevarname for
   stdvarname: sea_surface_wave_significant_height
   varalias: Hs
   !!! standard_name:  sea_surface_wave_significant_height  is not unique !!!
   The following variables have the same standard_name:
    ['VAVH', 'VAVH_UNFILTERED']
   Searching *_specs.yaml config file for definition
   Variable defined in *_specs.yaml is:
   Hs = VAVH
   100%|███████████████████████████████████████████| 26/26 [00:00<00:00, 91.45it/s]
   Concatenate ...
   ... done concatenating
   Total:  46661  footprints found
   Apply region mask
   Specified region: NorwegianSea
    --> Bounded by polygon:
   lons: [5.1, -0.8, -6.6, -9.6, -8.6, -7.5, 1.7, 8.5, 7.2, 16.8, 18.7, 22.6, 18.4, 14.7, 11.7, 5.1]
   lats: [62.1, 62.3, 63.2, 64.7, 68.5, 71.1, 72.6, 74.0, 76.9, 76.3, 74.5, 70.2, 68.3, 66.0, 64.1, 62.1]
   Values found for chosen region and time frame.
   Region mask applied
   For chosen region and time:  351 footprints found

   ## Summary:
   Time used for retrieving satellite data: 0.34 seconds
   Satellite object initialized including 351 footprints.
   # -----


Investigating the satellite_object you will find something like::

   >>> sco.
   sco.edate             sco.product           sco.units
   sco.get_item_child(   sco.provider          sco.varalias
   sco.get_item_parent(  sco.quicklook(        sco.varname
   sco.mission           sco.region            sco.vars
   sco.obstype           sco.sdate             sco.write_to_nc(
   sco.path_local        sco.stdvarname
   sco.processing_level  sco.twin

With the retrieved variables in sa_obj.vars::

   >>> sco.vars.keys()
   dict_keys(['sea_surface_wave_significant_height', 'time', 'time_unit', 'latitude', 'longitude', 'datetime', 'meta'])

Read pure L2 satellite data from eumetsat
*****************************************

.. note::

   There are currently problems with L2 from eumetsat/colhub which
   will be fixed again hopefully soon.

.. code-block:: python3

   >>> from wavy.satmod import satellite_class as sc
   >>> sd = "2020-11-1 12"
   >>> ed = "2020-11-1 12"
   >>> region = 'mwam4' # default
   >>> mission = 's3a' # default
   >>> twin = 30 # default
   >>> varalias = 'Hs' # default

   >>> sco = sc(sd,edate=ed,product="eumetsat_L2")

Retrieve pure L2 data and compare against L3
********************************************

.. note::

   There are currently problems with L2 from eumetsat/colhub which
   will be fixed again hopefully soon.

Having downloaded the altimetry data, you can do:

.. code-block:: python3

   >>> # imports
   >>> from wavy.satmod import satellite_class as sc

   >>> # settings
   >>> sd = "2020-11-1 12"
   >>> ed = "2020-11-1 12"
   >>> region = 'NorwegianSea'
   >>> mission = 's3a' # default
   >>> varalias = 'Hs' # default
   >>> twin = 30 # default

   >>> # retrievals
   >>> sco_e = sc(sd,edate=ed,region=region,product='eumetsat_L2')
   >>> sco_c = sc(sd,edate=ed,region=region,product='cmems_L3')

   >>> # plotting
   >>> import matplotlib.pyplot as plt
   >>> stdname = sco_e.stdvarname
   >>> fig = plt.figure(figsize=(9,3.5))
   >>> ax = fig.add_subplot(111)
   >>> ax.plot(sco_e.vars['datetime'],sco_e.vars[stdname],'r.',label='L2 eumetsat')
   >>> ax.plot(sco_c.vars['datetime'],sco_c.vars[stdname],'k.',label='L3 cmems')
   >>> plt.legend(loc='upper left')
   >>> plt.ylabel('Hs [m]')
   >>> plt.show()

This yields the following figure:

.. image:: ./docs_fig_L2_vs_L3.png
   :scale: 80

Appy basic filters to raw L2 data
*********************************

.. code-block:: python3

   >>> from wavy.satmod import satellite_class as sc
   >>> import matplotlib.pyplot as plt

   >>> sd = "2020-11-1 12"
   >>> ed = "2020-11-1 12"
   >>> region = 'mwam4' # default
   >>> mission = 's3a' # default
   >>> twin = 30 # default

   >>> # landmask filter
   >>> sco_lm = sc(sd,edate=ed,product='eumetsat_L2',land_mask=True,filterData=True)

.. note::

   More examples with filters are coming soon ...

5. access/read model data
#########################
Model output can be accessed and read using the modelmod module. The modelmod config file model_specs.yaml needs adjustments if you want to include a model that is not present as default. Given that the model output file you would like to read follows the cf-conventions and standard_names are unique, the minimum information you have to provide are usually:

.. code-block:: yaml

   modelname:
       path_template:
       file_template:
       init_times: []
       init_step:

Often there are ambiguities due to the multiple usage of standard_names. Any such problem can be solved here in the config-file by adding the specified variable name like:

.. code-block:: yaml

    vardef:
        Hs: VHM0
        time: time
        lons: lon
        lats: lat

The variable aliases (left hand side) need to be specified in the variable_info.yaml. Basic variables are already defined. All specs listed here are also used when **wavy** writes the retrieved values to netcdf.

.. code-block:: python3

   >>> from wavy.modelmod import model_class as mc
   >>> model = 'mwam4' # default
   >>> varalias = 'Hs' # default
   >>> sd = "2020-11-1"
   >>> sd = "2020-11-2"
   >>> mco = mc(sdate=sd) # one time slice
   >>> mco_p = mc(sdate=sd,edate=ed) # time period
   >>> mco_lt = mc(sdate=sd,leadtime=12) # time slice with lead time

Whenever the keyword "leadtime" is None, a best estimate is assumed and retrieved. The output will be something like::

   >>> mco = mc(sdate=sd)

   >>> mco.
   mco.edate             mco.leadtime          mco.units
   mco.fc_date           mco.model             mco.varalias
   mco.filestr           mco.quicklook(        mco.varname
   mco.get_item_child(   mco.sdate             mco.vars
   mco.get_item_parent(  mco.stdvarname

   >>> mco.vars.keys()
   dict_keys(['longitude', 'latitude', 'time', 'datetime', 'time_unit', 'sea_surface_wave_significant_height', 'meta', 'leadtime'])

.. note::

   Even though it is possible to access a time period, **wavy** is not yet optimized to do so and the process will be slow. The reason, being the ambiguous use of lead times, will be improved in future versions.

6. read in-situ observations (.d22 and netcdf/thredds)
######################################################

Currently two data types can be read .d22-files and netcdf-files.

read .d22 files
***************

.. note::

   The .d22-files used in these examples and specified here are only available to MET Norway stuff.

.d22-files can be read in by adjusting d22_var_dicts.yaml config file. Currently, there are wave related variables included. Other variables like wind are about to be included. Another config-file that needs adjustment is the insitu_specs.yaml. There you need to define specs related to the in-situ observation of choice as well as path and filename. A call for the retrieval of an in-situ time series could be like:

.. code-block:: python3

   >>> from wavy.insitumod import insitu_class as ic
   >>> varalias = 'Hs' # default
   >>> sd = "2020-1-1"
   >>> ed = "2020-1-5"
   >>> nID = 'ekofiskL'
   >>> sensor = 'waverider'
   >>> ico = ic(nID,sensor,sd,ed)

In contrast to the L3 satellite time series, in-situ time series are not filtered or underwent rigorous outlier detection. There are various operations that can be performed to massage the time series as you wish.It is in particular interesting to remove double reported values, which is often the case. This is done with setting unique=True.

.. code-block:: python3

   >>> ico = ic(nID,sensor,sd,ed,unique=True)

read .nc-files
**************

.. code-block:: python3

   >>> from wavy.insitumod import insitu_class as ic
   >>> varalias = 'Hs' # default
   >>> sd = "2020-1-1"
   >>> ed = "2020-1-5"
   >>> nID = 'D_Breisundet_wave'
   >>> sensor = 'wavescan'
   >>> ico = ic(nID,sensor,sd,ed)

Additionally, outliers can be removed, missing data can be treated, and super-observations can be formed. Below is a example:

.. code-block:: python3

   >>> # blockMean filter
   >>> ico_bm = ic(nID,sensor,sd,ed,unique=True,priorOp='square',postOp='root',smoother='blockMean',stwin=3,etwin=3,date_incr=1,filterData=True)

Now, let's check how this could look like:

.. code-block:: python3

   >>> import matplotlib.pyplot as plt
   >>> fig = plt.figure(figsize=(9,3.5))
   >>> ax = fig.add_subplot(111)
   >>> ax.plot(ico.vars['datetime'],ico.vars[ico.stdvarname],'ko',label='raw')
   >>> ax.plot(ico_bm.vars['datetime'],ico_bm.vars[ico.stdvarname],'r-',label='hourly blockMean')
   >>> plt.legend(loc='upper left')
   >>> plt.ylabel('Hs [m]')
   >>> plt.show()

.. image:: ./docs_fig_ts_insitu.png
   :scale: 80

7. collocating model and observations
#####################################
One of the main focus of **wavy** is to ease the collocation of observations and numerical wave models for the purpose of model validation. For this purpose there is the config-file collocation_specs.yaml where you can specify the name and path for the collocation file to be dumped if you wish to save them.

Collocation of satellite and wave model
****************************************

.. code-block:: python3

   >>> from wavy.satmod import satellite_class as sc
   >>> from wavy.collocmod import collocation_class as cc

   >>> model = 'mwam4' # default
   >>> mission = 's3a' # default
   >>> varalias = 'Hs' # default
   >>> sd = "2020-11-1 12"
   >>> sco = sc(sdate=sd,region=model,sat=mission,varalias=varalias)
   >>> cco = cc(model=model,obs_obj_in=sco,distlim=6,date_incr=1)

   >>> # plotting
   >>> import matplotlib.pyplot as plt
   >>> fig = plt.figure(figsize=(9,3.5))
   >>> ax = fig.add_subplot(111)
   >>> ax.plot(cco.vars['datetime'],cco.vars['obs_values'],color='gray',marker='o',linestyle='None',alpha=.4,label='obs')
   >>> ax.plot(cco.vars['datetime'],cco.vars['model_values'],'b.',label='model',lw=2)
   >>> plt.legend(loc='upper left')
   >>> plt.ylabel('Hs [m]')
   >>> plt.show()

.. image:: ./docs_fig_ts_sat.png
   :scale: 80

This can also be done for a time period:

.. code-block:: python3

   >>> sd = "2020-11-1"
   >>> ed = "2020-11-2"
   >>> sco = sc(sdate=sd,edate=ed,region=model,mission=mission,varalias=varalias)
   >>> cco = cc(model=model,obs_obj_in=sco,distlim=6,date_incr=1)

Collocation of in-situ data and wave model
******************************************

The following example may take a few minutes.

.. code-block:: python3

   >>> # imports
   >>> from wavy.insitumod import insitu_class as ic
   >>> from wavy.collocmod import collocation_class as cc

   >>> # settings
   >>> model = 'mwam4' # default
   >>> varalias = 'Hs' # default
   >>> sd = "2020-1-1 01"
   >>> ed = "2020-1-4 00"
   >>> nID = 'ekofiskL'
   >>> sensor = 'waverider'

   >>> # retrievals
   >>> ico_gam = ic(nID,sensor,sd,ed,smoother='linearGAM',cleaner='linearGAM',date_incr=1./6.,unique=True,filterData=True)
   >>> ico_raw = ic(nID,sensor,sd,ed)

   >>> # collocation
   >>> cco_gam = cc(model=model,obs_obj_in=ico_gam,distlim=6,date_incr=1)
   >>> cco_raw = cc(model=model,obs_obj_in=ico_raw,distlim=6,date_incr=1)

Let's plot the results:

.. code-block:: python3

   >>> import matplotlib.pyplot as plt
   >>> stdname = ico_raw.stdvarname

   >>> fig = plt.figure(figsize=(9,3.5))
   >>> ax = fig.add_subplot(111)
   >>> ax.plot(ico_raw.vars['datetime'],ico_raw.vars[stdname],color='gray',marker='o',label='raw',linestyle='None',alpha=.4)
   >>> ax.plot(cco_raw.vars['datetime'],cco_raw.vars['obs_values'],'ko',label='collocated obs')
   >>> ax.plot(ico_gam.vars['datetime'],ico_gam.vars[stdname],'b-',label='gam',lw=2)
   >>> ax.plot(cco_gam.vars['datetime'],cco_gam.vars['model_values'],'r-',label='mwam4',lw=2)
   >>> plt.legend(loc='upper left')
   >>> plt.ylabel('Hs [m]')
   >>> plt.show()

.. image:: ./docs_fig_col_insitu.png
   :scale: 80

8. dump collocation ts to a netcdf file
#######################################
The collocation results can now be dumped to a netcdf file. The path and filename can be entered as keywords but also predefined config settings can be used from collocation_specs.yaml:

.. code-block:: python3

   >>> cco_raw.write_to_nc()

9. validate the collocated time series
#######################################
Having collocated a quick validation can be performed using the validationmod. validation_specs.yaml can be adjusted.

.. code-block:: python3

   >>> val_dict = cco_raw.validate_collocated_values()

   # ---
   Validation stats
   # ---
   Correlation Coefficient: 0.95
   Mean Absolute Difference: 0.22
   Root Mean Squared Difference: 0.27
   Normalized Root Mean Squared Difference: 0.08
   Debiased Root Mean Squared Difference: 0.24
   Bias: -0.13
   Normalized Bias: -0.04
   Scatter Index: 8.05
   Mean of Model: 3.02
   Mean of Observations: 3.14
   Number of Collocated Values: 72

The entire validation dictionary will then be in val_dict.

10. quick look examples
#######################
The script "wavyQuick.py" is designed to provide quick and easy access to information regarding satellite coverage and basic validation. Checkout the help:

.. code-block:: bash

   $ cd ~/wavy/apps/standalone
   $ ./wavyQuick.py -h

Browsing for satellite data of a given satellite mission and show footprints on map for a given time step and region:

For a model domain, here mwam4

.. code-block:: bash

   $ ./wavyQuick.py -sat s3a -reg mwam4 -sd 2020110112 --show

.. image:: ./docs_fig_sat_quicklook_001.png
   :scale: 25

or for a user-defined polygon

.. code-block:: bash

   $ ./wavyQuick.py -sat s3a -reg NorwegianSea -sd 2020110112 --show

.. image:: ./docs_fig_sat_quicklook_002.png
   :scale: 25

Browsing for satellite data and show footprints on map for time period would be the same approach simply adding an ending date:

.. code-block:: bash

   $ ./wavyQuick.py -sat s3a -reg NorwegianSea -sd 2020110100 -ed 2020110300 --show

.. image:: ./docs_fig_sat_quicklook_003.png
   :scale: 25

The same could be done choosing 10m wind speed instead of significant wave height:

.. code-block:: bash

   $ ./wavyQuick.py -var U -sat s3a -reg NorwegianSea -sd 2020110100 -ed 2020110300 --show

.. image:: ./docs_fig_sat_quicklook_004.png
   :scale: 25

The -sat argument can also be a list of satellites (adding the -l argument) or simply all available satellites:

.. code-block:: bash

   $ ./wavyQuick.py -sat list -l s3a,s3b,al -mod mwam4 -reg mwam4 -sd 2020110112 -lt 30 -twin 30 --col --show
   $ ./wavyQuick.py -sat all -mod mwam4 -reg mwam4 -sd 2020110112 -lt 30 -twin 30 --col --show

Now, dump the satellite data to a netcdf-file for later use:

.. code-block:: bash

   $ ./wavyQuick.py -sat s3a -reg mwam4 -sd 2020110100 -ed 2020110300 -dump /home/patrikb/tmp_altimeter/quickdump/test.nc

Browse for satellite data, collocate with wave model output and show footprints and model output for one time step and a given lead time (-lt 0) and time constraint (-twin 30):

.. code-block:: bash

   $ ./wavyQuick.py -sat s3a -reg NorwegianSea -mod mwam4 -sd 2020110112 -lt 0 -twin 30 --col --show

This results in a validation summary based on the collocated values:

.. code::

   # ---
   Validation stats
   # ---
   Correlation Coefficient: 0.95
   Mean Absolute Difference: 0.62
   Root Mean Squared Difference: 0.70
   Normalized Root Mean Squared Difference: 0.13
   Debiased Root Mean Squared Difference: 0.67
   Bias: 0.22
   Normalized Bias: 0.04
   Scatter Index: 12.71
   Mean of Model: 5.26
   Mean of Observations: 5.04
   Number of Collocated Values: 237

And of course the figure:

.. image:: ./docs_fig_sat_quicklook_005.png
   :scale: 40

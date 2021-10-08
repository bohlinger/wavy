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
**wavy** obtains custom information from config files. There are default versions which can be adjusted to user needs. The following config files exist as per now:


.. code-block:: bash

   $ ls
   collocation_specs.yaml.default  region_specs.yaml.default
   d22_var_dicts.yaml.default      satellite_specs.yaml.default
   insitu_specs.yaml.default       validation_specs.yaml.default
   model_specs.yaml.default        variable_info.yaml.default
   quicklook_specs.yaml.default

Another config file containing customized information to adjust the automated plotting routines is in the making.

**wavy** browses the directory structure as follows:

    * check if env 'WAVY_CONFIG' is set or specified in .env
    * check if a config folder exists using xdg
    * fall back on default files within the package

For simplicity, let's assume that we have our custom config files in *~/wavy/config*.

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

   cmems:
      L3:
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

   $ ./wavyDownload.py -sat s3a -sd 2020110100 -ed 2020111000 -prov cmems -lev L3

You can find the downloaded files in your chosen download directory.

3. download L2 stallite altimetry data
######################################
L2 satellite data are obtained from eumetsat and colhub using the SentinelAPI. This requires user credentials for eumetsat and colhub, which are free of costs as well.
Enter your account credentials into the .netrc-file as you did for the L3 data. Your .netrc should have included the following:

.. code::

   machine https://colhub.met.no/ login {USER} password {PASSWORD}
   machine https://coda.eumetsat.int/search login {USER} password {PASSWORD}

Ammend the satellite config file for L2 data and add the download directory of your choice like:

.. code-block:: yaml

   eumetsat:
      L2:
         dst:
             path_template: /home/patrikb/tmp_altimeter/L2/mission

As you can see, this is customized to my username patrikb. Adjust this and continue with downloading some satellite altimeter data:

.. code-block:: bash

   $ ./wavyDownload.py -sat s3a -sd 2020110100 -ed 2020111000 -prov eumetsat -lev L2

4. read satellite data
######################
Once the satellite data is downloaded one can access and read the data for further use in python.

L3 data from cmems
******************

L3 data can be read like:

.. code-block:: python3

   >>> from datetime import datetime
   >>> from wavy.satmod import satellite_class as sc
   >>> region = 'NorwegianSea'
   >>> varalias = 'Hs' # default
   >>> mission = 's3a' # default
   >>> provider = 'cmems' # default
   >>> level = 'L3' # default
   >>> twin = 30 # default
   >>> sd = datetime(2020,11,1)
   >>> ed = datetime(2020,11,2)
   >>> sco = sc(sdate=sd,edate=ed,region=region)

This would results in a satellite_class object and the following output message::

   >>> sco = sc(sdate=sd,edate=ed,region=region)
   Total:  148425  footprints found
   In chosen time period:  46661  footprints found
   Specified region: NorwegianSea
    --> Bounded by polygon:
   lons: [5.1, -0.8, -6.6, -9.6, -8.6, -7.5, 1.7, 8.5, 7.2, 16.8, 18.7, 22.6, 18.4, 14.7, 11.7, 5.1]
   lats: [62.1, 62.3, 63.2, 64.7, 68.5, 71.1, 72.6, 74.0, 76.9, 76.3, 74.5, 70.2, 68.3, 66.0, 64.1, 62.1]
   Values found for chosen region and time frame.
   For chosen region and time:  351 footprints found
   Time used for retrieving satellite data: 2.22 seconds
   Satellite object initialized including 351 footprints.

Investigating the satellite_object you will find something like::


   >>> sco.
   sco.edate             sco.path_local        sco.twin
   sco.get_item_child(   sco.provider          sco.varalias
   sco.get_item_parent(  sco.region            sco.varname
   sco.mission           sco.sdate             sco.vars
   sco.obstype           sco.stdvarname        sco.write_to_nc(

With the retrieved variables in sa_obj.vars::

   >>> sco.vars.keys()
   dict_keys(['time', 'latitude', 'longitude', 'sea_surface_wave_significant_height', 'time_unit', 'datetime', 'meta'])

Read pure L2 satellite data from eumetsat
*****************************************

.. code-block:: python3

   >>> from wavy.satmod import satellite_class as sc
   >>> from datetime import datetime
   >>> sd = datetime(2020,11,1,12)
   >>> ed = datetime(2020,11,1,12)
   >>> region = 'mwam4' # default
   >>> mission = 's3a' # default
   >>> level = 'L2'
   >>> twin = 30 # default
   >>> varalias = 'Hs' # default

   >>> sco = sc(sd,edate=ed,provider='eumetsat',level=level)

Retrieve pure L2 data and compare against L3
********************************************

Having downloaded the altimetry data, you can do:

.. code-block:: python3

   >>> # imports
   >>> from wavy.satmod import satellite_class as sc
   >>> from datetime import datetime

   >>> # settings
   >>> sd = datetime(2020,11,1,12)
   >>> ed = datetime(2020,11,1,12)
   >>> region = 'NorwegianSea'
   >>> mission = 's3a' # default
   >>> varalias = 'Hs' # default
   >>> twin = 30 # default

   >>> # retrievals
   >>> sco_e = sc(sd,edate=ed,region=region,provider='eumetsat',level='L2')
   >>> sco_c = sc(sd,edate=ed,region=region,provider='cmems')

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
   >>> from datetime import datetime
   >>> import matplotlib.pyplot as plt

   >>> sd = datetime(2020,11,1,12)
   >>> ed = datetime(2020,11,1,12)
   >>> region = 'mwam4' # default
   >>> mission = 's3a' # default
   >>> twin = 30 # default

   >>> # landmask filter
   >>> sco_lm = sc(sd,edate=ed,provider='eumetsat',land_mask=True,filterData=True,level='L2')

.. note::

   More examples with filters are coming soon ...

5. access/read model data
#########################
Model output can be accessed and read using the modelmod module. The modelmod config file model_specs.yaml needs adjustments if you want to include a model that is not present as default. Given that the model output file you would like to read in follows the cf-conventions and standard_names are unique, the minimum information you have to provide are usually:

.. code-block:: yaml

   modelname:
       path_template:
       file_template:
       init_times: []
       init_step:

Often there are ambiguities due to the multiple usage of standard_names. Any such problem can be solved here in the config-file by adding a variable like:

.. code-block:: yaml

    vars:
        Hs: VHM0
        time: time
        lons: lon
        lats: lat

.. code-block:: python3

   >>> from datetime import datetime
   >>> from wavy.modelmod import model_class as mc
   >>> model = 'mwam4' # default
   >>> varalias = 'Hs' # default
   >>> sd = datetime(2020,11,1)
   >>> ed = datetime(2020,11,2)
   >>> mco = mc(sdate=sd) # one time slice
   >>> mco = mc(sdate=sd,edate=ed) # time period
   >>> mco = mc(sdate=sd,leadtime=12) # time slice with lead time

The output will be something like::

   >>> mco = mc(sdate=sd)
   Time used for retrieving model data: 1.88 seconds
    ### model_class object initialized ###
   >>> mco.
   mco.edate             mco.get_item_parent(  mco.stdvarname
   mco.fc_date           mco.leadtime          mco.varalias
   mco.filestr           mco.model             mco.varname
   mco.get_item_child(   mco.sdate             mco.vars
   >>> mco.vars.keys()
   dict_keys(['longitude', 'latitude', 'time', 'datetime', 'time_unit', 'sea_surface_wave_significant_height', 'meta', 'leadtime'])

.. note::

   Even though it is possible to access a time period, **wavy** is not yet optimized to do so and the process will be slow. The reason being the ambiguous use of lead times. Whenever the keyword "leadtime" is None, a best guess is assumed and retrieved.

6. read in-situ observations (.d22 and netcdf/thredds)
######################################################

Currently two data types can be read .d22-files and netcdf-files.

read .d22 files
***************

.d22-files can be read in by adjusting d22_var_dicts config file. Currently, there are wave related variables included. Other variables like wind are about to be included. Another config-file that needs adjustment is the insitu_specs.yaml. There you need to define specs related to the in-situ observation of choice as well as path and filename. A call for the retrieval of an in-situ time series could be like:

.. code-block:: python3

   >>> from datetime import datetime
   >>> from wavy.insitumod import insitu_class as ic
   >>> varalias = 'Hs' # default
   >>> sd = datetime(2020,1,1,0)
   >>> ed = datetime(2020,1,5,0)
   >>> nID = 'ekofiskL'
   >>> sensor = 'waverider'
   >>> ico = ic(nID,sensor,sd,ed)

In contrast to the L3 satellite time series, in-situ time series are not filtered or underwent rigorous outlier detection. There are various operations that can be performed to massage the time series as you wish.It is in particular interesting to remove double reported values, which is often the case. This is done with setting unique=True.

.. code-block:: python3

   >>> ico = ic(nID,sensor,sd,ed,unique=True)

read .nc-files
**************

.. code-block:: python3

   >>> from datetime import datetime
   >>> from wavy.insitumod import insitu_class as ic
   >>> varalias = 'Hs' # default
   >>> sd = datetime(2020,1,1,0)
   >>> ed = datetime(2020,1,5,0)
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
   >>> stdname = ico.stdvarname
   >>> fig = plt.figure(figsize=(9,3.5))
   >>> ax = fig.add_subplot(111)
   >>> ax.plot(ico.vars['datetime'],ico.vars[stdname],'ko',label='raw')
   >>> ax.plot(ico_bm.vars['datetime'],ico_bm.vars[stdname],'r-',label='hourly blockMean')
   >>> plt.legend(loc='upper left')
   >>> plt.ylabel('Hs [m]')
   >>> plt.show()

.. image:: ./docs_fig_ts_insitu.png
   :scale: 80

.. note::

   It is important to note that due to different sampling frequencies there are still amibiguities that will have to be removed in future fixes.

7. collocating model and observations
#####################################
One of the main focus of wavy is to ease the collocation of observations and numerical wave models for the purpose of model validation. For this purpose there is the config-file collocation_specs.yaml where you can specify the name and path for the collocation file to be dumped.

Collocation of satellite and wave model
****************************************

.. code-block:: python3

   >>> from datetime import datetime
   >>> from wavy.satmod import satellite_class as sc
   >>> from wavy.collocmod import collocation_class as cc

   >>> model = 'mwam4' # default
   >>> mission = 's3a' # default
   >>> varalias = 'Hs' # default
   >>> sd = datetime(2020,11,1,12)
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

   >>> sd = datetime(2020,11,1)
   >>> ed = datetime(2020,11,2)
   >>> sco = sc(sdate=sd,edate=ed,region=model,mission=mission,varalias=varalias)
   >>> cco = cc(model=model,obs_obj_in=sco,distlim=6,date_incr=1)

Collocation of in-situ data and wave model
******************************************

The following example may take a few minutes.

.. code-block:: python3

   >>> # imports
   >>> from datetime import datetime
   >>> from wavy.insitumod import insitu_class as ic
   >>> from wavy.collocmod import collocation_class as cc

   >>> # settings
   >>> model = 'mwam4' # default
   >>> varalias = 'Hs' # default
   >>> sd = datetime(2020,1,1,1)
   >>> ed = datetime(2020,1,4,0)
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

   $ ./wavyQuick.py -sat s3a -reg mwam4 -sd 2020110100 -ed 2020110300 --show

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

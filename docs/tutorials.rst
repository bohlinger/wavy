Tutorials
==========
The following tutorials describe how to setup the wavy config files, to retrieve satellite and insitu data, to perform collocation and validation tasks. **wavy** is intended to be fleksibel such that customization can be achieved with minimal changes in the code. The **wavy** config files build the fundament for this approach and serve a similar purpose as namelist files often used for runtime changes in numerical modelling.

In general, executable files usually have help function which can be read using e.g.:

.. code-block:: bash

   $ ./{Executable}.py -h

e.g.:

.. code-block:: bash

   $ cd ~/wavy/apps/standalone
   $ ./wavyDownload.py -h

.. toctree::
   :maxdepth: 3
   :glob:

   tutorials_config
   tutorials_satmod
   tutorials_regions
   tutorials_filters

6. access/read model data
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
   >>> ed = "2020-11-2"
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

For the modelclass objects a quicklook fct exists to depict a certain time step of what you loaded::

   >>> mco.quicklook() # for a map


.. note::

   Even though it is possible to access a time period, **wavy** is not yet optimized to do so and the process will be slow. The reason, being the ambiguous use of lead times, will be improved in future versions.

7. read in-situ observations (.d22 and netcdf/thredds)
######################################################

Currently two data types can be read .d22-files and netcdf-files. Edit the insity_specs.yaml file in your config folder and adjust the directories.

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

Again, for the insitu class there is also a quicklook fct available::

   >>> ico.quicklook()

8. collocating model and observations
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
   >>> sco = sc(sdate=sd,region=model,mission=mission,varalias=varalias)
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

For the collocation class object there is also a quicklook fct implemented which allows to view both the time series and a map as for the satellite class object::

   >>> cco.quicklook(ts=True)
   >>> cco.quicklook(m=True)

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

9. dump collocation ts to a netcdf file
#######################################
The collocation results can now be dumped to a netcdf file. The path and filename can be entered as keywords but also predefined config settings can be used from collocation_specs.yaml:

.. code-block:: python3

   >>> cco_raw.write_to_nc()

10. validate the collocated time series
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

11. quick look examples
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

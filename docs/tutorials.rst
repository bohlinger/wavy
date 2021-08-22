Tutorials
==========
The following tutorials describe how to setup the wavy config files, to perform collocation and validation tasks. **wavy** is intended to be fleksibel such that customization can be achieved with minimal changes in the code. The **wavy** config files build the fundament for this approach and serve a similar purpose as namelist files often used for runtime changes in numerical modelling.

In general, executable files usually have help function which can be read using e.g.:

.. code-block:: bash

   $ ./{Executable}.py -h

e.g.:

.. code-block:: bash

   $ cd ~/wavy/apps/standalone
   $ ./download.py -h

1. **wavy** config files, a brief overview.
###########################################
**wavy** obtains custom information from config files. There are default versions which can be adjusted to user needs. The following config files exist as per now:


.. code-block:: bash

   $ ls
   buoy_specs.yaml.default         region_specs.yaml.default
   collocation_specs.yaml.default  satellite_specs.yaml.default
   d22_var_dicts.yaml.default      station_specs.yaml.default
   model_specs.yaml.default        validation_specs.yaml.default
   quicklook_specs.yaml.default    variable_info.yaml.default

Another config file containing customized information to adjust the automated plotting routines is in the making.

**wavy** browses the directory structure as follows:

    * check if env 'WAVY_CONFIG' is set or specified in .env
    * check if a config folder exists using xdg
    * fall back to default files within the package

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

   altimeter:
       cmems:
           local:
               path_template: /home/patrikb/tmp_altimeter/L3/mission

You can proceed now and download L3 data using the download.py script:

.. code-block:: bash

   $ cd ~/wavy/apps/standalone

To get help check ...

.. code-block:: bash

   $ ./download.py -h

... then download some satellite altimeter data:

.. code-block:: bash

   $ ./download.py -sat s3a -sd 2020110100 -ed 2020111000 -provider cmems

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

    altimeter:
        eumetsat:
            local:
                path_template: /home/patrikb/tmp_altimeter/L2/mission

... then download some satellite altimeter data:

.. code-block:: bash

   $ ./download.py -sat s3a -sd 2020110100 -ed 2020111000 -provider eumetsat

4. read satellite data
######################
Once the satellite data is downloaded one can access and read the data for further use in python. L3 data can be read like:

.. code-block:: python3

   >>> from datetime import datetime
   >>> from wavy.satmod import satellite_class
   >>> region = 'NorwegianSea'
   >>> varalias = 'Hs' # default
   >>> sat = 's3a' # default
   >>> sd = datetime(2020,11,1)
   >>> ed = datetime(2020,11,2)
   >>> sa_obj = satellite_class(sdate=sd,edate=ed,region=region)

This would results in a satellite_class object and the following output message::

   >>> sa_obj = satellite_class(sdate=sd,edate=ed,region=region)
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

   >>> sa_obj.
   sa_obj.edate               sa_obj.sat
   sa_obj.get_local_filelst(  sa_obj.sdate
   sa_obj.matchregion(        sa_obj.stdvarname
   sa_obj.matchregion_poly(   sa_obj.twin
   sa_obj.matchregion_rect(   sa_obj.varalias
   sa_obj.path_local          sa_obj.varname
   sa_obj.read_local_files(   sa_obj.vars
   sa_obj.region

With the retrieved variables in sa_obj.vars::

   >>> sa_obj.vars.keys()
   dict_keys(['time', 'latitude', 'longitude', 'sea_surface_wave_significant_height', 'time_unit', 'datetime', 'meta'])

5. access/read model data
#########################
Model output can be accessed and read using the modelmod module:

.. code-block:: python3

   >>> from datetime import datetime
   >>> from wavy.modelmod import model_class
   >>> model = 'mwam4' # default
   >>> varalias = 'Hs' # default
   >>> sd = datetime(2020,11,1)
   >>> ed = datetime(2020,11,2)
   >>> mc_obj = model_class(sdate=sd) # one time slice
   >>> mc_obj = model_class(sdate=sd,edate=ed) # time period
   >>> mc_obj = model_class(sdate=sd,leadtime=12) # time slice with lead time

The output will be something like::

   >>> mc_obj = model_class(sdate=sd)
   Time used for retrieving model data: 1.88 seconds
    ### model_class object initialized ###
   >>> mc_obj.
   mc_obj.edate       mc_obj.leadtime    mc_obj.stdvarname  mc_obj.vars
   mc_obj.fc_date     mc_obj.model       mc_obj.varalias
   mc_obj.filestr     mc_obj.sdate       mc_obj.varname
   >>> mc_obj.vars.keys()
   dict_keys(['longitude', 'latitude', 'time', 'datetime', 'time_unit', 'sea_surface_wave_significant_height', 'meta', 'leadtime'])

.. note::

   Even though it is possible to access a time period, **wavy** is not yet optimized to do so and the process will be slow. The reason being the ambiguous use of lead times. Whenever the keyword "leadtime" is None, a best guess is assumed and retrieved.

6. read in-situ observations (.d22)
###################################

7. collocating model and observations
#####################################

8. dump collocation ts to a netcdf file
#######################################

9. validate the collocated time series
######################################

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

   $ ./wavyQuick.py -sat s3a -reg mwam4 -sd 2020110100 -ed 2020110300 -dump /home/patrikb/tmp_altimeter/quickdump/

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
   Debiased Root Mean Squared Difference: 0.67
   Bias: 0.22
   Scatter Index: 12.71
   Mean of Model: 5.26
   Mean of Observations: 5.04
   Number of Collocated Values: 237

And of course the figure:

.. image:: ./docs_fig_sat_quicklook_005.png
   :scale: 40


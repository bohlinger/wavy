Import Satellite data
#####################

As a next step, you can access these data everywhere with wavy when having set your .env file and the WAVY_CONFIG path therein. For illustration purposes, let's create a directory called ~/project_tmp. This is now your project directory. Let's assume your ammended satellite_specs.yaml file is in this directory. In this directory create an .env file with the content e.g.:

.. code-block:: bash

   WAVY_CONFIG=/home/patrikb/project_tmp/

The content of your directory looks then like:

.. code-block:: bash

   (base) patrikb@pc5591:~/project_tmp$ ls -la
   total 20
   drwxrwxr-x  2 patrikb patrikb 4096 Aug  3 12:31 .
   drwx------ 52 patrikb patrikb 4096 Aug  3 13:51 ..
   -rw-rw-r--  1 patrikb patrikb   31 Aug  3 12:31 .env
   -rwxr-xr-x  1 patrikb patrikb 6257 Aug  3 12:31 satellite_specs.yaml

Now, open python in the wavy conda environment:

.. code-block:: bash
   
   $ conda activate wavy
   $ python

First you initialize a satellite class object with chosen parameters. 

.. code-block:: python3

   >>> # imports
   >>> from wavy.satellite_module import satellite_class as sc

   >>> # settings
   >>> region = 'global'
   >>> varalias = 'Hs' # default
   >>> name = 's3a'
   >>> nID = 'cmems_L3_NRT'
   >>> twin = 30 # default
   >>> sd = "2022-2-1 11" # can also be datetime object
   >>> ed = "2022-2-1 12" # not necessary if twin is specified

   >>> # retrieval
   >>> sco = sc(sd=sd,ed=ed,region=region,nID=nID,name=name)
   
Then you can import the corresponding data with the .populate() method. If you have not downloaded data you can use the repo-data for this exercise which is located in the wavy/test directory, e.g. 

.. code-block:: python3

    >>> tmpdir = '/home/patrikb/wavy/tests/data/L3/s3a'

.. code-block:: python3

    >>> sco = sco.populate(path=tmpdir)

Or in one line:

.. code-block:: python3

   >>> sco = sc(sd="2022-2-1 11",ed="2022-2-1 12", region="global",
                nID="cmems_L3_NRT",name="s3a").populate(path=tmpdir)

You have now read in 1 hour of significant wave height from the satellite mission s3a. The stdout message looks like::

   >>> sco = sc(sd="2022-2-1 11",ed="2022-2-1 12", region="global",nID="cmems_L3_NRT",name="s3a").populate(path=tmpdir)
   # ----- 
    ### Initializing satellite_class object ###
 
    Given kwargs:
   {'sd': '2022-2-1 11', 'ed': '2022-2-1 12', 'region': 'global', 'nID': 'cmems_L3_NRT', 'name': 's3a'}
 
    ### satellite_class object initialized ###
   # ----- 
    ### Read files and populate satellite_class object
    ## Find and list files ...
   23 valid files found
   source template: /home/patrikb/tmp_altimeter/L3/name/%Y/%m

   Checking variables..
    Get filevarname for 
   stdvarname: sea_surface_wave_significant_height 
   varalias: Hs
    !!! standard_name:  sea_surface_wave_significant_height  is not unique !!! 
   The following variables have the same standard_name:
    ['VAVH', 'VAVH_UNFILTERED']
    Searching *_cfg.yaml config file for definition
    Variable defined in *_cfg.yaml is:
   Hs = VAVH

   Choosing reader..
   #
   Environmental variable for WAVY_DIR not defined
   Defaults are chosen
   #
   Chosen reader: satellite_readers.read_local_ncfiles

   Reading..
   Reading 25 chunks of files with chunk size 1
   Total of 23 files
   100%|███████████████████████████████████████████| 24/24 [00:00<00:00, 43.18it/s]
    changing variables to aliases
    Get filevarname for 
   stdvarname: sea_surface_wave_significant_height 
   varalias: Hs
    !!! standard_name:  sea_surface_wave_significant_height  is not unique !!! 
   The following variables have the same standard_name:
    ['VAVH', 'VAVH_UNFILTERED']
    Searching *_cfg.yaml config file for definition
    Variable defined in *_cfg.yaml is:
   Hs = VAVH
      VAVH is alreade named correctly and therefore not adjusted
    Get filevarname for 
   stdvarname: time 
   varalias: time
    Get filevarname for 
   stdvarname: longitude 
   varalias: lons
      lons is alreade named correctly and therefore not adjusted
    Get filevarname for 
   stdvarname: latitude 
   varalias: lats
      lats is alreade named correctly and therefore not adjusted
    enforcing lon max min = -180/180
 
    ## Summary:
   4287 footprints retrieved.
   Time used for retrieving data:
   0.58 seconds
 
    ### satellite_class object populated ###
   # ----- 

The satellite_class object has multiple class methods and class variables:

.. code-block:: python3

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
  sco.filter_linearGAM(                         

With the retrieved variables in sco.vars::

   >>> sco.vars
   <xarray.Dataset> Size: 137kB
   Dimensions:  (time: 4287)
   Coordinates:
     * time     (time) datetime64[ns] 34kB 2022-02-01T10:30:00 ... 2022-02-01T12...
   Data variables:
       Hs       (time) float64 34kB 3.905 4.011 4.096 4.163 ... 2.15 2.18 2.203
       lons     (time) float64 34kB 165.0 165.0 165.0 165.0 ... -14.05 -14.13 -14.2
       lats     (time) float64 34kB 40.97 41.03 41.08 41.14 ... 69.16 69.11 69.06
   Attributes:
       title:    wavy dataset


You can readily explore what you obtained utilizing the quicklook function.

.. code-block:: python3

   >>> sco.quicklook(ts=True) # for time series
   >>> sco.quicklook(m=True) # for a map
   >>> sco.quicklook(a=True) # for all


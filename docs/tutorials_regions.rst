Define regions of interest
##########################

To retrieve satellite data for your region of interest you will need the region_specs.yaml file. Copy this file to your project directory and ammend it accordingly. Regions can be defined as given in the .default file as rectangular lat/lon regions, polygons, or by pointing to geojson files. 

.. note::

   When choosing to create a polynom make sure that it is closed by repeating the first coordinates in the end.

Continue with python:

.. code-block:: python3

   >>> # imports
   >>> from wavy.satmod import satellite_class as sc

   >>> # settings
   >>> sd = "2020-11-1 12"
   >>> ed = "2020-11-1 12"
   >>> region = 'NorwegianSea'
   >>> product = 'cmems_L3_NRT' # default
   >>> mission = 's3a' # default
   >>> varalias = 'Hs' # default
   
   >>> # retrieval
   >>> sco = sc(sdate=sd,edate=ed,region=region)

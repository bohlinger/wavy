Define regions of interest
##########################

To retrieve satellite data for your region of interest you will need the region_cfg.yaml file. Copy this file to your project directory and ammend it accordingly. Regions can be defined as given in the .default file as rectangular lat/lon regions, polygons, or by pointing to geojson files. 

.. note::

   When choosing to create a polynom make sure that it is closed by repeating the first coordinates in the end.

Continue with python, ...here is an example using a polygon for the Norwegian Sea:

- You can define the region directly when importing the data with wavy

.. code-block:: python3

   >>> # imports
   >>> from wavy.satellite_module import satellite_class as sc

   >>> # settings
   >>> sd = "2022-02-01"
   >>> ed = "2022-02-03"
   >>> region = 'NorwegianSea'
   >>> nID = 'cmems_L3_NRT' 
   >>> name = 's3a'
   >>> varalias = 'Hs' # default
   >>> path = '/path/to/your/wavy/tests/data/L3/s3a'
   
   >>> # retrieval
   >>> sco = sc(sd=sd, ed=ed, nID=nID, name=name, region=region).populate(path=path)

- Or you can crop it afterwards

.. code-block:: python3

   >>> sco = sc(sd=sd, ed=ed, nID=nID, name=name).populate(path=path)
   >>> sco = sco.crop_to_region('NorwegianSea')

It is also possible to define a rectangular custom region, directly in the python code, without having to specify it in the region_cfg.yaml file : 

.. code-block:: python3

   >>> sd = "2022-2-01 01"
   >>> ed = "2022-2-03 23"
   >>> name = 's3a'
   >>> varalias = 'Hs'
   >>> nID = 'cmems_L3_NRT'
   >>> region_dict = {'name': 'custom',
   ...                'region': {
   ...                 'llcrnrlon': -20.,
   ...                 'llcrnrlat': 40.,
   ...                 'urcrnrlon': 20.,
   ...                 'urcrnrlat': 80.}}

   >>> # init satellite_object
   >>> sco = sc(sd=sd,ed=ed,nID=nID,name=name,varalias=varalias)
   >>> # read data
   >>> sco = sco.populate(path=path, region=region_dict)

.. _consolidate-label:

Consolidate multiple retrievals
###############################

The consolidate_class enables the consolidation of multiple satellite_class or insitu_class objects. This is used under the hood in the mulitsat_class and multiins_class but can also be used manually as follows:

.. code-block:: python3

   >>> from wavy.satellite_module import satellite_class as sc
   >>> from wavy.consolidate import consolidate_class as cs

   >>> # satellite consolidate 
   >>> sco1 = sc(sd="2022-2-1",ed ="2022-2-3",region="NordicSeas", nID="cmems_L3_NRT",
   ...           name='s3a').populate(path='/path/to/your/wavy/tests/data/L3/s3a')
   >>> sco2 = sc(sd="2022-2-1",ed ="2022-2-3",region="NordicSeas", nID="cmems_L3_NRT",
   ...           name='s3b').populate(path='/path/to/your/wavy/tests/data/L3/s3b')

   >>> cso = cs([sco1,sco2])

The same can be done with insitu_class objects.

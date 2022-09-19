Consolidate multiple retrievals
###############################

The consolidate_class enables the consolidation of multiple satellite_class or insitu_class objects. This is used under the hood in the mulitsat_class and multiins_class but can also be used manually as follows:

.. code-block:: python3

   >>> # imports
   >>> from wavy.consolidate import consolidate_class as cs
   >>> from wavy.satmod import satellite_class as sc
   >>> # get satellite data rom s3a and s3b
   >>> sco1 = sc( sdate="2020-11-1 12",region="NordicSeas",mission='s3a' )
   >>> sco2 = sc( sdate="2020-11-1 13",region="NordicSeas",mission='s3b' )
   >>> # consolidate obs
   >>> cso = cs([sco1,sco2])

The same can be done with insitu_class objects.

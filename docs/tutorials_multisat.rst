Use multiple satellite missions
###############################

Sometimes it is most practical to retrieve data from multiple satellite missions for the same time period simultaneously. For this there is a multisat_class that consolidates satellite missions to one object but also keeps each satellite_class object for later use.

Open python in the wavy conda environment:

.. code-block:: bash
   
   $ conda activate wavy
   $ python

Import the library and proceed as with the satellite_class. The only difference is now that you can give the multisat_class a list of satellite_missions. Download your satellite data and then start the retrieval:

.. code-block:: python3

   >>> # imports
   >>> from wavy.satmod import multisat_class as ms

   >>> # retrieval
   >>> mso = ms( sdate   = "2020-11-1",
                 edate   = "2020-11-3",
                 region  = "global",
                 mission = ['s3a','s3b'] )

You have now retrieved 24h of significant wave height from 2 satellite missions.
The satellite_class object has multiple class methods and class variables:

.. code-block:: python3

   >>> # have a look
   >>> mso.quicklook(a=True,mode='indiv')

In contrast to the quicklook for the satellite_class object, there is an additional keyword added *model='indiv'*. This keyword is optional and affects the time series quicklook by plotting all missions in different colors. If the keyword is not used all data will be plotted in the same color.

The content of the multisat_class object should look like this:

.. code-block:: python3

    >>> mso.
    mso.edate       mso.obstype     mso.quicklook(  mso.units
    mso.label       mso.ocos        mso.region      mso.varalias
    mso.mission     mso.product     mso.sdate       mso.varname
    mso.obsname     mso.provider    mso.stdvarname  mso.vars

With the retrieved variables in mso.vars::

   >>> sco.vars.keys()
   dict_keys(['sea_surface_wave_significant_height', 'time', 'time_unit', 'latitude', 'longitude', 'datetime', 'meta'])

The class variable mso.ocos contains the individual satellite_class objects.

.. note 

    A class function for writing the data to pickle and netcdf will follow soon.

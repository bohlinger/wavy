Use multiple insitu-sources
###########################

Sometimes it is most practical to retrieve data from multiple insitu sources for the same time period simultaneously. For this there is a multiins_class that consolidates multiple insitu_class objects to one object but also keeps each insitu_class object for later use.

Open python in the wavy conda environment:

.. code-block:: bash
   
   $ conda activate wavy
   $ python

Import the library and proceed as with the insitu_class. The only difference is now that you can give the multiins_class a list of tags. Download your satellite data and then start the retrieval:

.. code-block:: python3

   >>> # imports
   >>> from wavy.multiins import multiins_class as mi

   >>> # retrieval
   >>> mio = mi( sdate = '2022-1-1',edate='2022-1-2',
                 tags=['E39','wave','pytest'],
                 varalias='Hs' )

You have now retrieved 24h of significant wave height from all insitu sources tagged with the strings 'E39', 'wave', and 'pytest'. This alos demonstrates that multiple tags are possible. The multiins_class object has multiple class methods and class variables:

.. code-block:: python3

   >>> # have a look
   >>> mio.quicklook(a=True,mode='indiv')

Just like with the multisat_class, there is an additional keyword added *model='indiv'*. This keyword is optional and affects the time series quicklook by plotting all retrievals in different colors. If the keyword is not used all data will be plotted in the same color.

The content of the multiins_class object should look like this:

.. code-block:: python3

   >>> mio.
   mio.edate       mio.obstype     mio.quicklook(  mio.varalias
   mio.label       mio.ocos        mio.sdate       mio.varname
   mio.mission     mio.product     mio.stdvarname  mio.vars
   mio.obsname     mio.provider    mio.units

With the retrieved variables in mio.vars::

   >>> mio.vars.keys()
   dict_keys(['sea_surface_wave_significant_height', 'longitude', 'latitude', 'time', 'time_unit', 'datetime'])

The class variable mso.ocos contains the individual insitu_class objects.

.. note 

    A class function for writing the data to pickle and netcdf will follow soon.

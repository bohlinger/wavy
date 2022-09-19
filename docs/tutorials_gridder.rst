Gridded satellite observations and statistics
#############################################

Once satellite observations are retrieved or even collocated model data are available **wavy** can display this data in custom grids for your region of interest.

Gridding of satellite observations
**********************************

Retrieve satellite observations from multiple satellites:

.. code-block:: python3

   >>> from wavy.satmod import satellite_class as sc
   >>> region = 'global'
   >>> sd = "2020-11-1"
   >>> ed = "2020-11-2"
   >>> sco = sc(sdate=sd,edate=ed,region='global')


Apply the gridder:

.. code-block:: python3

   >>> from wavy.gridder import gridder_class as gc
   >>> from wavy.grid_stats import apply_metric
   >>> bb = (-179,179,-80,80) # lonmin,lonmax,latmin,latmax
   >>> res = (5,5) # lon/lat
   >>> gco = gc(oco=sco,bb=bb,res=res)
   >>> var_gridded_dict,lon_grid,lat_grid = apply_metric(gco=gco)
   >>> gco.quicklook(val_grid=var_gridded_dict,lon_grid=lon_grid,lat_grid=lat_grid,metric='mor')

.. image:: ./docs_fig_gridder_obs.png
   :scale: 80

Information of the grid and the values from observations and model can also be obtained directly from the gridder_class object:

.. code-block:: python3

   >>> ovals,mvals,Midx = gco.get_obs_grid_idx()

ovals represent observation values, mvals are model values, and Midx is the matrix of indices. *mvals* is empty since no model values have been retrieved yet.

Gridding of collocated data
***************************
We first need to collocate the data with the collocation_class

.. code-block:: python3

   >>> from wavy.collocmod import collocation_class as cc
   >>> # collocate
   >>> cco = cc(model='mwam4',obs_obj_in=sco,distlim=6,date_incr=1)
   >>> # reduce region to part of model domain for better visual
   >>> bb = (-20,20,50,80) # lonmin,lonmax,latmin,latmax
   >>> res = (5,5) # lon/lat
   >>> gco = gc(cco=cco,bb=bb,res=res)
   >>> var_gridded_dict,lon_grid,lat_grid = apply_metric(gco=gco)
   >>> # plot all validation metrics on grid
   >>> gco.quicklook(val_grid=var_gridded_dict,lon_grid=lon_grid,lat_grid=lat_grid,metric='all')

.. |ex1| image:: ./docs_fig_gridder_coll_nov.png
   :scale: 50
.. |ex2| image:: ./docs_fig_gridder_coll_rmse.png
   :scale: 50

+-------------------+------------------+
| |ex1|             | |ex2|            |
|                   |                  |
+-------------------+------------------+

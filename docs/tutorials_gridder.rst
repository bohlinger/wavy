Gridded satellite observations and statistics
#############################################

Once satellite observations are retrieved or even collocated model data are available **wavy** can display this data in custom grids for your region of interest.

Gridding of satellite observations
**********************************

Retrieve satellite observations from multiple satellites:

.. code-block:: python3

   >>> from wavy.satmod import satellite_class as sc
   >>> region = 'global'
   >>> sd = "2021-12-1"
   >>> ed = "2022-3-1"
   >>> missions = ['s3a','s3b','c2','al','h2b','cfo']
   >>> sco = sc(sdate=sd,edate=ed,region='global',mission=missions)


Apply the gridder:

.. code-block:: python3

   >>> from wavy.gridder import gridder_class as gc
   >>> from wavy.grid_stats import apply_metric
   >>> bb = (-179,179,-80,80) # lonmin,lonmax,latmin,latmax
   >>> res = (5,5) # lon/lat
   >>> gco = gc(oco=sco,bb=bb,res=res)
   >>> var_gridded_dict,lon_grid,lat_grid = apply_metric(gco=gco)
   >>> gco.quicklook(val_grid=var_gridded_dict,lon_grid=lon_grid,lat_grid=lat_grid,metric='mor')

.. image:: ./docs_fig_ts_sat.png
   :scale: 80

Gridding of collocated data
***************************

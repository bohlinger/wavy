Read model data
###############

Model output can be accessed and read using the modelmod module. The modelmod config file model_specs.yaml needs adjustments if you want to include a model that is not present as default. Given that the model output file you would like to read follows the cf-conventions and standard_names are unique, the minimum information you have to provide are usually:

.. code-block:: yaml

   modelname:
       path_template:
       file_template:
       init_times: []
       init_step:

Often there are ambiguities due to the multiple usage of standard_names. Any such problem can be solved here in the config-file by adding the specified variable name like:

.. code-block:: yaml

    vardef:
        Hs: VHM0
        time: time
        lons: lon
        lats: lat

The variable aliases (left hand side) need to be specified in the variable_info.yaml. Basic variables are already defined. All specs listed here are also used when **wavy** writes the retrieved values to netcdf.

.. code-block:: python3

   >>> from wavy.modelmod import model_class as mc
   >>> model = 'mwam4' # default
   >>> varalias = 'Hs' # default
   >>> sd = "2020-11-1"
   >>> ed = "2020-11-2"
   >>> mco = mc(sdate=sd) # one time slice
   >>> mco_p = mc(sdate=sd,edate=ed) # time period
   >>> mco_lt = mc(sdate=sd,leadtime=12) # time slice with lead time

Whenever the keyword "leadtime" is None, a best estimate is assumed and retrieved. The output will be something like::

   >>> mco = mc(sdate=sd)

   >>> mco.
   mco.edate             mco.leadtime          mco.units
   mco.fc_date           mco.model             mco.varalias
   mco.filestr           mco.quicklook(        mco.varname
   mco.get_item_child(   mco.sdate             mco.vars
   mco.get_item_parent(  mco.stdvarname

   >>> mco.vars.keys()
   dict_keys(['longitude', 'latitude', 'time', 'datetime', 'time_unit', 'sea_surface_wave_significant_height', 'meta', 'leadtime'])

For the modelclass objects a quicklook fct exists to depict a certain time step of what you loaded::

   >>> mco.quicklook() # for a map


.. note::

   Even though it is possible to access a time period, **wavy** is not yet optimized to do so and the process will be slow. The reason, being the ambiguous use of lead times, will be improved in future versions.


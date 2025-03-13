Retrieve AIS
############

**wavy** offers the possibility to retrieve AIS data which essentially is data on characteristics and location of marine vessels. This works as follows:

.. code-block:: python3

    >>> import wavy.ais_module as ais

    >>> bbox = ['5.89', '62.3', '6.5', '62.7']
    >>> sd = '2017-01-03 08'
    >>> ed = '2017-01-03 09'

    >>> ais_ds = ais.get_AIS_data(bbox, sd, ed)

The output is an xarray dataset and looks like the following:
   
.. code-block:: python3

    >>> ais_ds
    <xarray.Dataset> Size: 12MB
    Dimensions:       (time: 155292)
    Coordinates:
      * time          (time) datetime64[ns] 1MB 2017-01-02T09:00:00 ... 2017-01-0...
    Data variables:
        mmsi          (time) int64 1MB 259098000 219770000 ... 259097000 259257000
        lons          (time) float64 1MB 6.33 6.189 6.28 6.192 ... 6.279 6.337 6.124
        lats          (time) float64 1MB 62.38 62.56 62.6 ... 62.62 62.41 62.47
        grosstonnage  (time) float64 1MB 2.967e+03 218.0 769.0 ... 2.967e+03 380.0
        dwt           (time) float64 1MB 970.0 nan 250.0 27.0 ... nan nan 891.0 nan
        length        (time) float64 1MB 109.5 26.2 64.32 29.2 ... nan 109.5 34.6
        breadth       (time) float64 1MB 16.99 9.2 11.26 9.0 ... 10.6 nan 17.0 10.6
        draught       (time) float64 1MB 3.0 1.64 nan 1.5 1.64 ... nan nan 3.425 nan
        shiptypeeng   (time) object 1MB 'Vehicles Carrier' ... 'Passenger Ship'

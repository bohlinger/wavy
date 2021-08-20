Tutorials
==========
The following tutorials describe how to setup the wavy config files, to perform collocation and validation tasks. **wavy** is intended to be fleksibel such that customization can be achieved with minimal changes in the code. The **wavy** config files build the fundament for this approach and serve a similar purpose as namelist files often used for runtime changes in numerical modelling.

#. The **wavy** config files, a brief overview.

#. download L3 satellite altimetry data
   L3 satellite data is obtained from Copernicus with the product identifier WAVE_GLO_WAV_L3_SWH_NRT_OBSERVATIONS_014_001. User credentials are required for this task. So before you can start you have to get a Copernicus account (free of costs).

#. download L2 stallite altimetry data
   L2 satellite data are obtained from eumetsat using colhub and the SentinelAPI. This requires user credentials for eumetsat and colhub, which are free of costs as well.

#. read satellite data

#. access/read model data

#. read in-situ observations

#. collocating model and observations

#. dump collocation ts to a netcdf file

#. validate the collocated time series

#. quick look examples

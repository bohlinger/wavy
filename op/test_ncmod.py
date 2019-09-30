import os
import netCDF4
import numpy as np

"""
Additional information on variable retrieved from wave model
"""
var_dict = {
                'time':
                    {
                    'units':'seconds since',
                    'standard_name':'time',
                    },
                'Hs':
                    {
                    'units':'m',
                    'standard_name':'significant_height_from_wave_model',
                    'valid_range':[0., 30.]
                    },
                'Tp':
                    {
                    'units':'s',
                    'standard_name':'peak_period_from_wave_model',
                    'valid_range':[0., 30.]
                    },
                'Tm':
                    {
                    'units':'s',
                    'standard_name':'mean_period_from_wave_model',
                    'valid_range':[0., 30.]
                    },
                }


def dumptonc_ts_pos(outpath,filename,title,basetime,\
                    coll_dict,model,varname):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file
    """
    time = coll_dict['time']
    var_model = coll_dict[varname]
    lons_model = coll_dict['lons_model']
    lats_model = coll_dict['lats_model']
    lons_pos = coll_dict['lons_pos']
    lats_pos = coll_dict['lats_pos']
    idx = coll_dict['idx']
    idy = coll_dict['idy']
    hdist = coll_dict['hdist']
    fullpath = outpath + filename
    print ('Dump data to file: ' + fullpath)
    if os.path.isfile(fullpath):
        nc = netCDF4.Dataset(
                        fullpath,mode='a',
                        clobber=False
                        )
        # variables
        startidx = len(nc['time'])
        endidx = len(nc['time'])+len(time)
        nc.variables['time'][startidx:endidx] = time[:]
        nc.variables[varname][startidx:endidx] = var_model[:]
        # nc.variables['longitude_model'][startidx:endidx] = lons_model[:]
        # nc.variables['latitude_model'][startidx:endidx] = lats_model[:]
    else:
        os.system('mkdir -p ' + outpath)
        nc = netCDF4.Dataset(
                        fullpath,mode='w',
                        format='NETCDF4'
                        )
        # global attributes
        nc.title = title
        nc.netcdf_version = "4"
        nc.processing_level = "No post-processing performed"
        nc.static_position_station =  ("Latitude: "
                            + "{:.4f}".format(lats_pos[0])
                            + ", Longitude: "
                            + "{:.4f}".format(lons_pos[0]))
        nc.static_position_model =  ("Latitude: "
                            + "{:.4f}".format(lats_model[0])
                            + ", Longitude: "
                            + "{:.4f}".format(lons_model[0]))
        nc.static_collocation_idx =  ("idx: "
                            + str(idx[0])
                            + ", idy: "
                            + str(idy[0]))
        # dimensions
        dimsize = None
        dimtime = nc.createDimension(
                                'time',
                                size=dimsize
                                )
        # variables
        nctime = nc.createVariable(
                               'time',
                               np.float64,
                               dimensions=('time')
                               )
        ncvar_model = nc.createVariable(
                               varname,
                               np.float64,
                               dimensions=('time')
                               )
        # generate time for netcdf file
        # time
        nctime.standard_name = var_dict['time']['standard_name']
        nctime.units = var_dict['time']['units'] + ' ' + str(basetime)
        nctime[:] = time
        # var_model
        ncvar_model.standard_name = var_dict[varname]['standard_name']
        ncvar_model.units = var_dict[varname]['units']
        ncvar_model.valid_range = var_dict[varname]['valid_range'][0], \
                                  var_dict[varname]['valid_range'][1]
        ncvar_model[:] = var_model
    nc.close()

def test_dump ():
  import numpy as np
  from numpy import array
  import datetime

  out = 'out/'
  fn = 'test.nc'
  title = 'asdf'
  bt = datetime.datetime (1970, 1, 1)
  coll_dict = {'lats_pos': [56.54], 'lats_model': [56.55768265623981], 'basetime': datetime.datetime(1970, 1, 1, 0, 0), 'idx': [array([273])], 'time': [1546398000.0], 'Hs': [5.6404104], 'idy': [array([417])], 'lons_pos': [3.21], 'lons_model': [3.2104459637242972], 'hdist': [1.9651770580555135]}

  model = 'mwam4'
  varname = 'Hs'

  for i in range (1000):
    coll_dict['time'][0] += 1
    coll_dict['Hs'][0] += 2
    # print (coll_dict)
    dumptonc_ts_pos (out, fn, title, bt, coll_dict, model, varname)

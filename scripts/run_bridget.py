'''
 Run the Bridge Model for each gridpoint in the given netcdf file, yeah
'''
import netCDF4
import sys
import pytz
import datetime
import numpy as np

(USERWIS, USEMODEL) = range(2)

def find_initts( nc ):
    ''' Provided the given netcdf file object, figure out the start time '''
    tm = nc.variables['time']
    ts = datetime.datetime.strptime( tm.units[14:], "%Y-%m-%d %H:%M:%S")
    ts = ts.replace(tzinfo=pytz.timezone("UTC"))
    return ts

def make_output(nc, initts):
    ''' Generate an output file to hold our results '''
    fn = '%s_output.nc' % (initts.strftime("%Y%m%d%H%M"),)
    ncout = netCDF4.Dataset(fn, 'w')
    # Setup dimensions
    ncout.createDimension('i_cross', len(nc.dimensions['i_cross']))
    ncout.createDimension('j_cross', len(nc.dimensions['j_cross']))
    ncout.createDimension('time', 72)

    # Setup variables
    tm = ncout.createVariable('time', np.int, ('time',))
    tm.units = "minutes since %s" % (initts.strftime("%Y-%m-%d %H:%M:%S"),)
    
    i_cross = ncout.createVariable('i_cross', np.float, ('i_cross',))
    i_cross.units = "m"
    i_cross[:] = range(len(nc.dimensions['i_cross']))

    j_cross = ncout.createVariable('j_cross', np.float, ('j_cross',))
    j_cross.units = "m"
    j_cross[:] = range(len(nc.dimensions['j_cross']))

    lat = ncout.createVariable('lat', np.float, ('i_cross', 'j_cross'))
    lat.long_name = "latitude"
    lat.units = "degrees_north"
    lat[:] = nc.variables['latitcrs'][:]

    lon = ncout.createVariable('lon', np.float, ('i_cross', 'j_cross'))
    lon.long_name = "longitude"
    lon.units = "degrees_east"
    lon[:] = nc.variables['longicrs'][:]

    bdeckt = ncout.createVariable('bdeckt', np.float, 
                                  ('time', 'i_cross', 'j_cross'))
    bdeckt.coordinates = "lon lat"
    bdeckt.units = "K"
    bdeckt.long_name = 'Bridge Deck Temperature'

    h = ncout.createVariable('h', np.float, 
                                  ('time', 'i_cross', 'j_cross'))
    h.coordinates = "lon lat"
    #h.units = "m"
    #h.long_name = 'Depth of Frost'

    swout = ncout.createVariable('swout', np.float, 
                                  ('time', 'i_cross', 'j_cross'))
    swout.coordinates = "lon lat"
    swout.units = "W m s-2"
    swout.long_name = 'Shortwave outgoing'

    lwout = ncout.createVariable('lwout', np.float, 
                                  ('time', 'i_cross', 'j_cross'))
    lwout.coordinates = "lon lat"
    lwout.units = "W m s-2"
    lwout.long_name = 'Longwave outgoing'

    lf = ncout.createVariable('lf', np.float, 
                               ('time', 'i_cross', 'j_cross'))
    lf.coordinates = "lon lat"

    tmpk = ncout.createVariable('tmpk', np.float, 
                                 ('time', 'i_cross', 'j_cross'))
    tmpk.coordinates = "lon lat"
    tmpk.units = 'K'

    dwpk = ncout.createVariable('dwpk', np.float, 
                                 ('time', 'i_cross', 'j_cross'))
    dwpk.coordinates = "lon lat"

    wmps = ncout.createVariable('wmps', np.float, 
                                ('time', 'i_cross', 'j_cross'))
    wmps.coordinates = "lon lat"
  
    ifrost = ncout.createVariable('ifrost', np.int, 
                                  ('time', 'i_cross', 'j_cross'))
    ifrost.coordinates = "lon lat"
    ifrost.missing_value = 0

    frostd = ncout.createVariable('frostd', np.float, 
                                  ('time', 'i_cross', 'j_cross'))
    frostd.coordinates = "lon lat"
    frostd.missing_value = -99.


    ncout.close()
    return netCDF4.Dataset(fn)

def run_model(nc, initts, ncout):
    ''' Actually run the model, please '''
    t2 = nc.variables['t2']
    u10 = nc.variables['u10']
    v10 = nc.variables['v10']
    tm = nc.variables['time']
    lwdown = nc.variables['lwdown']
    swdown = nc.variables['swdown']
    q2 = nc.variables['q2']
    rc = nc.variables['rain_con']
    rn = nc.variables['rain_non']
    lats = nc.variables['latitcrs']
    lons = nc.variables['longicrs']

    for i in range(len(nc.dimensions['i_cross'])):
        for j in range(len(nc.dimensions['j_cross'])):
            pass

if __name__ == '__main__':
    ''' Do something please '''
    fn = sys.argv[1]
    nc = netCDF4.Dataset(fn)
    
    initts = find_initts( nc )
    ncout = make_output(nc, initts)
    run_model(nc, initts, ncout)
    
    nc.close()
    ncout.close()
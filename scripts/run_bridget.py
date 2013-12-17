'''
 Run the Bridge Model for each gridpoint in the given netcdf file, yeah
'''
import netCDF4
import sys
import pytz
import datetime
import numpy as np
import subprocess
from pyiem.datatypes import temperature
import os

IOFFSET = 62
JOFFSET = 70

def find_initts( nc ):
    ''' Provided the given netcdf file object, figure out the start time '''
    tm = nc.variables['time']
    ts = datetime.datetime.strptime( tm.units[14:], "%Y-%m-%d %H:%M:%S")
    ts = ts.replace(tzinfo=pytz.timezone("UTC"))
    return ts

def make_output(nc, initts):
    ''' Generate an output file to hold our results '''
    fn = 'output/%s_output.nc' % (initts.strftime("%Y%m%d%H%M"),)
    ncout = netCDF4.Dataset(fn, 'w')
    # Setup dimensions
    ncout.createDimension('i_cross', len(nc.dimensions['i_cross']))
    ncout.createDimension('j_cross', len(nc.dimensions['j_cross']))
    ncout.createDimension('time', 72*60)

    # Setup variables
    tm = ncout.createVariable('time', np.int, ('time',))
    tm.units = "minutes since %s" % (initts.strftime("%Y-%m-%d %H:%M:%S"),)
    tm[:] = range(72*60)
    
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
    return netCDF4.Dataset(fn, 'a')

def make_rwis(i, j, initts, oldncout):
    ''' Generate spinup file '''
    if oldncout is None:
        if not os.path.isfile('faux_rwis.txt'): # create one!
            o = open('faux_rwis.txt', 'w')
            for hr in range(-12,0,1):
                o.write("%s     44      44     10\n" % (
                 (initts + datetime.timedelta(hours=hr)).strftime("%Y%m%d%H%M"),))
            o.close()
        return 'faux_rwis.txt'

    i = i - IOFFSET
    j = j - JOFFSET
    # Generate the rwis.txt file
    ts0 = find_initts(oldncout)
    o = open('rwis.txt', 'w')
    for tstep in range(1, len(oldncout.dimensions['time']), 60):
        ts = ts0 + datetime.timedelta(
                                minutes=int(oldncout.variables['time'][tstep]))
        if ts >= initts:
            break
        o.write("%s %7.2f %7.2f %7.2f\n" % ( ts.strftime("%Y%m%d%H%M"), 
            temperature(oldncout.variables['tmpk'][tstep,i,j], 'K').value("F"), 
            temperature(oldncout.variables['bdeckt'][tstep,i,j], "K").value("F"), 
            (oldncout.variables['wmps'][tstep,i,j])*2.0 ) )
    
    o.close()
    return 'rwis.txt'

def run_model(nc, initts, ncout, oldncout):
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

    otmpk = ncout.variables['tmpk']

    #mini = 200
    #minj = 200
    #maxi = 0
    #maxj = 0
    for i in range(len(nc.dimensions['i_cross'])):
        for j in range(len(nc.dimensions['j_cross'])):
            lat = lats[i,j]
            lon = lons[i,j]
            '''Hey, we only care about Iowa data! -97 40 -90 43.75'''
            if lat < 40 or lat > 43.75 or lon < -97 or lon > -90:
                continue
            rwisfn = make_rwis(i, j, initts, oldncout)
            #mini = min(i, mini)
            #minj = min(j, minj)
            #maxi = max(i, maxi)
            #maxj = max(j, maxj)
            #continue
            modelfp = open('modeldata.txt', 'w')
            for t in range(1, len(nc.dimensions['time'])):
                ts = initts + datetime.timedelta(minutes=int(tm[t]))
                
                modelfp.write("%s %6.1f %6.2f %7.6f %7.2f %7.2f %7.2f %7.4f\n"%( 
                              ts.strftime("%Y-%m-%d_%H:%M:%S00"), tm[t],
                              t2[t,i,j], q2[t,i,j], 
                              (u10[t,i,j]**2 + v10[t,i,j]**2)**0.5,
                              swdown[t,i,j], lwdown[t,i,j],
                              rc[t,i,j] + rn[t,i,j]
                              ))
                
            modelfp.close()
            
            proc = subprocess.Popen("./bridgemodel %s %s %s" % (rwisfn,
                                                    "modeldata.txt", 
                                                    "propertyFile"),
                                    shell=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            sd = proc.stdout.read()

            # Process the output file! 
            for line in open('pavetemp.out'):
                tokens = line.split()
                if len(tokens) < 7:
                    continue
                ts = datetime.datetime.strptime(tokens[0], "%Y-%m-%d_%H:%M")
                ts = ts.replace(tzinfo=pytz.timezone("UTC"))
                t = int((ts - initts).days * 1400 + ((ts - initts).seconds / 60))
                otmpk[t,i,j] = float(tokens[1])

    # ncks -d i_cross,62,82 -d j_cross,70,98 201312131200_output.nc 201312131200_output2.nc
    #print mini, minj, maxi, maxj #62 70 82 98

def find_last_output(initts):
    ''' See if we have a previous run on file, that can be used to spin up
    our current run '''
    for i in range(-12,-73,-12):
        ts = initts + datetime.timedelta(hours=i)
        testfn = 'output/%s_iaoutput.nc' % (ts.strftime("%Y%m%d%H%M"),)
        if os.path.isfile(testfn):
            print '  Using %s as warmup values' % (testfn,)
            return netCDF4.Dataset(testfn, 'r')
    print 'Did not find a previous output, will use dummy RWIS data :('
    return None

def downsize_output(initts):
    ''' Subset the output file, so to save some space 66% actually '''
    fn1 = "output/%s_output.nc" % (initts.strftime("%Y%m%d%H%M"),)
    fn2 = "output/%s_iaoutput.nc" % (initts.strftime("%Y%m%d%H%M"),)
    if os.path.isfile(fn2):
        os.unlink(fn2)
    cmd = "ncks -d i_cross,%s,82 -d j_cross,%s,98 %s %s" % (IOFFSET, JOFFSET,
                                                            fn1, fn2)
    p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE,
                         stdout=subprocess.PIPE)
    data= p.stdout.read()
    os.unlink(fn1)

if __name__ == '__main__':
    ''' Do something please '''
    fn = sys.argv[1]
    nc = netCDF4.Dataset(fn)
    
    initts = find_initts( nc )
    ncout = make_output(nc, initts)
    oldncout = find_last_output(initts)
    run_model(nc, initts, ncout, oldncout)
    downsize_output( initts )
    if os.path.isfile('faux_rwis.txt'):
        os.unlink('faux_rwis.txt')
    
    nc.close()
    ncout.close()
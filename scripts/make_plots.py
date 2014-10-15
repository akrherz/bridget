'''
 Generate plots from the Frost Model output!
'''
import netCDF4
import sys
import datetime
import pytz
from pyiem.plot import MapPlot
from pyiem.datatypes import temperature
import numpy as np

PVARS = {
 'bdeckt': {'title': 'Bridge Deck Temperature [F]', 'units': 'F',
            'levels': np.arange(-30,61,5)},
 'frostd': {'title': 'Frost Depth [mm]', 'units': 'mm',
            'levels': np.arange(0.0001,0.05,0.005)},
 'ifrost': {'title': 'Frost Present (yes/no)', 'units': '',
            'levels': [0,1,2]}         
}

def compute_sts( nc ):
    ''' Figure out the start time of this netcdf file '''
    tm = nc.variables['time'].units.replace("minutes since ", "")
    dt = datetime.datetime.strptime(tm, "%Y-%m-%d %H:%M:%S")
    return dt.replace(tzinfo=pytz.timezone("UTC"))

def make_plots(nc):
    ''' Generate some plots '''
    sts = compute_sts(nc)
    lats = nc.variables['lat'][:]
    lons = nc.variables['lon'][:]
    rts = (sts.astimezone(pytz.timezone("America/Chicago"))).strftime(
                                                            "%d %b %Y %H %p")
    for i, tm in enumerate(nc.variables['time'][:]):
        dt = sts + datetime.timedelta(minutes=float(tm))
        if dt.minute != 0:
            continue
        fhour = int( tm / 60.0 )
        fts = (dt.astimezone(pytz.timezone("America/Chicago"))).strftime(
                                                            "%d %b %Y %H %p")
        for pvar in PVARS:
            m = MapPlot(title='ISUMM5/Bridget Modelled %s' % (
                                                    PVARS[pvar]['title'],),
                        subtitle='Model Run: %s Forecast Valid: %s' % (rts, fts))
            vals = nc.variables[pvar][i,:,:]
            if pvar == 'bdeckt':
                vals = temperature(vals, 'K').value('F')
            m.pcolormesh(lons, lats, vals, PVARS[pvar]['levels'], units='mm')
            pqstr = "plot c %s model/frost/bridget/%02i/%s_%02i_f%03i.png bogus png" % (
                                        sts.strftime("%Y%m%d%H%M"), sts.hour,
                                        pvar, sts.hour, fhour)
            m.postprocess(pqstr=pqstr)
            m.close()
                      
if __name__ == '__main__':
    ''' See how we are called '''
    fn = sys.argv[1]
    nc = netCDF4.Dataset(fn)
    make_plots(nc)
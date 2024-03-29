"""Run the Bridge Model for each gridpoint in the given netcdf file, yeah"""
import os
import shutil
import sys
import datetime
import subprocess

import pytz
import numpy as np
import numpy.ma as ma
import netCDF4
from pyiem.datatypes import temperature
from pyiem.util import logger

LOG = logger()
IOFFSET = 62
JOFFSET = 70
CONDITIONS = ["Dry", "frosty", "Icy/Snowy", "Melting", "Freezing", "Wet"]


def find_initts(nc):
    """Provided the given netcdf file object, figure out the start time"""
    tm = nc.variables["time"]
    ts = datetime.datetime.strptime(tm.units[14:], "%Y-%m-%d %H:%M:%S")
    ts = ts.replace(tzinfo=pytz.UTC)
    return ts


def make_output(nc, initts):
    """Generate an output file to hold our results"""
    fn = "output/%s_output.nc" % (initts.strftime("%Y%m%d%H%M"),)
    LOG.info("creating output: %s", fn)
    ncout = netCDF4.Dataset(fn, "w")
    # Setup dimensions
    ncout.createDimension("i_cross", len(nc.dimensions["i_cross"]))
    ncout.createDimension("j_cross", len(nc.dimensions["j_cross"]))
    ncout.createDimension("time", 72 * 4 + 1)

    # Setup variables
    tm = ncout.createVariable("time", int, ("time",))
    tm.units = "minutes since %s" % (initts.strftime("%Y-%m-%d %H:%M:%S"),)
    tm[:] = range(0, 72 * 60 + 1, 15)

    i_cross = ncout.createVariable("i_cross", float, ("i_cross",))
    i_cross.units = "m"
    i_cross[:] = range(len(nc.dimensions["i_cross"]))

    j_cross = ncout.createVariable("j_cross", float, ("j_cross",))
    j_cross.units = "m"
    j_cross[:] = range(len(nc.dimensions["j_cross"]))

    lat = ncout.createVariable("lat", float, ("i_cross", "j_cross"))
    lat.long_name = "latitude"
    lat.units = "degrees_north"
    lat[:] = nc.variables["latitcrs"][:]

    lon = ncout.createVariable("lon", float, ("i_cross", "j_cross"))
    lon.long_name = "longitude"
    lon.units = "degrees_east"
    lon[:] = nc.variables["longicrs"][:]

    # DEFINE THE VARS, PLEASE!
    dims = ("time", "i_cross", "j_cross")
    icond = ncout.createVariable("icond", np.byte, dims)
    icond.coordinates = "lon lat"
    icond.long_name = "Pavement Textual Condition"
    icond.value0 = "Dry"
    icond.value1 = "frosty"
    icond.value2 = "Icy/Snowy"
    icond.value3 = "Melting"
    icond.value4 = "Freezing"
    icond.value5 = "Wet"

    bdeckt = ncout.createVariable("bdeckt", float, dims)
    bdeckt.coordinates = "lon lat"
    bdeckt.units = "K"
    bdeckt.long_name = "Bridge Deck Temperature"
    bdeckt.missing_value = np.array(1e20, bdeckt.dtype)

    h = ncout.createVariable("h", float, dims)
    h.coordinates = "lon lat"
    h.missing_value = np.array(1e20, h.dtype)

    swout = ncout.createVariable("swout", float, dims)
    swout.coordinates = "lon lat"
    swout.units = "W m s-2"
    swout.long_name = "Shortwave outgoing"
    swout.missing_value = np.array(1e20, swout.dtype)

    lwout = ncout.createVariable("lwout", float, dims)
    lwout.coordinates = "lon lat"
    lwout.units = "W m s-2"
    lwout.long_name = "Longwave outgoing"
    lwout.missing_value = np.array(1e20, lwout.dtype)

    lf = ncout.createVariable("lf", float, dims)
    lf.coordinates = "lon lat"
    lf.missing_value = np.array(1e20, lf.dtype)

    tmpk = ncout.createVariable("tmpk", float, dims)
    tmpk.coordinates = "lon lat"
    tmpk.units = "K"
    tmpk.missing_value = np.array(1e20, tmpk.dtype)

    dwpk = ncout.createVariable("dwpk", float, dims)
    dwpk.coordinates = "lon lat"
    dwpk.missing_value = np.array(1e20, dwpk.dtype)

    wmps = ncout.createVariable("wmps", float, dims)
    wmps.coordinates = "lon lat"
    wmps.missing_value = np.array(1e20, wmps.dtype)

    ifrost = ncout.createVariable("ifrost", int, dims)
    ifrost.coordinates = "lon lat"
    ifrost.missing_value = 0
    ifrost.missing_value = -1

    frostd = ncout.createVariable("frostd", float, dims)
    frostd.coordinates = "lon lat"
    frostd.missing_value = -99.0
    frostd.missing_value = np.array(1e20, frostd.dtype)

    ncout.close()
    return netCDF4.Dataset(fn, "a")


def make_rwis(i, j, initts, oldncout, modeltemp):
    """Generate spinup file"""
    if oldncout is None:
        with open("faux_rwis.txt", "w", encoding="utf8") as fh:
            for hr in range(-12, 0, 1):
                fh.write(
                    "%s     %.3f      %.3f     10\n"
                    % (
                        (initts + datetime.timedelta(hours=hr)).strftime(
                            "%Y%m%d%H%M"
                        ),
                        temperature(modeltemp, "K").value("F") + 5,
                        temperature(modeltemp, "K").value("F") + 5,
                    )
                )
        return "faux_rwis.txt"

    i = i - IOFFSET
    j = j - JOFFSET
    # Generate the rwis.txt file
    ts0 = find_initts(oldncout)
    obs = 0
    o = open("rwis.txt", "w", encoding="utf8")
    for tstep in range(0, len(oldncout.dimensions["time"]), 4):
        ts = ts0 + datetime.timedelta(
            minutes=int(oldncout.variables["time"][tstep])
        )
        if ts >= initts:
            break
        tmpf = temperature(oldncout.variables["tmpk"][tstep, i, j], "K").value(
            "F"
        )
        if ma.is_masked(tmpf) or tmpf < -50 or tmpf > 150:
            continue
        o.write(
            ("%s %7.2f %7.2f %7.2f\n")
            % (
                ts.strftime("%Y%m%d%H%M"),
                tmpf,
                temperature(
                    oldncout.variables["bdeckt"][tstep, i, j], "K"
                ).value("F"),
                (oldncout.variables["wmps"][tstep, i, j]) * 2.0,
            )
        )
        obs += 1
    o.close()
    if obs == 0:
        LOG.info("  i=%s j=%s found no temperature data in old file", i, j)
        with open("faux_rwis.txt", "w", encoding="utf8") as fh:
            for hr in range(-12, 0, 1):
                fh.write(
                    "%s     %.3f      %.3f     10\n"
                    % (
                        (initts + datetime.timedelta(hours=hr)).strftime(
                            "%Y%m%d%H%M"
                        ),
                        temperature(modeltemp, "K").value("F") + 5,
                        temperature(modeltemp, "K").value("F") + 5,
                    )
                )
        return "faux_rwis.txt"
    return "rwis.txt"


def run_model(nc, initts, ncout, oldncout):
    """Actually run the model, please"""
    t2 = nc.variables["t2"]
    u10 = nc.variables["u10"]
    v10 = nc.variables["v10"]
    tm = nc.variables["time"]
    lwdown = nc.variables["lwdown"]
    swdown = nc.variables["swdown"]
    q2 = nc.variables["q2"]
    rc = nc.variables["rain_con"]
    rn = nc.variables["rain_non"]
    lats = nc.variables["latitcrs"]
    lons = nc.variables["longicrs"]

    # keep masking in-tact as we only write data below when we have it
    # https://github.com/numpy/numpy/blob/master/doc/release/1.11.0-notes.rst#assigning-to-slicesviews-of-maskedarray
    otmpk = ma.array(ncout.variables["tmpk"][:])
    otmpk._sharedmask = False
    owmps = ma.array(ncout.variables["wmps"][:])
    owmps._sharedmask = False
    oswout = ma.array(ncout.variables["swout"][:])
    oswout._sharedmask = False
    olwout = ma.array(ncout.variables["lwout"][:])
    olwout._sharedmask = False
    oh = ma.array(ncout.variables["h"][:])
    oh._sharedmask = False
    olf = ma.array(ncout.variables["lf"][:])
    olf._sharedmask = False
    obdeckt = ma.array(ncout.variables["bdeckt"][:])
    obdeckt._sharedmask = False
    oifrost = ma.array(ncout.variables["ifrost"][:])
    oifrost._sharedmask = False
    odwpk = ma.array(ncout.variables["dwpk"][:])
    odwpk._sharedmask = False
    ofrostd = ma.array(ncout.variables["frostd"][:])
    ofrostd._sharedmask = False
    oicond = ma.array(ncout.variables["icond"][:])
    oicond._sharedmask = False
    errorcount = 0
    for i in range(len(nc.dimensions["i_cross"])):
        if errorcount > 100:
            LOG.info("Too many errors, aborting....")
            sys.exit()
        for j in range(len(nc.dimensions["j_cross"])):
            lat = lats[i, j]
            lon = lons[i, j]
            # Hey, we only care about Iowa data! -97 40 -90 43.75
            if lat < 40 or lat > 43.75 or lon < -97 or lon > -90:
                continue
            rwisfn = make_rwis(i, j, initts, oldncout, t2[1, i, j])
            modelfp = open("modeldata.txt", "w", encoding="utf8")
            for t in range(1, len(nc.dimensions["time"])):
                ts = initts + datetime.timedelta(minutes=int(tm[t]))

                modelfp.write(
                    ("%s %6.1f %6.2f %7.6f %7.2f %7.2f " "%7.2f %7.4f\n")
                    % (
                        ts.strftime("%Y-%m-%d_%H:%M:%S00"),
                        tm[t],
                        t2[t, i, j],
                        q2[t, i, j],
                        (u10[t, i, j] ** 2 + v10[t, i, j] ** 2) ** 0.5,
                        swdown[t, i, j],
                        lwdown[t, i, j],
                        rc[t, i, j] + rn[t, i, j],
                    )
                )

            modelfp.close()

            proc = subprocess.Popen(
                ("./bridgemodel %s %s %s")
                % (rwisfn, "modeldata.txt", "propertyFile"),
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            se = proc.stderr.read().decode("utf-8")
            if se != "":
                errorcount += 1
                LOG.info(
                    "bridgemodel error i:%03i j:%03i stderr:|%s|",
                    i,
                    j,
                    se.strip(),
                )
                continue
            # Process the output file!
            for line in open("pavetemp.out", encoding="utf8"):
                tokens = line.split()
                if len(tokens) < 7:
                    continue
                ts = datetime.datetime.strptime(tokens[0], "%Y-%m-%d_%H:%M")
                ts = ts.replace(tzinfo=pytz.timezone("UTC"))
                if ts.minute not in (0, 15, 30, 45) or ts < initts:
                    continue
                tidx = int(
                    int(
                        (ts - initts).days * 1400
                        + ((ts - initts).seconds / 60)
                    )
                    / 15
                )
                otmpk[tidx, i, j] = float(tokens[1])
                owmps[tidx, i, j] = float(tokens[2])
                oswout[tidx, i, j] = float(tokens[3])
                olwout[tidx, i, j] = float(tokens[4])
                oh[tidx, i, j] = float(tokens[5])
                olf[tidx, i, j] = float(tokens[6])
                if tokens[7] != "nan":
                    obdeckt[tidx, i, j] = float(tokens[7])
                oifrost[tidx, i, j] = 1 if float(tokens[8]) > 0 else 0
                ofrostd[tidx, i, j] = float(tokens[8])
                odwpk[tidx, i, j] = float(tokens[9])
                scond = tokens[-1].strip()
                if scond in CONDITIONS:
                    oicond[tidx, i, j] = CONDITIONS.index(scond)

    ncout.variables["tmpk"][:] = otmpk
    ncout.variables["wmps"][:] = owmps
    ncout.variables["swout"][:] = oswout
    ncout.variables["lwout"][:] = olwout
    ncout.variables["h"][:] = oh
    ncout.variables["lf"][:] = olf
    ncout.variables["bdeckt"][:] = obdeckt
    ncout.variables["ifrost"][:] = oifrost
    ncout.variables["frostd"][:] = ofrostd
    ncout.variables["dwpk"][:] = odwpk
    ncout.variables["icond"][:] = oicond


def find_last_output(initts):
    """See if we have a previous run on file, that can be used to spin up
    our current run"""
    for i in range(-12, -73, -12):
        ts = initts + datetime.timedelta(hours=i)
        testfn = "output/%s_iaoutput.nc" % (ts.strftime("%Y%m%d%H%M"),)
        if os.path.isfile(testfn):
            LOG.info("  Using %s as warmup values", testfn)
            return netCDF4.Dataset(testfn, "r")
    LOG.info("Did not find a previous output, will use dummy RWIS data :(")
    return None


def downsize_output(initts):
    """Subset the output file, so to save some space 66% actually"""
    fn1 = "output/%s_output.nc" % (initts.strftime("%Y%m%d%H%M"),)
    fn2 = "output/%s_iaoutput.nc" % (initts.strftime("%Y%m%d%H%M"),)
    fn3 = "/mesonet/share/frost/bridget/%s_iaoutput.nc" % (
        initts.strftime("%Y%m%d%H%M"),
    )
    if os.path.isfile(fn2):
        os.unlink(fn2)
    cmd = "ncks -d i_cross,%s,82 -d j_cross,%s,98 %s %s" % (
        IOFFSET,
        JOFFSET,
        fn1,
        fn2,
    )
    proc = subprocess.Popen(
        cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE
    )
    proc.stdout.read()
    # Make sure fn2 exists before deleting the old one
    if os.path.isfile(fn2):
        os.unlink(fn1)
        shutil.copyfile(fn2, fn3)


def main():
    """GO MAIN GO"""
    fn = sys.argv[1]
    LOG.debug("Running for %s", fn)
    nc = netCDF4.Dataset(fn)

    initts = find_initts(nc)
    ncout = make_output(nc, initts)
    oldncout = find_last_output(initts)
    run_model(nc, initts, ncout, oldncout)
    ncout.close()
    nc.close()
    downsize_output(initts)
    if os.path.isfile("faux_rwis.txt"):
        os.unlink("faux_rwis.txt")


if __name__ == "__main__":
    main()

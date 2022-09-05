#!/usr/bin/env python3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Copyright 2019 California Institute of Technology. ALL RIGHTS RESERVED.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# United States Government Sponsorship acknowledged. This software is subject to
# U.S. export control laws and regulations and has been classified as 'EAR99 NLR'
# (No [Export] License Required except when exporting to an embargoed country,
# end user, or in support of a prohibited end use). By downloading this software,
# the user agrees to comply with all applicable U.S. export laws and regulations.
# The user has the responsibility to obtain export licenses, or other export
# authority as may be required before exporting this software to any 'EAR99'
# embargoed foreign country or citizen of those countries.
#
# Authors: Piyush Agram, Yang Lei
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import math
import datetime
import numpy as np
from osgeo import gdal

def cmdLineParse():
    '''
    Command line parser.
    '''
    import argparse

    parser = argparse.ArgumentParser(description='Output geo grid')
    parser.add_argument('-m', '--input_m', dest='indir_m', type=str, required=True,
            help='Input folder with ISCE swath files for master image or master image file name (in GeoTIFF format and Cartesian coordinates)')
    parser.add_argument('-s', '--input_s', dest='indir_s', type=str, required=True,
            help='Input folder with ISCE swath files for slave image or slave image file name (in GeoTIFF format and Cartesian coordinates)')
#    parser.add_argument('-o', '--output', dest='outfile', type=str, default='geogrid.csv',
#            help='Output grid mapping')
    parser.add_argument('-d', '--dem', dest='demfile', type=str, required=True,
            help='Input DEM')
    parser.add_argument('-sx', '--dhdx', dest='dhdxfile', type=str, default="",
            help='Input slope in X')
    parser.add_argument('-sy', '--dhdy', dest='dhdyfile', type=str, default="",
            help='Input slope in Y')
    parser.add_argument('-vx', '--vx', dest='vxfile', type=str, default="",
            help='Input velocity in X')
    parser.add_argument('-vy', '--vy', dest='vyfile', type=str, default="",
            help='Input velocity in Y')
    parser.add_argument('-srx', '--srx', dest='srxfile', type=str, default="",
            help='Input search range in X')
    parser.add_argument('-sry', '--sry', dest='sryfile', type=str, default="",
            help='Input search range in Y')
    parser.add_argument('-csminx', '--csminx', dest='csminxfile', type=str, default="",
            help='Input chip size min in X')
    parser.add_argument('-csminy', '--csminy', dest='csminyfile', type=str, default="",
            help='Input chip size min in Y')
    parser.add_argument('-csmaxx', '--csmaxx', dest='csmaxxfile', type=str, default="",
            help='Input chip size max in X')
    parser.add_argument('-csmaxy', '--csmaxy', dest='csmaxyfile', type=str, default="",
            help='Input chip size max in Y')
    parser.add_argument('-ssm', '--ssm', dest='ssmfile', type=str, default="",
            help='Input stable surface mask')
    parser.add_argument('-fo', '--flag_optical', dest='optical_flag', type=bool, required=False, default=0,
            help='flag for reading optical data (e.g. Landsat): use 1 for on and 0 (default) for off')

    return parser.parse_args()

class Dummy(object):
    pass


def loadProduct(xmlname):
    '''
    Load the product using Product Manager.
    '''
    import isce
    from iscesys.Component.ProductManager import ProductManager as PM

    pm = PM()
    pm.configure()

    obj = pm.loadProduct(xmlname)

    return obj


def getMergedOrbit(product):
    import isce
    from isceobj.Orbit.Orbit import Orbit

    ###Create merged orbit
    orb = Orbit()
    orb.configure()

    burst = product[0].bursts[0]
    #Add first burst orbit to begin with
    for sv in burst.orbit:
        orb.addStateVector(sv)


    for pp in product:
        ##Add all state vectors
        for bb in pp.bursts:
            for sv in bb.orbit:
                if (sv.time< orb.minTime) or (sv.time > orb.maxTime):
                    orb.addStateVector(sv)

    return orb


def loadMetadata(indir):
    '''
    Input file.
    '''
    import os
    import numpy as np

    frames = []
    for swath in range(1,4):
        inxml = os.path.join(indir, 'IW{0}.xml'.format(swath))
        if os.path.exists(inxml):
            ifg = loadProduct(inxml)
            frames.append(ifg)

    info = Dummy()
    info.sensingStart = min([x.sensingStart for x in frames])
    info.sensingStop = max([x.sensingStop for x in frames])
    info.startingRange = min([x.startingRange for x in frames])
    info.farRange = max([x.farRange for x in frames])
    info.prf = 1.0 / frames[0].bursts[0].azimuthTimeInterval
    info.rangePixelSize = frames[0].bursts[0].rangePixelSize
    info.lookSide = -1
    info.numberOfLines = int( np.round( (info.sensingStop - info.sensingStart).total_seconds() * info.prf)) + 1
    info.numberOfSamples = int( np.round( (info.farRange - info.startingRange)/info.rangePixelSize)) + 1
    info.orbit = getMergedOrbit(frames)
    info.indir = indir

    return info


def coregisterLoadMetadataOptical(indir_m, indir_s):
    '''
    Input file.
    '''
    import os
    import numpy as np

    from osgeo import gdal, osr
    import struct
    import re

    import isce
    from components.contrib.geo_autoRIFT.geogrid import GeogridOptical
#    from geogrid import GeogridOptical

    obj = GeogridOptical()

    x1a, y1a, xsize1, ysize1, x2a, y2a, xsize2, ysize2, trans = obj.coregister(indir_m, indir_s)

    DS = gdal.Open(indir_m, gdal.GA_ReadOnly)

    info = Dummy()
    info.startingX = trans[0]
    info.startingY = trans[3]
    info.XSize = trans[1]
    info.YSize = trans[5]

    if re.findall("L[CO]08_",DS.GetDescription()).__len__() > 0:
        nameString = os.path.basename(DS.GetDescription())
        info.time = nameString.split('_')[3]
    elif re.findall("L[EO]07_",DS.GetDescription()).__len__() > 0:
         #pdb.set_trace()
        nameString = os.path.basename(DS.GetDescription())
        info.time = nameString.split('_')[3]
    elif re.findall("S2._",DS.GetDescription()).__len__() > 0:
        info.time = DS.GetDescription().split('_')[2]
    else:
        raise Exception('Optical data NOT supported yet!')

    info.numberOfLines = ysize1
    info.numberOfSamples = xsize1

    info.filename = indir_m

    DS1 = gdal.Open(indir_s, gdal.GA_ReadOnly)

    info1 = Dummy()

    if re.findall("L[CO]08_",DS1.GetDescription()).__len__() > 0:
        nameString1 = os.path.basename(DS1.GetDescription())
        info1.time = nameString1.split('_')[3]
    elif re.findall("L[EO]07_",DS1.GetDescription()).__len__() > 0:
        nameString1 = os.path.basename(DS1.GetDescription())
        info1.time = nameString1.split('_')[3]
    elif re.findall("S2._",DS1.GetDescription()).__len__() > 0:
        info1.time = DS1.GetDescription().split('_')[2]
    else:
        raise Exception('Optical data NOT supported yet!')

    return info, info1

def getProjectionSystem(filename):
    '''
    Testing with Greenland.
    '''
    import subprocess
    import re

    if not filename:
        raise Exception('File {0} does not exist'.format(filename))

    from osgeo import gdal, osr
    ds = gdal.Open(filename, gdal.GA_ReadOnly)
    srs = osr.SpatialReference()
    srs.ImportFromWkt(ds.GetProjection())
    srs.AutoIdentifyEPSG()
    ds = None

    if srs.IsGeographic():
        epsgstr = srs.GetAuthorityCode('GEOGCS')
    elif srs.IsProjected():
        epsgstr = srs.GetAuthorityCode('PROJCS')
    elif srs.IsLocal():
        raise Exception('Local coordinate system encountered')
    else:
        raise Exception('Non-standard coordinate system encountered')
    if not epsgstr:  #Empty string->use shell command gdalsrsinfo for last trial
        cmd = 'gdalsrsinfo -o epsg {0}'.format(filename)
        epsgstr = subprocess.check_output(cmd, shell=True)
        epsgstr = re.findall("EPSG:(\d+)", str(epsgstr))[0]
    if not epsgstr:  #Empty string
        raise Exception('Could not auto-identify epsg code')
    epsgcode = int(epsgstr)
    return epsgcode


def determineBbox(self, zrange=[-200,4000]):
    '''
    Dummy.
    '''
    import numpy as np
    import datetime
    from osgeo import osr,gdal
    
#        import pdb
#        pdb.set_trace()

#        rng = self.startingRange + np.linspace(0, self.numberOfSamples, num=21)
    rng = self.startingRange + np.linspace(0, self.numberOfSamples-1, num=21) * self.rangePixelSize
    deltat = np.linspace(0, 1., num=21)[1:-1]

    lonlat = osr.SpatialReference()
    lonlat.ImportFromEPSG(4326)

    coord = osr.SpatialReference()
    coord.ImportFromEPSG(getProjectionSystem(self.demname))
    
    trans = osr.CoordinateTransformation(lonlat, coord)

    llhs = []
    xyzs = []


    ###First range line
    for rr in rng:
        for zz in zrange:
            llh = self.orbit.rdr2geo(self.sensingStart, rr, side=self.lookSide, height=zz)
            llhs.append(llh)
            if gdal.__version__[0] == '2':
                x,y,z = trans.TransformPoint(llh[1], llh[0], llh[2])
            else:
                x,y,z = trans.TransformPoint(llh[0], llh[1], llh[2])
            xyzs.append([x,y,z])

    ##Last range line
    sensingStop = self.sensingStart + datetime.timedelta(seconds = (self.numberOfLines-1) / self.prf)
    for rr in rng:
        for zz in zrange:
            llh = self.orbit.rdr2geo(sensingStop, rr, side=self.lookSide, height=zz)
            llhs.append(llh)
            if gdal.__version__[0] == '2':
                x,y,z = trans.TransformPoint(llh[1], llh[0], llh[2])
            else:
                x,y,z = trans.TransformPoint(llh[0], llh[1], llh[2])
            xyzs.append([x,y,z])


    ##For each line in middle, consider the edges
    for frac in deltat:
        sensingTime = self.sensingStart + datetime.timedelta(seconds = frac * (self.numberOfLines-1)/self.prf)
#            print('sensing Time: %f %f %f'%(sensingTime.minute,sensingTime.second,sensingTime.microsecond))
        for rr in [rng[0], rng[-1]]:
            for zz in zrange:
                llh = self.orbit.rdr2geo(sensingTime, rr, side=self.lookSide, height=zz)
                llhs.append(llh)
                if gdal.__version__[0] == '2':
                    x,y,z = trans.TransformPoint(llh[1], llh[0], llh[2])
                else:
                    x,y,z = trans.TransformPoint(llh[0], llh[1], llh[2])
                xyzs.append([x,y,z])


    llhs = np.array(llhs)
    xyzs = np.array(xyzs)

    xlim = [np.min(xyzs[:,0]), np.max(xyzs[:,0])]
    ylim = [np.min(xyzs[:,1]), np.max(xyzs[:,1])]
    return xlim, ylim

def getazimuthAngle(filename):
    from zipfile import ZipFile
    
    with ZipFile(filename, 'r') as zipObj:
        listOfiles = zipObj.namelist()
        for l in listOfiles:
            if l.endswith('xml') and l.split('/')[-2]=='annotation':
                f = zipObj.open(l)
                fl = str(f.read())
                azimuth_angle = float(fl.split('platformHeading>')[1][:-2])
                break
    return azimuth_angle


def generateFile(obj, xlim, ylim):
    from iscesys import DateTimeUtil as DTU
    import xml.etree.ElementTree as ET
    
    tree = ET.parse(obj.ref_fn[::-1].split('/', 1)[-1][::-1] + '/reference.xml')
    ref_dir = str(tree.getroot().findall(".//*[@name='safe']")[0].text[1:-1])
    azimuth_angle = getazimuthAngle(ref_dir)
    grd_res = obj.rangePixelSize / math.sin(obj.incidenceAngle)
    sensingStart = DTU.seconds_since_midnight(obj.sensingStart)
    tmid = obj.sensingStart + datetime.timedelta(seconds = (np.floor(obj.numberOfLines/2)-1) / obj.prf)
    satvmid = obj.orbit.interpolateOrbit(tmid)._velocity
    azm_res = np.linalg.norm(satvmid) / obj.prf
    lines = []
    lines.append("Radar parameters:")
    lines.append("Range: " + "{:.0f}".format(obj.startingRange) + " {:.5f}".format(obj.rangePixelSize))
    lines.append("Azimuth: " + "{:.3f}".format(sensingStart) + " {:.3f}".format(obj.prf))
    lines.append("Dimensions: " + str(obj.numberOfSamples) + " " + str(obj.numberOfLines))
    lines.append("Incidence Angle: " + "{:.4f}".format(obj.incidenceAngle*180/np.pi))
    lines.append("Azimuth angle: " + "{:.4f}".format(360-azimuth_angle))

    lines.append("\nMap inputs: ")
    lines.append("EPSG: " + str(getProjectionSystem(obj.demname)))
    lines.append("Smallest Allowable Chip Size in m: " + str(obj.chipSizeX0))
    lines.append("Grid spacing in m: " + "{:.4f}".format(obj.gridSpacingX))
    lines.append("Repeat Time: " + "{:.4e}".format(obj.repeatTime))
    lines.append("XLimits: " + "{:.0f}".format(xlim[0]) + "  " + "{:.0f}".format(xlim[1]))
    lines.append("YLimits: " + "{:.5e}".format(ylim[0]) + "  " + "{:.5e}".format(ylim[1]))
    lines.append("Extent in km: " + "{:.3f}".format((xlim[1]-xlim[0])/1000.0) + "  " + "{:.3f}".format((ylim[1]-ylim[0])/1000.0))
    lines.append("DEM: " + obj.demname)
    lines.append("Slopes: " + obj.dhdxname + " " + obj.dhdyname)
    lines.append("Velocities: /DATA/glacier-vel/geogrid_req/ref_velx.tif  /DATA/glacier-vel/geogrid_req/ref_vely.tif")
    lines.append("Stable Surface Mask: /DATA/glacier-vel/geogrid_req/ssm.tif")

    lines.append("\nOutputs:")
    lines.append("Window locations: window_location.tif")
    lines.append("Window rdr_off2vel_x vector: window_rdr_off2vel_x_vec.tif")
    lines.append("Window rdr_off2vel_y vector: window_rdr_off2vel_y_vec.tif")
    lines.append("Window chip size min: window_chip_size_min.tif")
    lines.append("Window chip size max: window_chip_size_max.tif")
    lines.append("Window stable surface mask: window_stable_surface_mask.tif")
    lines.append("Output Nodata Value: " + str(obj.nodata_out) + "\n")

    demDS = gdal.Open(obj.demname, gdal.GA_ReadOnly)

    geoTrans = demDS.GetGeoTransform()
    demXSize = demDS.RasterXSize
    demYSize = demDS.RasterYSize

    lOff = int(np.max( [np.floor((ylim[1] - geoTrans[3])/geoTrans[5]), 0.]))
    lCount = int(np.min([ np.ceil((ylim[0] - geoTrans[3])/geoTrans[5]), demYSize-1.]) - lOff)

    pOff = int(np.max([ np.floor((xlim[0] - geoTrans[0])/geoTrans[1]), 0.]))
    pCount = int(np.min([ np.ceil((xlim[1] - geoTrans[0])/geoTrans[1]), demXSize-1.]) - pOff)

    lines.append("Starting processing ....")
    lines.append("Xlimits : " + str(round(geoTrans[0] + pOff * geoTrans[1])) +  "  " + str(round(geoTrans[0] + (pOff + pCount) * geoTrans[1])))
    lines.append("Ylimits : " + "{:.5e}".format(geoTrans[3] + (lOff + lCount) * geoTrans[5]) +  "  " + "{:.5e}".format(round(geoTrans[3] + lOff * geoTrans[5])))
    lines.append("Origin index (in DEM) of geogrid: " + str(pOff) + "   " + str(lOff))

    lines.append("Dimensions of geogrid: " + str(pCount) + " x " + str(lCount))
    lines.append("Ground range pixel size: " + "{:.6}".format(grd_res))
    lines.append("Azimuth pixel size: " + "{:.6}".format(azm_res))
    with open('testGeogrid.txt', 'w') as f:
        for line in lines:
            f.write(line)
            f.write('\n')


def runGeogrid(info, info1, dem, dhdx, dhdy, vx, vy, srx, sry, csminx, csminy, csmaxx, csmaxy, ssm, **kwargs):
    '''
    Wire and run geogrid.
    '''

    import isce
    from components.contrib.geo_autoRIFT.geogrid import Geogrid
#    from geogrid import Geogrid

    from osgeo import gdal
    dem_info = gdal.Info(dem, format='json')

    obj = Geogrid()
    obj.configure()

    obj.startingRange = info.startingRange
    obj.rangePixelSize = info.rangePixelSize
    obj.sensingStart = info.sensingStart
    obj.prf = info.prf
    obj.lookSide = info.lookSide
    print(obj.startingRange, obj.rangePixelSize, obj.sensingStart)
    obj.repeatTime = (info1.sensingStart - info.sensingStart).total_seconds()
    obj.numberOfLines = info.numberOfLines
    obj.numberOfSamples = info.numberOfSamples
    obj.nodata_out = -32767
    obj.chipSizeX0 = 240
    obj.gridSpacingX = dem_info['geoTransform'][1]
    obj.orbit = info.orbit
    obj.demname = dem
    obj.dhdxname = dhdx
    obj.dhdyname = dhdy
    obj.vxname = vx
    obj.vyname = vy
    obj.srxname = srx
    obj.sryname = sry
    obj.csminxname = csminx
    obj.csminyname = csminy
    obj.csmaxxname = csmaxx
    obj.csmaxyname = csmaxy
    obj.ssmname = ssm
    obj.winlocname = "window_location.tif"
    obj.winoffname = "window_offset.tif"
    obj.winsrname = "window_search_range.tif"
    obj.wincsminname = "window_chip_size_min.tif"
    obj.wincsmaxname = "window_chip_size_max.tif"
    obj.winssmname = "window_stable_surface_mask.tif"
    obj.winro2vxname = "window_rdr_off2vel_x_vec.tif"
    obj.winro2vyname = "window_rdr_off2vel_y_vec.tif"
    obj.ref_fn = info.indir
    ##dt-varying search range scale (srs) rountine parameters
#    obj.srs_dt_unity = 5
#    obj.srs_max_scale = 10
#    obj.srs_max_search = 20000
#    obj.srs_min_search = 0

    obj.getIncidenceAngle()
    xlim, ylim = determineBbox(obj)
    generateFile(obj, xlim, ylim)
    obj.geogrid()
    with open('testGeogrid.txt', 'a') as f:
        line = "Scene-center lat/lon: " + str(obj.cen_lat) + "  " + str(obj.cen_lon)
        f.write(line)

    run_info = {
        'chipsizex0': obj.chipSizeX0,
        'gridspacingx': obj.gridSpacingX,
        'vxname': vx,
        'vyname': vy,
        'sxname': kwargs.get('dhdxs'),
        'syname': kwargs.get('dhdys'),
        'maskname': kwargs.get('sp'),
        'xoff': obj.pOff,
        'yoff': obj.lOff,
        'xcount': obj.pCount,
        'ycount': obj.lCount,
        'dt': obj.repeatTime,
        'epsg': kwargs.get('epsg'),
        'XPixelSize': obj.X_res,
        'YPixelSize': obj.Y_res,
        'cen_lat': obj.cen_lat,
        'cen_lon': obj.cen_lon,
    }


    return run_info


def runGeogridOptical(info, info1, dem, dhdx, dhdy, vx, vy, srx, sry, csminx, csminy, csmaxx, csmaxy, ssm, **kwargs):
    '''
    Wire and run geogrid.
    '''

    import isce
    from components.contrib.geo_autoRIFT.geogrid import GeogridOptical
#    from geogrid import GeogridOptical

    from osgeo import gdal
    dem_info = gdal.Info(dem, format='json')

    obj = GeogridOptical()

    obj.startingX = info.startingX
    obj.startingY = info.startingY
    obj.XSize = info.XSize
    obj.YSize = info.YSize
    from datetime import date
    import numpy as np
    d0 = date(np.int(info.time[0:4]),np.int(info.time[4:6]),np.int(info.time[6:8]))
    d1 = date(np.int(info1.time[0:4]),np.int(info1.time[4:6]),np.int(info1.time[6:8]))
    date_dt_base = d1 - d0
    obj.repeatTime = date_dt_base.total_seconds()
#    obj.repeatTime = (info1.time - info.time) * 24.0 * 3600.0
    obj.numberOfLines = info.numberOfLines
    obj.numberOfSamples = info.numberOfSamples
    obj.nodata_out = -32767
    obj.chipSizeX0 = 240
    obj.gridSpacingX = dem_info['geoTransform'][1]

    obj.dat1name = info.filename
    obj.demname = dem
    obj.dhdxname = dhdx
    obj.dhdyname = dhdy
    obj.vxname = vx
    obj.vyname = vy
    obj.srxname = srx
    obj.sryname = sry
    obj.csminxname = csminx
    obj.csminyname = csminy
    obj.csmaxxname = csmaxx
    obj.csmaxyname = csmaxy
    obj.ssmname = ssm
    obj.winlocname = "window_location.tif"
    obj.winoffname = "window_offset.tif"
    obj.winsrname = "window_search_range.tif"
    obj.wincsminname = "window_chip_size_min.tif"
    obj.wincsmaxname = "window_chip_size_max.tif"
    obj.winssmname = "window_stable_surface_mask.tif"
    obj.winro2vxname = "window_rdr_off2vel_x_vec.tif"
    obj.winro2vyname = "window_rdr_off2vel_y_vec.tif"
    ##dt-varying search range scale (srs) rountine parameters
#    obj.srs_dt_unity = 32
#    obj.srs_max_scale = 10
#    obj.srs_max_search = 20000
#    obj.srs_min_search = 0

    obj.runGeogrid()

    run_info = {
        'chipsizex0': obj.chipSizeX0,
        'gridspacingx': obj.gridSpacingX,
        'vxname': vx,
        'vyname': vy,
        'sxname': kwargs.get('dhdxs'),
        'syname': kwargs.get('dhdys'),
        'maskname': kwargs.get('sp'),
        'xoff': obj.pOff,
        'yoff': obj.lOff,
        'xcount': obj.pCount,
        'ycount': obj.lCount,
        'dt': obj.repeatTime,
        'epsg': kwargs.get('epsg'),
        'XPixelSize': obj.X_res,
        'YPixelSize': obj.Y_res,
        'cen_lat': obj.cen_lat,
        'cen_lon': obj.cen_lon,
    }

    return run_info

def main():
    '''
    Main driver.
    '''

    inps = cmdLineParse()

    if inps.optical_flag == 1:
        metadata_m, metadata_s = coregisterLoadMetadataOptical(inps.indir_m, inps.indir_s)
        runGeogridOptical(metadata_m, metadata_s, inps.demfile, inps.dhdxfile, inps.dhdyfile, inps.vxfile, inps.vyfile, inps.srxfile, inps.sryfile, inps.csminxfile, inps.csminyfile, inps.csmaxxfile, inps.csmaxyfile, inps.ssmfile)
    else:
        metadata_m = loadMetadata(inps.indir_m)
        metadata_s = loadMetadata(inps.indir_s)
        runGeogrid(metadata_m, metadata_s, inps.demfile, inps.dhdxfile, inps.dhdyfile, inps.vxfile, inps.vyfile, inps.srxfile, inps.sryfile, inps.csminxfile, inps.csminyfile, inps.csmaxxfile, inps.csmaxyfile, inps.ssmfile)


if __name__ == '__main__':
    main()

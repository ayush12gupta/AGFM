# Glacier Flow Analysis using SAR Data

This repo contains the code for setting up the folder structure as well as some of the postprocessing steps. For feature tracking we have used [autoRIFT](https://github.com/nasa-jpl/autoRIFT.git) and for co-registration of SAR Images we have used [ISCE](https://github.com/isce-framework/isce2), we have provided instuctions on how to [install ISCE](docs/ISCE Installation.pdf).

Follow instuctions given in [*INSTRUCTION.md*](docs/INSTRUCTION.md) for setting up the ISCE file structure for running ***topsApp.py***. Once the co-registeration has been completed it would generate ***./merged/reference.slc.full*** and ***./merged/secondary.slc.full***, which are merged and co-registered reference and secondary images respectively.

### **Running Geogrid for Convertion to Image coordinates** 
Now, since the pre-processing is completed, we use Geogrid for mapping data from geographic coordinates to image coordinates, which takes in following parameters as input:-

* Reference (-m) and Coregistered Secondary Image (-s)
* SRTM DEM i.e. Digital Elevation Model (-d)
* Local surface slopes along X and Y directions separately (-sx and -sy)
* Pixel wise maximum and minimum chip size for the region [OPTIONAL] (-csmaxx/y and -csminx/y)
* Reference velocity i.e. velocity results from some other study [OPTIONAL] (-vx and -vy)

```
usage: geogrid_autorift/testGeogrid_ISCE.py [-h] -m INDIR_M -s INDIR_S -d DEMFILE
                           [-sx DHDXFILE] [-sy DHDYFILE] [-vx VXFILE]
                           [-vy VYFILE] [-csminx CSMINXFILE] 
                           [-csminy CSMINYFILE] [-csmaxx CSMAXXFILE]
                           [-csmaxy CSMAXYFILE] [-ssm SSMFILE] 
                           [-fo OPTICAL_FLAG]

Output geo grid

optional arguments:
  -h, --help            show this help message and exit
  -m INDIR_M, --input_m INDIR_M
                        Input folder with ISCE swath files for master image or
                        master image file name (in GeoTIFF format and
                        Cartesian coordinates)
  -s INDIR_S, --input_s INDIR_S
                        Input folder with ISCE swath files for slave image or
                        slave image file name (in GeoTIFF format and Cartesian
                        coordinates)
  -d DEMFILE, --dem DEMFILE
                        Input DEM
  -sx DHDXFILE, --dhdx DHDXFILE
                        Input slope in X
  -sy DHDYFILE, --dhdy DHDYFILE
                        Input slope in Y
  -vx VXFILE, --vx VXFILE
                        Input velocity in X
  -vy VYFILE, --vy VYFILE
                        Input velocity in Y
  -csminx CSMINXFILE, --csminx CSMINXFILE
                        Input chip size min in X
  -csminy CSMINYFILE, --csminy CSMINYFILE
                        Input chip size min in Y
  -csmaxx CSMAXXFILE, --csmaxx CSMAXXFILE
                        Input chip size max in X
  -csmaxy CSMAXYFILE, --csmaxy CSMAXYFILE
                        Input chip size max in Y
  -ssm SSMFILE, --ssm SSMFILE
                        Input stable surface mask
  -fo OPTICAL_FLAG, --flag_optical OPTICAL_FLAG
                        flag for reading optical data (e.g. Landsat): use 1
                        for on and 0 (default) for off

```

Here we would have to convert all the inputs from WGS84 to UTM projection, it can be done using this command:
```bash
> gdalwarp -s_srs "EPSG:4326" -t_srs '+proj=utm +zone=43 +datum=WGS84 +units=m +no_defs' -of GTIFF INPUT_PATH OUTPUT_TIF_PATH
```
For the case of Himalayas it lies in zone 43 of UTM projection, adjust that value accordingly.

We can directly use the DEM that is automatically downloaded by ISCE during pre-processing step, and convert it to UTM projection. And the local slopes along X and Y directions seperately, can be computed using ***DirectionalSlope*** plugin of **QGIS**.

And for chip size if you have a reference velocity ( i.e. results from some other study ) for that region using, then you can allocate the chip size according to it, i.e. setting smaller chip min and max chip size for region having higher reference velocity, and for stable region ( velocity <15 m/yr ) we can set larger min and max chip sizes. However, if you don't have a reference velocity, it can be ignored as in that case, default min and max chip size would be used uniformy for all the regions

For running the geogrid script, ***geogrid_autorift/testGeogrid_ISCE.py*** is used. An example command calling testGeogrid_ISCE.py has been given below.
```bash
> python geogrid_autorift/testGeogrid_ISCE.py -m fine_coreg -s secondary -d demtest_roi.tif -sx demtest_roi_X.tif -sy demtest_roi_Y.tif -vx ref_velx.tif -vy ref_vely.tif -ssm ssm.tif
```

Output files may include all or few
```
"winlocname":        range/azimuth pixel indices (2-band; in units of integer image pixels), 
"winoffname":        downstream search (expected) range/azimuth pixel displacement (2-band; in units of integer image pixels), 
"winsrname":         range/azimuth search range (2-band; in units of integer image pixels), 
"wincsminname":      range/azimuth chip size minimum (2-band; in units of integer image pixels), 
"wincsmaxname":      range/azimuth chip size maximum (2-band; in units of integer image pixels), 
"winssmname":        stable surface mask (boolean), 
"winro2vxname":      2-by-1 conversion coefficients from radar range/azimuth displacement to x-direction motion velocity (3-band; 3rd band is conversion coefficient from range pixel displacement to range motion velocity), 
"winro2vyname":      2-by-1 conversion coefficients from radar range/azimuth displacement to y-direction motion velocity (3-band; 3rd band is conversion coefficient from azimuth pixel displacement to azimuth motion velocity). 
```

***NOTE:*** After execution of testGeogrid_ISCE.py copy the console outputs into testGeogrid.txt file, as the data from that file would be used by autoRIFT.

### **Running Geogrid for Convertion to Image coordinates** 

Now the only step left is feature tracking using autoRIFT (Autonomous Repeat Image Feature Tracking), it takes the reference and secondary images from ISCE merged output, as well as other local slope image etc,generated from Geogrid step as input.

```
usage: geogrid_autorift/testautoRIFT_ISCE.py [-h] -m INDIR_M -s INDIR_S [-g GRID_LOCATION]
                            [-o INIT_OFFSET] [-sr SEARCH_RANGE]
                            [-csmin CHIP_SIZE_MIN] [-csmax CHIP_SIZE_MAX]
                            [-vx OFFSET2VX] [-vy OFFSET2VY]
                            [-ssm STABLE_SURFACE_MASK] [-fo OPTICAL_FLAG]
                            [-nc NC_SENSOR] [-mpflag MPFLAG] [-ncname NCNAME]

Output geo grid

optional arguments:
  -h, --help            show this help message and exit
  -m INDIR_M, --input_m INDIR_M
                        Input master image file name (in ISCE format and radar
                        coordinates) or Input master image file name (in
                        GeoTIFF format and Cartesian coordinates)
  -s INDIR_S, --input_s INDIR_S
                        Input slave image file name (in ISCE format and radar
                        coordinates) or Input slave image file name (in
                        GeoTIFF format and Cartesian coordinates)
  -g GRID_LOCATION, --input_g GRID_LOCATION
                        Input pixel indices file name
  -csmin CHIP_SIZE_MIN, --input_csmin CHIP_SIZE_MIN
                        Input chip size min file name
  -csmax CHIP_SIZE_MAX, --input_csmax CHIP_SIZE_MAX
                        Input chip size max file name
  -vx OFFSET2VX, --input_vx OFFSET2VX
                        Input pixel offsets to vx conversion coefficients file
                        name
  -vy OFFSET2VY, --input_vy OFFSET2VY
                        Input pixel offsets to vy conversion coefficients file
                        name
  -ssm STABLE_SURFACE_MASK, --input_ssm STABLE_SURFACE_MASK
                        Input stable surface mask file name
  -fo OPTICAL_FLAG, --flag_optical OPTICAL_FLAG
                        flag for reading optical data (e.g. Landsat): use 1
                        for on and 0 (default) for off
  -nc NC_SENSOR, --sensor_flag_netCDF NC_SENSOR
                        flag for packaging output formatted for Sentinel ("S")
                        and Landsat ("L") dataset; default is None
  -mpflag MPFLAG, --mpflag MPFLAG
                        number of threads for multiple threading (default is
                        specified by 0, which uses the original single-core
                        version and surpasses the multithreading routine)
  -ncname NCNAME, --ncname NCNAME
                        User-defined filename for the NetCDF output to which
                        the ROI percentage and the production version will be
                        appended

```

Here -mpflag sets the no. of threads for systems with multiple cores in order to reduce inference time using parallel computing. And -ncname specifies the name of the ***NetCDF*** file that would be generated at the end of feature tracking algorithm, which would store all the meta data of the reference and secondary image, as well as the computed velocity in azimuthal and range direction.

An example command running ***geogrid_autorift/testautoRIFT_ISCE.py*** for performing feature tracking is given below:
```bash
> python ./geogrid_autorift/testautoRIFT_ISCE.py -m merged/reference.slc.full -s merged/secondary.slc.full -g window_location.tif -o window_offset.tif -vx window_rdr_off2vel_x_vec.tif -vy window_rdr_off2vel_y_vec.tif -ssm window_stable_surface_mask.tif -mpflag 128 -ncname exp2_nc
```

### **Post-processing**

For post-processing, we have used velocity in LOS and azimuth direction for computing velocity in flow direction, and rate of change of thickness of glacier, and visualised them in order to get a idea of spatial variation of flow.

The postprocessing function can be changed from ***geogrid_autorift/util.py***

We have performed all the steps in this [jupyter notebook](Post-Processing.ipynb).


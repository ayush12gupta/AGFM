# Glacier Flow Analysis using SAR Data

This repo contains the code for fully automated SAR based Offset tracking which can be used in various applicaitons we have some demos showing Offset tracking for glacier velocity estimation. For feature tracking we have used [autoRIFT](https://github.com/nasa-jpl/autoRIFT.git) and for co-registration of SAR Images we have used [ISCE](https://github.com/isce-framework/isce2), we have provided instuctions on how to [install ISCE](docs/ISCE Installation.pdf).

We support processing a batch of SAR Image for velocity estimation as well.

### **Running Single Image Pair velocity estimation** 
The overall pipeline of single task processing is shown in figure below. We just need to provide the reference and secondary image pair URL, which can be obtained from [ASF Platform](https://search.asf.alaska.edu/#). Script ***single_process.py*** is used for single image pair velocity estimation, it takes in following parameters as input:-

* Reference Image URL (-reference) and Secondary Image URL(--secondary)
* Save path i.e. directory in which files will be saved (--save_path)
* Velocity output file names (--netCDF_out)
* DATA Config file containing coregistration parameters [DEFAULT](./configs/isce_config.json) (--config)

An example command calling single_process.py has been given below.
       
    python single_process.py --reference REFERENCE_URL --secondary SECONDARY_URL --save_path OUT_PATH --netCDF_out POST_FILENAME

```
usage: single_process.py [-h] [--reference REFERENCE] [--secondary SECONDARY]
                         [--save_path SAVE_PATH] [--netCDF_out NETCDF_OUT]
                         [--config CONFIG]

optional arguments:
  -h, --help            show this help message and exit
  --reference REFERENCE
                        URL of reference zip file
  --secondary SECONDARY
                        URL of secondary zip file
  --save_path SAVE_PATH
                        directory in which orbit file needs to be saved
  --netCDF_out NETCDF_OUT
                        Output netCDF file
  --config CONFIG       Data config file

```

***NOTE:*** The offset tracking chipsize and num of threads used can be changed from data_config.json, and the parameters for coregisteration can be changed from isce_config.json changing the Region of Interest. 

During coregisteration ISCE automatically downloads the DEM files which are then processed later automatically to generate DEM and Slope files for ROI.

After autoRIFT algorithm is performed on the coregistered Master and Slave image, the offset maps are then used to compute velocities in radar amd azimuth directions. After, which postprocessing can then be applied using a custom post-processing function defined in ***geogrid_autorift/util.py***. 


Output files may include all or few
```
"velocity.tif":           velocity in azimuth and range directions,
"offset.tif":             displacement offsets in azimuth and range directions, 
"testGeogrid.txt":        metadata of SAR image pairs, 
"Post-processed Output":  outputs generated from custom post-processing function, 
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


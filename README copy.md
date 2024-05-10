# Automated Glacier Flow Monitoring using SAR Data

This repo contains the code for fully automated SAR based glacier monitoring pipeline which generates a 12 days separated time-series velocity maps.

This pipeline currently only supports the use of Sentinel-1 data. 

## Introduction

To comprehend glacier dynamics for a region, a time-series study of glacier change is essential. However, generation of a large time series data often requires a substantial amount of computation and time. To address these limitations, we developed an efficient pipeline for processing of such large time-series Sentinel-1 imagery, generating extensive time series data for 3D glacier flow velocities.

For feature tracking we developed a robut offset tracking module, built on top of [autoRIFT](https://github.com/nasa-jpl/autoRIFT.git). It has been modified for performing offset tracking based on NCC stacking, using the time-series of coregistered SAR imagery. For co-registration of SAR Images we have used stack processing mode of [ISCE](https://github.com/isce-framework/isce2), please go through the instructions on how to [install ISCE](https://github.com/isce-framework/isce2/blob/main/README.md).

The overall pipeline of single task processing is shown in figure below. 

![Pipeline for Velocity Estimation](./docs/overall_pipeline.png)


## Running the pipeline

We just need to provide the list of scene names for both ascending and descending track Sentinel-1 images. They can be generated using SSARA as demonstrated in [`Pre-processing.ipynb`](/notebook/Pre-processing.ipynb). And other parameters are given through [configs files](/configs/). 

An example for [pipeline config](/configs/pipeline_config.json) is given,
```json
{
    "save_path": "/DATA/run_2017_19", # Workflow run directory
    "SAR_dir": "/DATA/S2_Data/",  # Sentinel-1 file save directiory
    "shapefile_dir": "/DATA/shapefiles/CB_glacier_buffer.shp",  # Glacier region shapefile
    "config_path": "/DATA/Automated_Offset_Tracking/configs/data_config.json",  # Data config file path
    "polarisation": "vv"  # Sentinel-1 polarisation to be used
}
```

Similarly, an example for [data config](/configs/data_config.json) is given,
```json
{
    "num_threads": 64,  # No. of cores to be used
    "chip_min": 240,    # Minimum chip size for offset tracking
    "chip_max": 960,    # Maximum chip size for offset tracking
    "Orbit_dir": "/DATA/S2_Data/orbit/",  # Orbit file save directory
    "aux_dir": "/DATA/S2_Data/aux/",  # Auxilary file save directory
    "cred_config": "config/credentials.json",  # Credentials json file
    "ROI": "[32.06, 32.60, 77.09, 77.82]"  # Area of interest
}

```
Script [***stack_pipeline.py***](/stack_pipeline.py) is used for running the pipeline, it takes in following parameters as input:-
* Ascending track scene list (-t_asc)
* Descending track scene list (-t_des)
* Pipeline config file (--config)

An example command calling single_process.py has been given below.

```bash
    python single_process.py --reference REFERENCE_URL --secondary SECONDARY_URL --save_path OUT_PATH --netCDF_out POST_FILENAME
```

```
usage: stack_pipeline.py [-h] -t_asc DOWNLOAD_ASC_TXT -t_des DOWNLOAD_DES_TXT [--config CONFIG]

optional arguments:
  -h, --help            show this help message and exit
  -t_asc DOWNLOAD_ASC_TXT, --download_asc_txt DOWNLOAD_ASC_TXT
                        Data Ascending txt file
  -t_des DOWNLOAD_DES_TXT, --download_des_txt DOWNLOAD_DES_TXT
                        Data Descending txt file
  --config CONFIG       Pipeline config file

```

### Running Single Image Pair velocity estimation
We just need to provide the reference and secondary image pair URL, which can be obtained from [ASF Platform](https://search.asf.alaska.edu/#). Script ***single_process.py*** is used for single image pair velocity estimation, it takes in following parameters as input:-
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
                        Sentinel1 scene name of reference zip file
  --secondary SECONDARY
                        Sentinel1 scene name of secondary zip file
  --save_path SAVE_PATH
                        directory in which orbit file needs to be saved
  --netCDF_out NETCDF_OUT
                        Output netCDF file
  --config CONFIG       Data config file

```

***NOTE:*** The offset tracking chipsize and num of threads used can be changed from data_config.json, and the parameters for coregisteration can be changed from isce_config.json changing the Region of Interest. 

During coregisteration ISCE automatically downloads the DEM files which are then processed later automatically to generate DEM and Slope files for ROI.

After autoRIFT algorithm is performed on the coregistered Master and Slave image, the offset maps are then used to compute velocities in radar amd azimuth directions. After, which postprocessing can then be applied using a custom post-processing function defined in ***offset_tracking/util.py***. 


Output files may include all or few
```
"velocity.tif":           velocity in azimuth and range directions,
"offset.tif":             displacement offsets in azimuth and range directions, 
"testGeogrid.txt":        metadata of SAR image pairs, 
"Post-processed Output":  outputs generated from custom post-processing function, 
```

***NOTE:*** After execution of testGeogrid_ISCE.py copy the console outputs into testGeogrid.txt file, as the data from that file would be used by autoRIFT.

### **Running Batch of Image Pair for velocity estimation** 

For estimating velocity maps for a batch of images, we take in a [CSV file](./data/data_download1.csv) as input which contains:

* Pair Name
* Master and Slave Sentinel1 scene name 
* Region of Interest (lat/lon) --> Eg: [32.06, 32.77, 76.86, 77.82]
* Status

Here the batch process is divided in 2 process and can be run in parallel using script ***isce_batch.py*** and ***offset_batch.py***. The "status" in CSV provides the information regarding at which stage the process is at, Status=1 when coregisteration is completed and its waiting fot velocity estimation, at the same time the ones which are not able to complete have Status=-1. Similarly once velocity estimation is completed Status is changed to 2, and if its unsucessful its changed to -2.

For starting the execution of coregistration pipeline ***isce_batch.py*** should be used as shown below.

    python3 isce_batch.py --download_csv CSV_DATA --save_path OUTPUT_PATH --config CONFIG_PATH

```
usage: isce_batch.py [-h] [--save_path SAVE_PATH]
                     [--download_csv DOWNLOAD_CSV] [--config CONFIG]

optional arguments:
  -h, --help            show this help message and exit
  --save_path SAVE_PATH
                        directory in which orbit file needs to be saved
  --download_csv DOWNLOAD_CSV
                        Data CSV file
  --config CONFIG       Data config file
```

Parallel to executing coregisteration pipeline, the offset tracking pipeline can also be started which will check for Status=1 Image pairs and start computing velocity maps for them.

For starting the execution of offset tracking pipeline ***offset_batch.py*** should be used as shown below.

    python3 offset_batch.py --download_csv CSV_DATA --save_path OUTPUT_PATH --config CONFIG_PATH


```
usage: offset_batch.py [-h] [--save_path SAVE_PATH]
                       [--download_csv DOWNLOAD_CSV] [--config CONFIG]

optional arguments:
  -h, --help            show this help message and exit
  --save_path SAVE_PATH
                        directory in which orbit file needs to be saved
  --download_csv DOWNLOAD_CSV
                        Data CSV file
  --config CONFIG       Data config file
```


### **Running Stack Processing of Image Pair along with velocity correction**

For performing time-series analysis over a region, a stack of SAR images are choosen in [TXT file](./data/stack_data.txt). All the SAR images are coregistered with respect to the first time-step image, and we perform the the offset tracking on the image pairs that we selected automatically in ***image_pair.csv***. Once offset tracking is completed, the velocity maps generated is modelled into an overdetermined linear system by minimising the least square using SVD method.

The stack processing is performed using ***stack_process.py***. For starting the execution it should be called as shown below:

    python stack_process.py --save_path SAVE_DIR --aux AUX_DIR --data_path DATA_PATH --config CONFIG_PATH --download_txt DOWNLOAD_TXT

```
usage: stack_process.py [-h] [--save_path SAVE_PATH] [--aux AUX]
                        [--data_path DATA_PATH] [--download_txt DOWNLOAD_TXT]
                        [--config CONFIG]

optional arguments:
  -h, --help            show this help message and exit
  --save_path SAVE_PATH
                        directory in which orbit file needs to be saved
  --aux AUX             Sentinel-1 AUX file directory where aux file is saved
  --data_path DATA_PATH
                        directory in which orbit file needs to be saved
  --download_txt DOWNLOAD_TXT
                        Data CSV file
  --config CONFIG       Data config file
```

### **Post-processing**

For post-processing, we have used velocity in LOS and azimuth direction for computing velocity in flow direction, and rate of change of thickness of glacier, and visualised them in order to get a idea of spatial variation of flow.

The postprocessing function can be changed from ***geogrid_autorift/util.py***

We have performed all the steps in this [jupyter notebook](Post-Processing.ipynb).

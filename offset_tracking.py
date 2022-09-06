import os
from utils import *

def offset_tracking(config, cwd, netCDF_out):

    os.makedirs("offset_tracking", exist_ok=True)
    execute('rm -rf geom_reference/ fine_offsets/ fine_interferogram/ coarse_* ESD/ PICKLE/')
    execute('rm merged/z.rdr.full.* merged/lat.rdr.full.* merged/los.rdr.full.* merged/lon.rdr.full.* merged/pwr.bil* demLat_N31_N34_Lon_E076_E079.dem*')
    os.chdir("./offset_tracking")
    if not os.path.exists('window_location.tif'):
        execute(f'python {cwd}/geogrid_autorift/testGeogrid_ISCE.py -m ../fine_coreg -s ../secondary -d {config["dem"]["DEM"]} -sx {config["dem"]["X"]} -sy {config["dem"]["Y"]} -ssm {config["ssm"]}')
    else:
        print("./window_location.tif  ---> Already exists")
    
    if os.path.exists('../merged/reference.slc.full')&os.path.exists('../merged/secondary.slc.full'):
        cmd = f'time python {cwd}/geogrid_autorift/testautoRIFT_ISCE.py -m ../merged/reference.slc.full -s ../merged/secondary.slc.full -g window_location.tif -vx window_rdr_off2vel_x_vec.tif -vy window_rdr_off2vel_y_vec.tif -ssm window_stable_surface_mask.tif -mpflag {str(config["num_threads"])} -config {cwd}/configs/velocity_config.json -ncname {netCDF_out}'   # -nc S'
        execute(cmd)
    else:
        print("Coregistered images not available check isce.log for issues with ISCE coregisteration")
    
    execute('rm window*')
    print("Offset Tracking completed")
    
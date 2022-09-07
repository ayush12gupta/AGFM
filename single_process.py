import os
import json
import argparse

from utils import *
from offset_tracking import offset_tracking


'''
Make sure that reference is of earlier date as compared to secondary
'''
parser = argparse.ArgumentParser()
parser.add_argument('--reference', type=str, default=None, help="URL of reference zip file")
parser.add_argument('--secondary', type=str, default=None, help="URL of secondary zip file")
parser.add_argument('--save_path', type=str, default="./outptut", help="directory in which orbit file needs to be saved")
parser.add_argument('--netCDF_out', type=str, default="out", help="Output netCDF file")
parser.add_argument('--config', type=str, default="./configs/data_config.json", help="Data config file")

args = parser.parse_args()


def main():

    # Reading config file
    with open(args.config, 'r') as f:
        config = json.load(f)

    # Reading isce config file
    with open('configs/isce_config.json', 'r') as f:
        isce_config = json.load(f)

    os.makedirs(args.save_path, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(args.save_path)

    # Setting up ISCE enivonment for coregisteration
    execute(f'python3 {cwd}/preprocessing/setup_env.py --config {cwd}/configs/isce_config.json --reference {args.reference} --secondary {args.secondary} --orbit_path {config["Orbit_data"]} --data_pathR {config["SAR_dir"]} --data_pathS {config["SAR_dir"]}')
    ## NOTE: Document the changes in mergebursts to incorporate multilooking of .slc.full

    if not os.path.exists('merged/secondary.slc.full'):
        execute('time topsApp.py topsApp.xml --start=startup --end=mergebursts')

    dem_dirs = [d for d in os.listdir('./') if (d[-5:]=='wgs84')&(d[:6]=='demLat')]
    generate_dem_products(dem_dirs[0], isce_config['ROI'])
    
    if os.path.exists('merged/secondary.slc.full'):
        execute(f'python {cwd}/geogrid_autorift/topsinsar_filename.py')
        offset_tracking(config, cwd, args.netCDF_out)
        os.chdir(cwd)
    else:
        os.chdir(cwd)
        print("Some issue with ISCE coregisteration. Please look at log. Exiting.....")
        raise SystemExit()

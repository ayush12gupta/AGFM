import os
import json
import argparse
import pandas as pd

from utils import execute
from offset_tracking import offset_tracking


parser = argparse.ArgumentParser()
parser.add_argument('--save_path', type=str, default="./outptut", help="directory in which orbit file needs to be saved")
parser.add_argument('--download_csv', type=str, default="./data/data_download.csv", help="Data CSV file")
parser.add_argument('--config', type=str, default="./configs/data_config.json", help="Data config file")

args = parser.parse_args()


def main():

    # Reading config file
    with open(args.config, 'r') as f:
        config = json.load(f)

    os.makedirs(args.save_path, exist_ok=True)
    cwd = os.getcwd()

    data_pairs = pd.read_csv(args.download_csv, header=0)
    for i, row in data_pairs.iterrows():
        os.chdir(args.save_path)
        year = str(row['Start Date'][:4])
        os.makedirs(year, exist_ok=True)
        os.chdir(year)
        os.makedirs(str(row['Start Date']), exist_ok=True)
        os.chdir(str(row['Start Date']))
        execute(f'python3 {cwd}/preprocessing/setup_env.py --reference {config["SAR_dir"]}/{year} --secondary {config["SAR_dir"]}/{year} --orbit_path {config["Orbit_data"]} --data_pathR {config["SAR_dir"]}/{year} --data_pathS {config["SAR_dir"]}/{year}')
        
        if not os.path.exists('merged/secondary.slc.full'):
            execute('time topsApp.py topsApp.xml --start=startup --end=mergebursts')
        
        if os.path.exists('merged/secondary.slc.full'):
            print(os.getcwd(), "Coregisteration complete")
            data_pairs.at[i,'Status'] = 1
            data_pairs.to_csv(args.download_csv, index=False)
        else:
            data_pairs.at[i,'Status'] = -1  # Error in ISCE
            print("Some issue with ISCE coregisteration. Check for status -1 in CSV file")

        data_pairs.to_csv(args.download_csv, index=False)

    os.chdir(cwd)
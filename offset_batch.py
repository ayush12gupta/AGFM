import os, glob
import json
import time
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

    cwd = os.getcwd()

    data_pairs = pd.read_csv(args.download_csv, header=0)
    while ((data_pairs['Status']==0).sum()+(data_pairs['Status']==1).sum()):

        data_pairs = pd.read_csv(args.download_csv, header=0)
        
        if len(data_pairs[data_pairs['Status']==1])==0:
            time.sleep(60)
            continue
        
        row = data_pairs[data_pairs['Status']==1].iloc[0]
        index = data_pairs[data_pairs['Start Date']==row['Start Date']].index[0]        
        os.chdir(args.save_path)
        year = str(row['Start Date'][:4])
        os.chdir(os.path.join(year, str(row['Start Date'])))
        print(os.getcwd())
        if os.path.exists('merged/secondary.slc.full'):
            offset_tracking(config, cwd, str(row['Start Date']).replace('-',''))
        
        data_pairs = pd.read_csv(args.download_csv, header=0)
        if os.path.exists('velocity.tif'):
            data_pairs.at[index,'Status'] = 2
        else:
            data_pairs.at[index,'Status'] = -2
            print("Some issues with Offset Tracking")
        
        data_pairs.to_csv(args.download_csv, index=False)
    
    os.chdir(cwd)


if __name__=="__main__":
    main()
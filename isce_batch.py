import os
import glob
import json
import argparse
import pandas as pd

from utils import execute, generate_dem_products


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
    # for i, row in data_pairs[100:].iterrows():
    while ((data_pairs['Status']==0).sum()):
        
        data_pairs = pd.read_csv(args.download_csv, header=0)
        row = data_pairs[data_pairs['Status']==0].iloc[0]
        index = data_pairs[data_pairs['Start Date']==row['Start Date']].index[0]

        if row['Status']!=0:
            continue

        os.chdir(args.save_path)
        year = str(row['Start Date'][:4])
        os.makedirs(year, exist_ok=True)
        os.chdir(year)
        os.makedirs(str(row['Start Date']), exist_ok=True)
        os.chdir(str(row['Start Date']))
        execute(f'python3 {cwd}/preprocessing/setup_env.py --roi "{row["ROI"]}" --reference {row["Master"]} --secondary {row["Slave"]} --orbit_path {config["Orbit_dir"]} --data_pathR {config["SAR_dir"]}/{year} --data_pathS {config["SAR_dir"]}/{year}')
        
        if not os.path.exists('merged/secondary.slc.full'):
            execute('cp -r /DATA/glacier-vel/geogrid_req/dem/demLat_N31_N34_Lon_E076_E079* ./')
            out = execute(f'time {cwd}/topsApp.py topsApp.xml --start=startup --end=mergebursts')
            if out==0:
                data_pairs = pd.read_csv(args.download_csv, header=0)
                data_pairs.at[index,'Status'] = -1  # Error in ISCE
                print("Some issue with ISCE coregisteration. Check for status -1 in CSV file")
                continue

        if os.path.exists('merged/secondary.slc.full'):
            print(os.getcwd(), "Coregisteration complete")
            dem_dirs = glob.glob('./demLat*wgs84')
            # dem_dirs = [d for d in os.listdir('./') if (d[-5:]=='wgs84')&(d[:6]=='demLat')]
            generate_dem_products(dem_dirs[0], row['ROI'])
            execute(f'python3 {cwd}/geogrid_autorift/topsinsar_filename.py')
                
            data_pairs = pd.read_csv(args.download_csv, header=0)
            data_pairs.at[index,'Status'] = 1
            data_pairs.to_csv(args.download_csv, index=False)
        else:
            data_pairs = pd.read_csv(args.download_csv, header=0)
            data_pairs.at[index,'Status'] = -1  # Error in ISCE
            print("Some issue with ISCE coregisteration. Check for status -1 in CSV file")

        data_pairs.to_csv(args.download_csv, index=False)

    os.chdir(cwd)


if __name__=="__main__":
    main()
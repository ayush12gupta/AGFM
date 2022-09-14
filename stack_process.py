import os, shutil
import glob
import time
import json
import argparse
import pandas as pd
import xml.etree.ElementTree as ET
from preprocessing.get_orbit import get_orbit_fl
from utils import execute, generate_dem_products


parser = argparse.ArgumentParser()
parser.add_argument('--save_path', type=str, default="./output", help="directory in which orbit file needs to be saved")
parser.add_argument('--aux', type=str, default="/DATA/glacier-vel/data/aux/", help="Sentinel-1 AUX file directory where aux file is saved")
parser.add_argument('--data_path', type=str, default="/DATA/glacier-vel/stack/", help="directory in which orbit file needs to be saved")
parser.add_argument('--download_txt', type=str, default="./data/data_download.csv", help="Data CSV file")
parser.add_argument('--config', type=str, default="./configs/data_config.json", help="Data config file")

args = parser.parse_args()


def download_data(username, password, url, data_path, orbit_path):
    print(url.split('/')[-1].split('_'))
    dmy= url.split('/')[-1].split('_')[5]
    sentinel_type = url.split('/')[-1].split('_')[0]
    date, month, yr = dmy[6:8], dmy[4:6], dmy[:4]
    date = "%02d" % (int(date))
    rt = os.system("wget " + url + " -P " + data_path + " --user=" + username + " --password=" + password + " -nc")
    if rt!=0:
        print("Failed downloading ", url)
    get_orbit_fl(sentinel_type, date, month, yr, orbit_path)
    return url.split('/')[-1]


def offset_tracking(csv_file):

    pairs = pd.read_csv(csv_file, index=0)


def main():

    # Reading config files
    with open(args.config, 'r') as f:
        config = json.load(f)
    with open('./configs/isce_config.json', 'r') as f:
        config_isce = json.load(f)

    urls = []
    with open(args.download_txt, 'r') as f:
        for line in f:
            urls.append(line[:-1])

    # if os.path.exists(args.data_path):
    #     shutil.rmtree(args.data_path)
    # os.mkdir(args.data_path)
    cwd = os.getcwd()
    
    # for url in urls:
    #     download_data(config_isce['ASF_user'], config_isce['ASF_password'], url, args.data_path, config['Orbit_dir'])

    if os.path.exists(args.save_path):
        shutil.rmtree(args.save_path)
    os.mkdir(args.save_path)
    os.chdir(args.save_path)
    if len(glob.glob('merged/SLC/*/*.slc.full'))==0:
        print(config_isce["ROI"][1:-1].replace(',',''))
        execute('cp -r /DATA/glacier-vel/geogrid_req/dem/demLat_N31_N34_Lon_E076_E079* ./')
        dem = glob.glob('*.wgs84')
        execute(f'python3 {cwd}/stack/topsStack/stackSentinel.py -s {args.data_path} -d {dem[0]} -a {args.aux} -o {config["Orbit_dir"]} -b "{config_isce["ROI"][1:-1].replace(",","")}" -t "python3 {cwd}/stack/topsStack/" -W slc')

        run_files = glob.glob('run_files/*')
        flag = 0
        strt_time = time.time()
        for file in sorted(run_files):
            if flag:
                print("Coregisteration complete")
                break
            if 'merge' in file:
                flag = 1 

            execute(f'bash ./{file}')
            print(f"Time taken till {file}: ", time.time()-strt_time)
        print("Time taken for stack coregisteration: ", time.time()-strt_time)

    # Preparing for geogrid
    dem = glob.glob('*.wgs84')
    generate_dem_products(dem[0], config_isce["ROI"])
    xmlCoreg = glob.glob('coreg_secondarys/*/IW*.xml')
    for xml in xmlCoreg:
        tree = ET.parse(xml)
        root = tree.getroot()
        for state in root.findall('component'):
            if state.find('factoryname').text=='coregSwathSLCProduct':
                if state.find('factorymodule').text=='coregSwathSLCProduct':
                    state.find('factorymodule').text = 'contrib.stack.topsStack.coregSwathSLCProduct'

        tree.write(xml)
    
    os.makedirs('./offset_tracking', exist_ok=True)
    os.chdir("./offset_tracking")

    pair_fn = '../image_pairs.csv'
    pairs = pd.read_csv('../image_pairs.csv', index=0)
    # for i, row in pairs.iterrows():
    #     if os.path.exists(f'merged/SLC/{row[]}')

    secondarys = sorted(os.listdir('../coreg_secondarys'))

    if not os.path.exists('window_location.tif'):
        execute(f'python {cwd}/geogrid_autorift/testGeogrid_ISCE.py -m ../reference -s ../coreg_secondarys/{secondarys[1]} -d ../dem_crop.tif -sx ../dem_x.tif -sy ../dem_x.tif')   #  -ssm {config["ssm"]}
    else:
        print("./window_location.tif  ---> Already exists")

    print("Starting Offset Tracking")

    # offset_tracking(pair_fn)

if __name__=='__main__':
    main()
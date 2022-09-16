import os, shutil
import glob
import time
import json
import argparse
import pandas as pd
import numpy as np
from osgeo import gdal
import xml.etree.ElementTree as ET
from preprocessing.get_orbit import get_orbit_fl
from geogrid_autorift.util import numpy_array_to_raster
from utils import execute, generate_dem_products, get_deltaT, read_vals, get_DT


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


def add_elements(elements, acquisitionDates, gap):
    N = len(acquisitionDates)
    id = [f'stack_{str(elements[i])}_{str(elements[i+gap])}' for i in range(N-gap)]
    master = [acquisitionDates[i] for i in range(N-gap)]
    slave = [acquisitionDates[i+gap] for i in range(N-gap)]
    return id, master, slave


def generate_data(elements, acquisitionDates, max_gap=3):
    id = []
    master = []
    slave = []
    for gap in range(1, max_gap+1):
        idn, mastern, slaven = add_elements(elements, acquisitionDates, gap=gap)
        id.extend(idn)
        master.extend(mastern)
        slave.extend(slaven)

    return id, master, slave


def offset_tracking(config, cwd, master, slave, netCDF_out):

    execute('rm -rf ../coreg_secondarys ../ESD ../misreg')
    execute('cp ../testGeogrid.txt testGeogrid.txt')
    if os.path.exists(f'../../merged/SLC/{master}/{master}.slc.full')&os.path.exists(f'../../merged/SLC/{slave}/{slave}.slc.full'):
        cmd = f'time python {cwd}/geogrid_autorift/testautoRIFT_ISCE.py -m ../../merged/SLC/{master}/{master}.slc.full -s ../../merged/SLC/{slave}/{slave}.slc.full -g ../window_location.tif -chipmax {config["chip_max"]} -chipmin {config["chip_min"]} -mpflag {str(config["num_threads"])} -config {cwd}/configs/data_config.json -ncname {netCDF_out} -vx ../window_rdr_off2vel_x_vec.tif -vy ../window_rdr_off2vel_y_vec.tif -post' #-ssm window_stable_surface_mask.tif -nc S'
        execute(cmd)
    else:
        print("Coregistered images not available check isce.log for issues with ISCE coregisteration")
    
    print("Offset Tracking completed")


def offset_compute(csv_file, config, cwd):

    data_pairs = pd.read_csv(csv_file, header=0)
    while (data_pairs['Status']==0).sum():

        row = data_pairs[data_pairs['Status']==0].iloc[0]
        index = data_pairs[data_pairs['Id']==row['Id']].index[0]
        id  = row['Id']
        master, slave = row['Master'], row['Slave']
        os.makedirs(f'./{id}', exist_ok=True)
        os.chdir(f'./{id}')
        if os.path.exists(f'../../merged/SLC/{slave}/{slave}.slc.full'):
            offset_tracking(config, cwd, str(master), str(slave), id)

        data_pairs = pd.read_csv(os.path.join('../', csv_file), header=0)
        if os.path.exists('velocity.tif'):
            data_pairs.at[index,'Status'] = 1
        else:
            data_pairs.at[index,'Status'] = -1
            print("Some issues with Offset Tracking")

        data_pairs.to_csv(os.path.join('../', csv_file), index=False)
        os.chdir('../')

    print("Offset Tracking completed")


def compute_svd(A, L, deltaT, dT, nodata, num, N):
    shp = L[0].shape
    A = (A*deltaT)
    L = np.array(L).reshape(num, -1)*deltaT[0]
    U, s, VT = np.linalg.svd(A, full_matrices=True)
    S_inv = np.zeros(A.shape).T
    S_inv[:N-1,:N-1] = np.diag(1/s[:N-1])
    x_svd = ((VT.T.dot(S_inv).dot(U.T))@L) #.reshape(N, shp[0], shp[1])
    La = (A@x_svd).reshape(-1, num)/dT
    La = La.reshape(num, shp[0], shp[1])
    La[:,nodata] = -32767
    return La


def velocity_correction_band(csv_fn, tif_dir, dates, band=1, max_gap=3):
    data = pd.read_csv(csv_fn, header=0)
    num = len(data)
    sg = max_gap*(max_gap+1)/2
    assert (num + sg)%max_gap==0
    N = ((num + sg)//max_gap) - 1
    deltaT = get_deltaT(sorted(dates))
    A = np.zeros((num, N))
    nodata = None
    L = []
    Ids = []
    dT = []

    # till = 3
    # N = till
    # num = (3*(till+1))-6
    # A = np.zeros((num, till))
    # print(A.shape)
    # deltaT = deltaT[:till]
    idx = 0
    for i, row in data.iterrows():
        id = row['Id'].split('_')[1:]
        # if int(id[1:])>till:
        #     continue
        st, ed = int(id[0]), int(id[1])
        A[idx, st:ed] = 1
        val, nodata = read_vals(os.path.join(tif_dir, row['Id'], 'velocity.tif'), nodat=nodata, band=band)
        dT.append(get_DT(row['Master'], row['Slave']))
        L.append(val)
        Ids.append(row['Id'])
        idx += 1
    
    dT = np.array(dT)
    La = compute_svd(A, L, deltaT, dT, nodata, num, N)
    # To get projections and transformations
    ds = gdal.Open(os.path.join(tif_dir, data['Id'][0], 'velocity.tif'))
    projs = ds.GetProjection()
    geo = ds.GetGeoTransform()
    ds = None

    for i in range(num):
        ds = numpy_array_to_raster(f'{Ids[i]}_{band}.tif', np.array([La[i]]), projs, geo, nband=1)
        ds.FlushCache()
        ds = None


def velocity_correction(csv_fn, tif_dir, dates):
    '''
    Performs SVD for least square estimation of velocity correction of 
    offset tracking results

    csv_file: Path to csv file 'image_pair.csv'
    tif_dir: Path to offset_tracking tif file results
    dates: list of acquisition interval
    '''
    velocity_correction_band(csv_fn, tif_dir, dates, band=1)
    velocity_correction_band(csv_fn, tif_dir, dates, band=2)
    # Merging both the bands and additionally generating resulting velocity band
    data = pd.read_csv(csv_fn, header=0)

    # To get projections and transformations
    ds = gdal.Open(f"{data['Id'][0]}_1.tif")
    projs = ds.GetProjection()
    geo = ds.GetGeoTransform()
    ds = None

    for i, row in data.iterrows():
        val1, _ = read_vals(f"{row['Id']}_1.tif", nodat=None, band=1)
        val2, _ = read_vals(f"{row['Id']}_2.tif", nodat=None, band=1)
        nodata1 = (val1==-32767)
        nodata2 = (val2==-32767) 
        val3 = np.sqrt(val1**2 + val2**2)
        val3[nodata1|nodata2] = -32767

        ds = numpy_array_to_raster(f"{row['Id']}.tif", np.array([val3, val1, val2]), projs, geo, nband=3)
        ds.FlushCache()
        ds = None

        if os.path.exists(f"{row['Id']}.tif"):
            execute(f"rm {row['Id']}_*.tif")

    print("Velocity Correction completed")


def main():
    # Reading config files
    with open(args.config, 'r') as f:
        config = json.load(f)
    with open('./configs/isce_config.json', 'r') as f:
        config_isce = json.load(f)
    bbox = config_isce["ROI"][1:-1].replace(",","")

    urls = []
    with open(args.download_txt, 'r') as f:
        for line in f:
            urls.append(line[:-1])

    # if os.path.exists(args.data_path):
    #     shutil.rmtree(args.data_path)
    os.makedirs(args.data_path, exist_ok=True)
    cwd = os.getcwd()
    
    for url in urls:
        download_data(config_isce['ASF_user'], config_isce['ASF_password'], url, args.data_path, config['Orbit_dir'])

    # if os.path.exists(args.save_path):
    #     shutil.rmtree(args.save_path)
    os.makedirs(args.save_path, exist_ok=True)
    os.chdir(args.save_path)

    acquisitionDates = [url.split('_')[2] for url in urls]
    N = len(acquisitionDates)
    elements = [i for i in range(N)]
    id, master, slave = generate_data(elements, acquisitionDates, 2)
    df = pd.DataFrame({'Id':id, 'Master':master, 'Slave':slave, 'Status': [0]*len(id), 'ROI':[str(f"[{bbox.replace(' ', ', ')}]")]*len(id)})
    df.to_csv('image_test.csv')

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
    generate_dem_products(dem[0], config_isce["ROI"], config)
    
    dates = os.listdir('merged/SLC/')
    os.makedirs('./offset_tracking', exist_ok=True)
    os.chdir("./offset_tracking")

    pair_fn = '../image_pairs.csv'
    secondarys = sorted(os.listdir('../secondarys'))
    if not os.path.exists('window_location.tif'):
        execute(f'python {cwd}/geogrid_autorift/testGeogrid_ISCE.py -m ../reference -s ../secondarys/{secondarys[0]} -d ../dem_crop.tif -sx ../dem_x.tif -sy ../dem_x.tif')   #  -ssm {config["ssm"]}
    else:
        print("./window_location.tif  ---> Already exists")

    print("Starting Offset Tracking")
    offset_compute(pair_fn, config, cwd)
    
    print('Starting Velocity Correction ...')
    os.chdir('../')
    if os.path.exists('./velocity_corrected'):
        shutil.rmtree('./velocity_corrected')
        os.makedirs('./velocity_corrected')
    os.chdir('./velocity_corrected')
    velocity_correction(pair_fn, '../offset_tracking/', dates)
    os.chdir('../')

    execute('rm -rf ./geom_reference ./merged/geom_reference')


if __name__=='__main__':
    main()
    # urls = []
    # with open('./data/stack_data.txt', 'r') as f:
    #     for line in f:
    #         urls.append(line[:-1])
    
    # acquisitionDates = [url.split('_')[2] for url in urls]
    # N = len(acquisitionDates)
    # elements = [i for i in range(N)]
    # id, master, slave = generate_data(elements, acquisitionDates, 2)
    # df = pd.DataFrame({'Id':id, 'Master':master, 'Slave':slave, 'Status': [0]*len(id), 'ROI':[str(f"[{inps.bbox.replace(' ', ', ')}]")]*len(id)})

    # print(id, master, slave, len(id), N)
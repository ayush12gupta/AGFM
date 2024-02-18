import os, shutil
import glob
import time
import json
import argparse
import pandas as pd
import numpy as np
import subprocess, sys
from osgeo import gdal
import xml.etree.ElementTree as ET
from preprocessing.get_orbit import get_orbit_fl
from geogrid_autorift.util import numpy_array_to_raster
from utils import execute, generate_dem_products, get_deltaT, read_vals, get_DT, process_check_running
from datetime import datetime

VELOCITY_CORR_GAP = 3 ## HARDCODED the number of redundant observations
DIR_PATH = os.path.realpath(os.path.dirname(__file__))

def cmdLineParse():
    '''
    Command line parser.
    '''
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--save_path', type=str, default="./output", help="directory in which orbit file needs to be saved")
    parser.add_argument('--aux', type=str, default="/DATA/S2_Data/aux/", help="Sentinel-1 AUX file directory where aux file is saved")
    parser.add_argument('-d', '--data_path', type=str, default="/DATA/stack_data/", help="directory in which orbit file needs to be saved")
    parser.add_argument('-t', '--download_txt', type=str, default="./data/data_download.csv", help="Data CSV file")
    parser.add_argument('--config', type=str, default="./configs/data_config.json", help="Data config file")
    parser.add_argument('-m', '--mode', dest="mode", type=str, default='offset', help="Option: ['download', 'coreg', 'offset', 'post']")
    parser.add_argument('-shp', '--glacier_shp', dest='shpfile', type=str, default=None, help='Glacier Shapefile path')
    parser.add_argument('-msk', '--mask', dest='mask', type=str, default=None, help='Glacier Shapefile path')
    parser.add_argument('-p', '--polarization', type=str, default="vv", help="Polarisation to be used")
# parser.add_argument('--dem', dest="dem_path", type=str, required=True, help="Path for DEM file")
    return parser.parse_args()


def to_datetime(date):
    date = f'{date[:4]}/{date[4:6]}/{date[-2:]}'
    date = datetime.strptime(date, '%Y/%m/%d')
    return date

def execute_write(cmd, file, file_err):
    f = open(file, "w")
    f_err = open(file_err, "w")
    subprocess.check_call(cmd, shell=True, stdout=f, stderr=f_err)


def write_time(filename, exec_time):
    with open(filename, 'w') as f:
        for i, time in enumerate(exec_time):
            f.write(f'{i}: {time/60}\n')


def download_data(username, password, id, data_path, orbit_path):
    url = 'https://sentinel1.asf.alaska.edu/SLC/S{0}/{1}.zip'.format(id[2], id)
    print(url.split('/')[-1].split('_'))
    dmy= url.split('/')[-1].split('_')[5]
    sentinel_type = url.split('/')[-1].split('_')[0]
    date, month, yr = dmy[6:8], dmy[4:6], dmy[:4]
    date = "%02d" % (int(date))
    
    sentinel_file = f'{data_path}/{url.split("/")[-1]}'
    if not os.path.exists(f'{data_path}/{url.split("/")[-1]}'):
        rt = os.system("wget " + url + " -P " + data_path + " --user=" + username + " --password=" + password + " -nc")
        if rt!=0:
            print("Failed downloading ", url)
            
    # get_orbit_fl(sentinel_type, date, month, yr, orbit_path)
    os.system(f'eof --sentinel-file {sentinel_file} --save-dir {orbit_path} --asf-user {username} --asf-password {password}')
    return url.split('/')[-1]


def generate_data(elements, acquisitionDates, max_gap=3):
    id = []
    master = []
    slave = []
    steps = []
    
    for gap in range(1, max_gap+1):
        idn, mastern, slaven, stepsn = add_elements(elements, acquisitionDates, gap=gap)
        id.extend(idn)
        master.extend(mastern)
        slave.extend(slaven)
        steps.extend(stepsn)

    return id, master, slave, steps
    

def add_elements(elements, acquisitionDates, gap, max_step=3):
    N = len(acquisitionDates)
    id, master, slave, steps = [], [], [], []
    for i in range(N-gap):
        master_dt = to_datetime(acquisitionDates[i])
        slave_dt = to_datetime(acquisitionDates[i+gap])
        delta = (slave_dt - master_dt).days
        if (gap>1)&(delta>50):
            continue
        
        id.append(f'stack_{str(elements[i])}_{str(elements[i+gap])}')
        master.append(acquisitionDates[i])
        slave.append(acquisitionDates[i+gap])
        step0 = [acquisitionDates[i], acquisitionDates[i+gap]]
        dt1 = to_datetime(acquisitionDates[i])
        dt2 = to_datetime(acquisitionDates[min(i+max_step, N-1)])
        delt = (dt2 - dt1).days
        if delt==36:
            step0.extend(acquisitionDates[i:i+max_step+1])

        step0 = sorted(list(set(step0)))
        if (gap==max_step)&(delt==36):
            steps.append(','.join(step0[::3]))
        else:
            steps.append(','.join(step0))
            
    return id, master, slave, steps


def offset_tracking(config, cwd, master, slave, steps, mask, deltaT, step_dT=1):

    input_list = ','.join([f'../../merged/SLC/{st}/{st}.slc.full' for st in steps.split(',')])
    execute('cp ../testGeogrid.txt testGeogrid.txt')
    if os.path.exists(f'../../merged/SLC/{master}/{master}.slc.full')&os.path.exists(f'../../merged/SLC/{slave}/{slave}.slc.full'):
        # if (master2=='')|(slave2==''):
        #     cmd = f'time python {DIR_PATH}/offset_tracking/testautoRIFT_ISCE2.py -m ../../merged/SLC/{master}/{master}.slc.full -s ../../merged/SLC/{slave}/{slave}.slc.full -g ../window_location.tif -vx ../window_rdr_off2vel_x_vec.tif -vy ../window_rdr_off2vel_y_vec.tif -mpflag {str(config["num_threads"])} --mask {mask} -chipmin 240 -chipmax 960 -dT {int(deltaT)} -step_dT {int(step_dT)}'
        # else:
        #     cmd = f'time python {DIR_PATH}/offset_tracking/testautoRIFT_ISCE2.py -m ../../merged/SLC/{master}/{master}.slc.full,../../merged/SLC/{master2}/{master2}.slc.full -s ../../merged/SLC/{slave}/{slave}.slc.full,../../merged/SLC/{slave2}/{slave2}.slc.full -g ../window_location.tif -vx ../window_rdr_off2vel_x_vec.tif -vy ../window_rdr_off2vel_y_vec.tif -mpflag {str(config["num_threads"])} --mask {mask} -chipmin 240 -chipmax 960 -dT {int(deltaT)} -step_dT {int(step_dT)}'
        
        cmd = f'time python {DIR_PATH}/offset_tracking/testautoRIFT_ISCE2.py -i {input_list} -g ../window_location.tif -vx ../window_rdr_off2vel_x_vec.tif -vy ../window_rdr_off2vel_y_vec.tif -mpflag {str(config["num_threads"])} --mask {mask} -chipmin 240 -chipmax 960 -dT {int(deltaT)} -step_dT {int(step_dT)}'
        print('PO Command', cmd)
        cmd_file = open(f'offset_cmd.txt', 'a')
        cmd_file.write(cmd)
        cmd_file.close()
        execute(cmd)
    else:
        print("Coregistered images not available check isce.log for issues with ISCE coregisteration")
    
    print("Offset Tracking completed")


def offset_compute(csv_file, config, cwd, mask):

    # execute('rm -rf ../coreg_secondarys ../ESD ../misreg ../coreg_secondarys ../coarse_interferograms ../geom_reference')
    data_pairs = pd.read_csv(csv_file, header=0)
    print("Starting offset tracking estimation using: ", csv_file)
    while (data_pairs['Status']==0).sum():
        data_pairs = data_pairs.fillna('')
        row = data_pairs[data_pairs['Status']==0].iloc[0]
        index = data_pairs[data_pairs['Id']==row['Id']].index[0]
        id  = row['Id']
        id1, id2 = id.split('_')[1:]
        print(id)
        
        master, slave, steps = row['Master'], row['Slave'], row['Steps']
        
        date_m = str(master)
        date_s = str(slave)
        date_m = f'{date_m[:4]}/{date_m[4:6]}/{date_m[-2:]}'
        date_m = datetime.strptime(date_m, '%Y/%m/%d')
        date_s = f'{date_s[:4]}/{date_s[4:6]}/{date_s[-2:]}'
        date_s = datetime.strptime(date_s, '%Y/%m/%d')
        deltaT = abs((date_s-date_m).days)
        # Skip offset tracking for the file which have already been computed
        if os.path.exists(f'./{id}/velocity.tif'):
            data_pairs = pd.read_csv(csv_file, header=0)
            data_pairs.at[index,'Status'] = 1
            data_pairs.to_csv(csv_file, index=False)
            continue

        os.makedirs(f'./{id}', exist_ok=True)
        os.chdir(f'./{id}')
        
        # If coregistered image exists
        if os.path.exists(f'../../merged/SLC/{slave}/{slave}.slc.full'):
            offset_tracking(config, cwd, str(master), str(slave), str(steps), mask, deltaT, int(id2)-int(id1))
            data_pairs.at[index,'Status'] = 1
        data_pairs = pd.read_csv(os.path.join('../', csv_file), header=0)
        
        if os.path.exists('velocity.tif'):
            data_pairs.at[index,'Status'] = 1
        else:
            data_pairs.at[index,'Status'] = -1
            print("Some issues with Offset Tracking")

        data_pairs.to_csv(os.path.join('../', csv_file), index=False)
        os.chdir('../')

    print("Offset Tracking completed")


def velocity_postprocess(csv_file, shpfile):

    # execute('rm -rf ../coreg_secondarys ../ESD ../misreg')
    data_pairs = pd.read_csv(csv_file, header=0)
    
    print("Starting velocity postprocess estimation using: ", csv_file)
    print("Remaining Pairs:", (data_pairs['Status']==1).sum())
    while (data_pairs['Status']==1).sum():
        row = data_pairs[data_pairs['Status']==1].iloc[0]
        index = data_pairs[data_pairs['Id']==row['Id']].index[0]
        id  = row['Id']
        master, slave = row['Master'], row['Slave']
        # Skip offset tracking for the file which have already been computed
        if os.path.exists(f'./{id}/snr_res_crop.tif') and os.path.exists(f'./{id}/velocity_res_crop.tif'):
            data_pairs = pd.read_csv(csv_file, header=0)
            data_pairs.at[index,'Status'] = 2
            data_pairs.to_csv(csv_file, index=False)
            continue

        os.makedirs(f'./{id}', exist_ok=True)
        os.chdir(f'./{id}')
        print(os.getcwd())
        print(os.listdir('./'), os.path.exists(f'./velocity.tif'))
        # If coregistered image exists
        if os.path.exists(f'./velocity.tif'):
            date_m = str(master)
            date_s = str(slave)
            date_m = f'{date_m[:4]}/{date_m[4:6]}/{date_m[-2:]}'
            date_m = datetime.strptime(date_m, '%Y/%m/%d')
            date_s = f'{date_s[:4]}/{date_s[4:6]}/{date_s[-2:]}'
            date_s = datetime.strptime(date_s, '%Y/%m/%d')
            factor = int(abs((date_m - date_s).days))//12
            # if shpfile is not None:
            #     cmd = f'python {DIR_PATH}/geogrid_autorift/postprocessing.py -d ./ -f {factor} -shp {shpfile}'
            # else:
            cmd = f'python {DIR_PATH}/geogrid_autorift/postprocessing.py -d ./ -f {factor}'
            print(cmd)
            execute(cmd)
            data_pairs.at[index,'Status'] = 2
        data_pairs = pd.read_csv(os.path.join('../', csv_file), header=0)
        if os.path.exists('velocity_res.tif'):
            data_pairs.at[index,'Status'] = 2
        else:
            data_pairs.at[index,'Status'] = -2
            print("Some issues with Velocity Postprocessing")

        data_pairs.to_csv(os.path.join('../', csv_file), index=False)
        os.chdir('../')

    print("Offset Tracking completed")


def compute_svd(A, L, deltaT, dT, nodata, num, N):
    shp = L[0].shape
    A = (A*deltaT)
    L = np.array(L).reshape(num, -1)*deltaT[0]
    U, s, VT = np.linalg.svd(A, full_matrices=True)
    S_inv = np.zeros(A.shape).T
    # S_inv[:N-1,:N-1] = np.diag(1/s[:N-1])
    S_inv[:N,:N] = np.diag(1/s[:N])
    x_svd = ((VT.T.dot(S_inv).dot(U.T))@L) #.reshape(N, shp[0], shp[1])
    La = (A@x_svd).reshape(-1, num)/dT
    La = La.reshape(num, shp[0], shp[1])
    La[:,nodata] = -32767
    return La


def velocity_correction_band(csv_fn, tif_dir, dates, band=1, max_gap=3):
    '''
    Apply velocity corrections on a single band
    '''
    data = pd.read_csv(csv_fn, header=0)
    num = len(data)
    sg = max_gap*(max_gap+1)/2
    assert (num + sg)%max_gap==0
    N = int(((num + sg)//max_gap) - 1)
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


def velocity_correction(csv_fn, tif_dir, dates, max_gap=3):
    '''
    Performs SVD for least square estimation of velocity correction of 
    offset tracking results

    csv_file: Path to csv file 'image_pair.csv'
    tif_dir: Path to offset_tracking tif file results
    dates: list of acquisition interval
    '''
    velocity_correction_band(csv_fn, tif_dir, dates, band=1, max_gap=max_gap)
    velocity_correction_band(csv_fn, tif_dir, dates, band=2, max_gap=max_gap)
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


def stack_offset_tracking(config_path, urls, data_path, polarization, mask, save_path, mode):
    # Reading config files
    with open(config_path, 'r') as f:
        config = json.load(f)
    with open(f'{DIR_PATH}/configs/isce_config.json', 'r') as f:
        config_isce = json.load(f)
    bbox = config_isce["ROI"][1:-1].replace(",","")

    # urls = []
    exec_time = []

    # if os.path.exists(data_path):
    #     shutil.rmtree(data_path)
    os.makedirs(data_path, exist_ok=True)
    cwd = os.getcwd()
    num_images = len(urls)
    
    for url in urls:
        print(data_path+url+'.zip')
        if not os.path.exists(data_path+url+'.zip'):
            download_data(config_isce['ASF_user'], config_isce['ASF_password'], url, data_path, config['Orbit_dir'])
    
    # Downloading orbit data
    # os.system(f'eof --search-path {data_path} --save-dir {config["Orbit_dir"]} --asf-user {config_isce["ASF_user"]} --asf-password {config_isce["ASF_password"]}')

    # Removing the extra SAFE files
    safe_files = [url+'.zip' for url in urls]
    for fn in os.listdir(data_path):
        if fn not in safe_files:
            os.remove(os.path.join(data_path, fn))
    assert len(urls)==len(os.listdir(data_path))
    print("Dataset Download Completed")
    
    os.makedirs(save_path, exist_ok=True)
    os.chdir(save_path)

    print("Creating SAR pairs CSV sheet")
    if not os.path.exists('image_pairs.csv'):
        acquisitionDates = [url.split('_')[5][:8] for url in urls]
        N = len(acquisitionDates)
        elements = [i for i in range(N)]
        id, master, slave, steps = generate_data(elements, acquisitionDates, VELOCITY_CORR_GAP)
        df = pd.DataFrame({'Id':id, 'Master':master, 'Slave':slave, 'Steps':steps, 'Status': [0]*len(id), 
                           'ROI':[str(f"[{bbox.replace(' ', ', ')}]")]*len(id)})
        df.to_csv('image_pairs.csv')

    if mode == 'download':
        print("Completed. Exiting ......")
        return

    start_time = time.time()
    # If coregistration completed then move on to next step
    if len(glob.glob('merged/SLC/*/*.slc.full'))==0:
        print(config_isce["ROI"][1:-1].replace(',',''))
        # execute('cp -r /DATA/glacier-vel/geogrid_req/dem/demLat_N31_N34_Lon_E076_E079* ./')
        
        ## HARDCODED !!!!!!!!!!!!!!!!!!!!!!
        execute('cp -r /DATA/S2_Data/Dem/demLat_N31_N34_Lon_E076_E080.dem* ./')
        dem = glob.glob('*.wgs84')
        print("Using DEM file", dem)
        if not os.path.exists('./run_files'):
            execute(f'python3 {DIR_PATH}/stack/topsStack/stackSentinel.py -s {data_path} -d {dem[0]} -o {config["Orbit_dir"]} -a {config["aux_dir"]} -p {polarization} -e 0.6 --num_proc 64 -b "{config_isce["ROI"][1:-1].replace(",","")}" -t "python3 {DIR_PATH}/stack/topsStack/" -W slc')
            
        run_files = glob.glob('run_files/*')
        log_folder = 'log_files'
        os.makedirs(log_folder, exist_ok=True)
        if os.path.exists(f'{log_folder}/time_taken.txt'):
            os.remove(f'{log_folder}/time_taken.txt')
            
        flag = 0
        strt_time = time.time()
        for file in sorted(run_files):
            txtfile = open(f'{log_folder}/time_taken.txt', 'a')
            txtfile.write(f'Started process {file}\n')
            txtfile.close()
            
            if flag:
                print("Coregisteration complete")
                break
            
            if 'merge' in file:
                txtfile = open(f'{log_folder}/time_taken.txt', 'a')
                txtfile.write(f'Exiting at step {file}\n')
                txtfile.close()
                flag = 1 

            print("Starting Command:", f'bash ./{file}')
            # execute_write(f'bash ./{file}',f'{log_folder}/{file.split("/")[-1]}.txt')
            try:
                num = int(file.split('/')[-1].split('_')[1])
                execute_write(f'bash ./{file}', f'{log_folder}/{file.split("/")[-1]}.txt', f'{log_folder}/{file.split("/")[-1]}_error.txt')
                while process_check_running(num):
                    print("Waiting for the completion of Step:", num)
                    time.sleep(60)
                    
            except Exception as e:
                print('Error occured at', file)
                print(e)
                txtfile = open(f'{log_folder}/exception.txt', 'a')
                txtfile.write(f'Error occured in {file.split("/")[-1]}\n{e}\n')
                txtfile.close()
            
            time_log = f"Time taken till {file}: {time.time()-strt_time}"
            print(time_log)
            txtfile = open(f'{log_folder}/time_taken.txt', 'a')
            txtfile.write(f'{time_log}\n')
            txtfile.close()
            
        print("Time taken for stack coregisteration: ", time.time()-strt_time)

    prev = time.time()
    exec_time.append(prev-start_time)
    write_time('runtime.txt', exec_time)
    start_time = prev
    
    merged_files = glob.glob('merged/SLC/*/*.slc.full')
    if len(merged_files)==num_images:
        print("Complete. Exiting ......")
    else:
        print("Some issue occured. Exiting......")
        return
    
    if mode == 'coreg':
        print("Complete. Exiting ...... coreg")
        os.chdir(cwd)
        return

    # Preparing for geogrid
    dem = glob.glob('*.wgs84')
    dem_vrt = glob.glob('*.dem.vrt')
    generate_dem_products(dem[0], config_isce["ROI"], config)
    
    dates = os.listdir('merged/SLC/')
    os.makedirs('./offset_tracking', exist_ok=True)
    os.chdir("./offset_tracking")

    pair_fn = '../image_pairs.csv'
    secondarys = sorted(os.listdir('../secondarys'))
    if not os.path.exists('window_location.tif'):
        execute(f'python {DIR_PATH}/geogrid_autorift/testGeogrid_ISCE.py -m ../reference -s ../secondarys/{secondarys[0]} -d ../{dem_vrt[0]} -r "{config_isce["ROI"]}" --stackfl ../configs/config_reference')   #  -ssm {config["ssm"]}
    else:
        print("./window_location.tif  ---> Already exists")
    
    print("Starting Offset Tracking")
    offset_compute(pair_fn, config, cwd, mask)
    os.chdir("../")
    
    if mode == 'offset':
        print("Complete. Exiting ...... offset")
        os.chdir(cwd)
        return 

    prev = time.time()
    exec_time.append(prev-start_time)
    write_time('runtime.txt', exec_time)
    start_time = prev
    
    print('Starting Velocity Postprocessing ...')
    # print(os.getcwd())
    # velocity_postprocess(pair_fn, args.shpfile)
    
    if mode == 'post':
        print("Complete. Exiting ...... post")
        os.chdir(cwd)
        return
    
    # print('Starting Velocity Correction ...')
    # os.chdir('../')
    # if os.path.exists('./velocity_corrected'):
    #     shutil.rmtree('./velocity_corrected')
    # os.makedirs('./velocity_corrected')
    # os.chdir('./velocity_corrected')
    # velocity_correction(pair_fn, '../offset_tracking/', dates, VELOCITY_CORR_GAP)
    # os.chdir('../')

    prev = time.time()
    exec_time.append(prev-start_time)
    write_time('runtime.txt', exec_time)
    start_time = prev
    os.chdir(cwd)

    # execute('rm -rf ./geom_reference ./merged/geom_reference')


if __name__=='__main__':
    args = cmdLineParse()
    print(args.config, args.download_txt, args.mode, args.aux)
    urls = []
    with open(args.download_txt, 'r') as f:
        for line in f:
            urls.append(line[:-1])
    
    stack_offset_tracking(args.config, urls, args.data_path, args.polarization, args.mask, args.save_path, args.mode)

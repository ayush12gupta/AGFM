import os, glob, time
import json
import pandas as pd
import numpy as np
import subprocess
from preprocessing.setup_env import download_DEM
from utils import execute, generate_dem_products, process_check_running, cropRaster
from datetime import datetime


VELOCITY_CORR_GAP = 1 ## HARDCODED the number of redundant observations
STACK_MAX_STEP = 3 # Max step size for NCC stacking [default: 3 (36 days)] 
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
    # parser.add_argument('-shp', '--glacier_shp', dest='shpfile', type=str, default=None, help='Glacier Shapefile path')
    parser.add_argument('-msk', '--mask', dest='mask', type=str, default=None, help='Glacier Shapefile path')
    parser.add_argument('-p', '--polarization', type=str, default="vv", help="Polarisation to be used")
    # parser.add_argument('--dem', dest="dem_path", type=str, required=True, help="Path for DEM file")
    return parser.parse_args()


def to_datetime(date):
    date = f'{date[:4]}/{date[4:6]}/{date[-2:]}'
    date = datetime.strptime(date, '%Y/%m/%d')
    return date

def execute_write(cmd, file, file_err):
    # f = open(file, "w")
    f_err = open(file_err, "w") #, stdout=f
    subprocess.check_call(cmd, shell=True, stderr=f_err)


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


def generate_data(elements, acquisitionDates, max_gap=3, max_step=3):
    """
    To create list of image pairs even for redundant cases
    """
    id = []
    master = []
    slave = []
    steps = []
    
    for gap in range(1, max_gap+1):
        idn, mastern, slaven, stepsn = add_elements(elements, acquisitionDates, gap=gap, max_step=max_step)
        id.extend(idn)
        master.extend(mastern)
        slave.extend(slaven)
        steps.extend(stepsn)

    return id, master, slave, steps
    

def add_elements(elements, acquisitionDates, gap, max_step=3):
    """
    To create image pairs considering NCC stacking

    params:
        gap: No. of images between both timesteps in image pair
        max_step:

    Currently hardcoded for max_delta = 36
    """
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
        delt = (to_datetime(acquisitionDates[i+gap]) - to_datetime(acquisitionDates[i])).days
        
        if max_step > 1:
            # Computing deltaT between i and i + max_step
            dt1 = to_datetime(acquisitionDates[i])
            dt2 = to_datetime(acquisitionDates[min(i+max_step, N-1)])
            delt = (dt2 - dt1).days
            
            if delt==(max_step*12):
                step0.extend(acquisitionDates[i:i+max_step+1])
                
            elif delt>(max_step*12):
                deltb = delt
                j = i+max_step+1
                while (deltb>(max_step*12))|((j-i)*12!=deltb):
                    j-=1
                    dt1b = to_datetime(acquisitionDates[i])
                    dt2b = to_datetime(acquisitionDates[min(j, N-1)])
                    deltb = (dt2b - dt1b).days
                
                step0.extend(acquisitionDates[i:j+1])
    
            elif (i+max_step)>N-1:
                j = N-1
                deltb = delt
                while ((j-i)*12!=deltb):
                    j-=1
                    dt1b = to_datetime(acquisitionDates[i])
                    dt2b = to_datetime(acquisitionDates[min(j, N-1)])
                    deltb = (dt2b - dt1b).days
    
                step0.extend(acquisitionDates[i:j+1])
            
        step0 = sorted(list(set(step0)))
        if (gap==max_step)&(delt==(max_step*12)):
            steps.append(','.join(step0[::max_step]))
        else:
            steps.append(','.join(step0))
    return id, master, slave, steps


def offset_tracking(config, cwd, master, slave, steps, mask, deltaT, step_dT=1):

    input_list = ','.join([f'../../merged/SLC/{st}/{st}.slc.full' for st in steps.split(',')])
    execute('cp ../testGeogrid.txt testGeogrid.txt')
    if os.path.exists(f'../../merged/SLC/{master}/{master}.slc.full')&os.path.exists(f'../../merged/SLC/{slave}/{slave}.slc.full'):
        cmd = f'time python {DIR_PATH}/offset_tracking/testautoRIFT_ISCE.py -i {input_list} -g ../window_location.tif -vx ../window_rdr_off2vel_x_vec.tif -vy ../window_rdr_off2vel_y_vec.tif -mpflag {str(config["num_threads"])} --mask {mask} -chipmin {str(config["chip_min"])} -chipmax {str(config["chip_max"])} -dT {int(deltaT)} -step_dT {int(step_dT)}'
        print('PO Command', cmd)
        cmd_file = open(f'offset_cmd.txt', 'a')
        cmd_file.write(cmd)
        cmd_file.close()
        execute(cmd)
    else:
        print("Coregistered images not available check isce.log for issues with ISCE coregisteration")
    
    print("Offset Tracking completed")


def offset_compute(csv_file, config, cwd, shpfile):

    execute('rm -rf ../coreg_secondarys ../ESD ../misreg ../coreg_secondarys ../coarse_interferograms ../geom_reference')
    data_pairs = pd.read_csv(csv_file, header=0)
    mask_img = cropRaster(shpfile, './window_location.tif', './mask.tif')
    
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
            offset_tracking(config, cwd, str(master), str(slave), str(steps), '../mask.tif', deltaT, int(id2)-int(id1))
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


def stack_offset_tracking(config_path, urls, data_path, polarization, shapefile, save_path, mode):
    # Reading config files
    with open(config_path, 'r') as f:
        config = json.load(f)
    with open(config["cred_config"], 'r') as f:
        config_cred = json.load(f)
    
    bbox = config["ROI"][1:-1].replace(",","")

    exec_time = []
    os.makedirs(data_path, exist_ok=True)
    cwd = os.getcwd()
    num_images = len(urls)
    
    merged_files = glob.glob(f'{save_path}/merged/SLC/*/*.slc.full')
    if len(merged_files)!=num_images:
        for url in urls:
            print(data_path+url+'.zip')
            if not os.path.exists(data_path+url+'.zip'):
                download_data(config_cred['ASF_user'], config_cred['ASF_password'], url, data_path, config['Orbit_dir'])
        
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
        id, master, slave, steps = generate_data(elements, acquisitionDates, VELOCITY_CORR_GAP, STACK_MAX_STEP)
        df = pd.DataFrame({'Id':id, 'Master':master, 'Slave':slave, 'Steps':steps, 'Status': [0]*len(id), 
                           'ROI':[str(f"[{bbox.replace(' ', ', ')}]")]*len(id)})
        df.to_csv('image_pairs.csv')

    if mode == 'download':
        print("Completed. Exiting ......")
        return

    start_time = time.time()
    # If coregistration completed then move on to next step
    if len(glob.glob('merged/SLC/*/*.slc.full'))==0:
        print(config["ROI"][1:-1].replace(',',''))
        # execute('cp -r /DATA/glacier-vel/geogrid_req/dem/demLat_N31_N34_Lon_E076_E079* ./')
        
        ## HARDCODED !!!!!!!!!!!!!!!!!!!!!!
        # execute('cp -r /DATA/S2_Data/Dem/demLat_N31_N34_Lon_E076_E080.dem* ./')
        # dem = glob.glob('*.wgs84')

        # Downloading DEM file and converting to ISCE file
        roi = np.array(config["ROI"][1:-1].replace(' ','').split(',')).astype('float')
        download_DEM(roi, out_path='./dem_roi.tif')
        dem = glob.glob('*.dem')
        print("Using DEM file", dem)

        if not os.path.exists('./run_files'):
            execute(f'python3 {DIR_PATH}/stack/topsStack/stackSentinel.py -s {data_path} -d {dem[0]} -o {config["Orbit_dir"]} -a {config["aux_dir"]} -p {polarization} -e 0.6 --num_proc 64 -b "{config["ROI"][1:-1].replace(",","")}" -t "python3 {DIR_PATH}/stack/topsStack/" -W slc')
            
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
    # dem = glob.glob('*.wgs84')
    dem = glob.glob('*.dem')
    dem_vrt = glob.glob('*.dem.vrt')
    generate_dem_products(dem[0], config["ROI"], config)
    
    dates = os.listdir('merged/SLC/')
    os.makedirs('./offset_tracking', exist_ok=True)
    os.chdir("./offset_tracking")

    pair_fn = '../image_pairs.csv'
    secondarys = sorted(os.listdir('../secondarys'))
    if not os.path.exists('window_location.tif'):
        execute(f'python {DIR_PATH}/offset_tracking/testGeogrid_ISCE.py -m ../reference -s ../secondarys/{secondarys[0]} -d ../{dem_vrt[0]} -r "{config["ROI"]}" --stackfl ../configs/config_reference')   #  -ssm {config["ssm"]}
    else:
        print("./window_location.tif  ---> Already exists")
    
    print("Starting Offset Tracking")
    offset_compute(pair_fn, config, cwd, shapefile)
    os.chdir("../")
    
    if mode == 'offset':
        print("Complete. Exiting ...... offset")
        os.chdir(cwd)
        return 

    prev = time.time()
    exec_time.append(prev-start_time)
    write_time('runtime.txt', exec_time)
    start_time = prev
    
    # print('Starting Velocity Postprocessing ...')
    # # print(os.getcwd())
    # # velocity_postprocess(pair_fn, args.shpfile)
    
    # if mode == 'post':
    #     print("Complete. Exiting ...... post")
    #     os.chdir(cwd)
    #     return
    
    # prev = time.time()
    # exec_time.append(prev-start_time)
    # write_time('runtime.txt', exec_time)
    # start_time = prev
    # os.chdir(cwd)


if __name__=='__main__':
    args = cmdLineParse()
    print(args.config, args.download_txt, args.mode, args.aux)
    urls = []
    with open(args.download_txt, 'r') as f:
        for line in f:
            urls.append(line[:-1])
    
    stack_offset_tracking(args.config, urls, args.data_path, args.polarization, args.mask, args.save_path, args.mode)

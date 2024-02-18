import os, shutil
import glob
import time
import json
from datetime import datetime
from stack_process import stack_offset_tracking
from batch3Dinversion import batch_inversion


def cmdLineParse():
    '''
    Command line parser.
    '''
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-t_asc', '--download_asc_txt', type=str, required=True, help="Data Ascending CSV file")
    parser.add_argument('-t_des', '--download_des_txt', type=str, required=True, help="Data Descending CSV file")
    parser.add_argument('--config', type=str, default="./configs/pipeline_config.json", help="Data config file")
    # parser.add_argument('--dem', dest="dem_path", type=str, required=True, help="Path for DEM file")
    return parser.parse_args()


def file_pair_selection(download_asc_txt, download_des_txt):
    file_asc = open(download_asc_txt, "r")
    data_asc = file_asc.read().split('\n')[:-1]
    date_asc = [fl.split('_')[5][:8] for fl in data_asc]
    date_asc = [f'{fl[:4]}/{fl[4:6]}/{fl[6:8]}' for fl in date_asc]
    date_asc = [datetime.strptime(dt, '%Y/%m/%d').date() for dt in date_asc]
    
    file_des = open(download_des_txt, "r")
    data_des = file_des.read().split('\n')[:-1]
    date_des = [fl.split('_')[5][:8] for fl in data_des]
    date_des = [f'{fl[:4]}/{fl[4:6]}/{fl[6:8]}' for fl in date_des]
    date_des = [datetime.strptime(dt, '%Y/%m/%d').date() for dt in date_des]
    
    des_fl, asc_fl = [], []
    des_idx, asc_idx = 0, 0
    while (des_idx<len(data_des))&(asc_idx<len(data_asc)):
        diff = (date_des[des_idx] - date_asc[asc_idx]).days
        if diff<0:
            if abs(diff)>1:
                des_idx+=1
            else:
                des_fl.append(data_des[des_idx])
                asc_fl.append(data_asc[asc_idx])
                des_idx+=1
                asc_idx+=1
        elif diff>0:
            if abs(diff)>1:
                asc_idx+=1
            else:
                des_fl.append(data_des[des_idx])
                asc_fl.append(data_asc[asc_idx])
                des_idx+=1
                asc_idx+=1
                
    return asc_fl, des_fl


def main(inps):
    with open(inps.config, 'r') as f:
        config = json.load(f)
        
    os.makedirs(config['save_path'], exist_ok=True)
    os.chdir(config['save_path'])
    
    asc_txt, des_txt = file_pair_selection(inps.download_asc_txt, inps.download_des_txt)
        
    # # Performing coregistration + offset tracking for ascending images
    print("Current working directory:", os.getcwd())
    stack_offset_tracking(config['config_path'], asc_txt, f'{config["SAR_dir"]}/ascending/', config['polarisation'].lower(), \
                          config['mask'], f'{config["save_path"]}/stack_asc/', 'post')
    
    # Performing coregistration + offset tracking for descending images
    print("Current working directory:", os.getcwd())
    stack_offset_tracking(config['config_path'], des_txt, f'{config["SAR_dir"]}/descending/', config['polarisation'].lower(), \
                          config['mask'], f'{config["save_path"]}/stack_des/', 'post')

    # 3D Inversion
    print("Current working directory:", os.getcwd())
    batch_inversion(f'{config["save_path"]}/stack_asc/', f'{config["save_path"]}/stack_des/', f'{config["save_path"]}/3D_Velocities/')



if __name__=='__main__':
    inps = cmdLineParse()
    main(inps)

from turtle import down
import os
import json
import argparse

from get_orbit import get_orbit_fl


parser = argparse.ArgumentParser()
parser.add_argument('--reference', type=str, required=True, help="URL of reference zip file")
parser.add_argument('--secondary', type=str, required=True, help="URL of secondary zip file")
parser.add_argument('--orbit_path', type=str, default="./orbits/", help="directory in which orbit file needs to be saved")
parser.add_argument('--data_pathR', type=str, default="./data/", help="directory in which reference data files needs to be saved")
parser.add_argument('--data_pathS', type=str, default="./data/", help="directory in which secondary data files needs to be saved")
parser.add_argument('--config', type=str, default="./config/isce_config.json", help="ISCE config file")

args = parser.parse_args()


def download_data(args, username, password, url, data_path):
    dmy= url.split('/')[-1].split('_')[5]
    sentinel_type = url.split('/')[-1].split('_')[0]
    date, month, yr = dmy[6:8], dmy[4:6], dmy[:4]
    date = "%02d" % (int(date))
    rt = os.system("wget " + url + " -P " + data_path + " --user=" + username + " --password=" + password + " -nc")
    if rt!=0:
        print("Failed downloading ", url)
    get_orbit_fl(sentinel_type, date, month, yr, args.orbit_path)
    return url.split('/')[-1]


def setup(args):
    """
    Setting up the topsApp.xml ISCE input file and downloading SAR SLC images and their orbit file
    """

    # Reading config file
    with open("isce_config.json", 'r') as f:
        config = json.load(f)

    if not os.path.exists(args.data_pathR):
        os.mkdir(args.data_pathR)
    if not os.path.exists(args.data_pathS):
        os.mkdir(args.data_pathS)
    if not os.path.exists(args.orbit_path):
        os.mkdir(args.orbit_path)
    # Downloading Reference and secondary images
    ref_fn = download_data(args, args.reference, args.data_pathR)
    sec_fn = download_data(args, args.secondary, args.data_pathS)

    # Setup secondary.xml and reference.xml
    reference_xml = '''<component name="reference">
    <property name="orbit directory">{0}</property>
    <property name="output directory">reference</property>
    <property name="polarization">vv</property>
    <property name="safe">['{1}']</property>
</component>'''.format(args.orbit_path, args.data_pathR + ref_fn)
    with open("reference.xml", "w") as fid:
            fid.write(reference_xml)
    
    secondary_xml = '''<component name="secondary">
    <property name="orbit directory">{0}</property>
    <property name="output directory">secondary</property>
    <property name="polarization">vv</property>
    <property name="safe">['{1}']</property>
</component>'''.format(args.orbit_path, args.data_pathS + sec_fn)
    with open("secondary.xml", "w") as fid:
            fid.write(secondary_xml)
    
    # Setup topsApp.xml
    tops_xml = '''<?xml version="1.0" encoding="UTF-8"?>
    <topsApp>
    <component name="topsinsar">
        <property name="Sensor name">SENTINEL1</property>
        <component name="reference">
            <catalog>reference.xml</catalog>
        </component>
        <component name="secondary">
            <catalog>secondary.xml</catalog>
        </component>
        <property name="swaths">{0}</property>
        <property name="region of interest">{1}</property>
        <property name="do ESD">{2}</property>
        <property name="do interferogram">{3}</property>
        <property name="ESD coherence threshold">0.75</property>	    
        <property name="ampcor skip width">{4}</property>
        <property name="ampcor skip height">{4}</property>
        <property name="ampcor search window width">{5}</property>
        <property name="ampcor search window height">{5}</property>
        <property name="ampcor window width">{6}</property>
        <property name="ampcor window height">{6}</property>
        <property name="do ionosphere correction">{7}</property>
        <property name="apply ionosphere correction">{7}</property>
        <property name="consider burst properties in ionosphere computation">{8}</property>
        <property name="start ionosphere step">{9}</property>
        <property name="end ionosphere step">{10}</property>
        <property name="height of ionosphere layer in km">{11}</property>
        <property name="apply polynomial fit before filtering ionosphere phase">{12}</property>
        <property name="maximum window size for filtering ionosphere phase">{13}</property>
        <property name="minimum window size for filtering ionosphere phase">{14}</property>
        <property name="maximum window size for filtering ionosphere azimuth shift">{15}</property>
        <property name="minimum window size for filtering ionosphere azimuth shift">{16}</property>
        <property name="demfilename">demLat_N31_N34_Lon_E076_E079.dem.wgs84</property>
        <!--<property name="geocode demfilename">path_to_your_dem</property>-->
        <!--property name="geocode list">['merged/phsig.cor', 'merged/filt_topophase.unw', 'merged/los.rdr', 'merged/topophase.flat', 'merged/filt_topophase.flat','merged/topophase.cor','merged/filt_topophase.unw.conncomp']</property>-->
    </component>
    </topsApp>'''.format(config['swath'],config['ROI'],config['do_ESD'],config['do_interferogram'],config['ampcor_skip'],config['ampcor_search'],config['ampcor_window'],
                        config['do_ionosphere_correct'],config['ionosphere']['consider_burst'],config['ionosphere']['start'],config['ionosphere']['end'],config['ionosphere']['height'],
                        config['ionosphere']['apply_polyn_fit'],config['ionosphere']['max_window_phase'],config['ionosphere']['min_window_phase'],config['ionosphere']['max_window_azi'],config['ionosphere']['min_window_azi'])
    with open("topsApp.xml", "w") as fid:
        fid.write(tops_xml)


if __name__=='__main__':
    setup(args)
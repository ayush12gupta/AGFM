from urllib.request import Request, urlopen, urlretrieve
from bs4 import BeautifulSoup
import os
import datetime


def prev_day(day, month, yr):
    date = datetime.date(int(yr), int(month), int(day))
    date = date - datetime.timedelta(days=1)
    return date.day, date.month, date.year


def read_url(url, date):
    url = url.replace(" ","%20")
    req = Request(url)
    try:
        a = urlopen(req).read()
    except:
        return 0
    soup = BeautifulSoup(a, 'html.parser')
    x = (soup.find_all('a'))
    for i in x:
        file_name = i.extract().get_text()
        url_new = url + file_name
        url_new = url_new.replace(" ","%20")
        if(file_name[-1]=='/' and file_name[0]!='.'):
            read_url(url_new)
        if(url_new[-3:]=='zip'):
            dt = url_new.split('/')[-1].split('_')[-2][1:]
            if dt[6:8]==date:
                return url_new


def get_orbit_fl(sentinel_type, date, month, yr, orbit_direct):
    # base_url = "http://step.esa.int/auxdata/orbits/Sentinel-1/POEORB"
    datep, monthp, yrp = prev_day(date, month, yr)
    datep = "%02d" % (int(datep)); date = "%02d" % (int(date))
    monthp = "%02d" % (int(monthp)); month = "%02d" % (int(month))
    yrp = "%02d" % (int(yrp)); yr = "%02d" % (int(yr))
    print(f'sentinelsat --gnss -s {yrp}{monthp}{datep} -e {yr}{month}{date} --producttype AUX_POEORB --query="platformserialidentifier={sentinel_type[1:]}" -d --path {orbit_direct}')
    os.system(f'sentinelsat --gnss -s {yrp}{monthp}{datep} -e {yr}{month}{date} --producttype AUX_POEORB --query="platformserialidentifier={sentinel_type[1:]}" -d --path {orbit_direct}')
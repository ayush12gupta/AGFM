
def create_configISCE(config, roi, dem_path, type):

    if type=='topsApp':
        create_topsApp(config, roi, dem_path, type)
    elif type=='demApp':
        create_demApp(config, roi, dem_path, type)


def create_topsApp(config, roi, dem_path, name='topsApp'):
    # Setup topsApp.xml
    xml = '''<?xml version="1.0" encoding="UTF-8"?>
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
        <property name="demfilename">{17}</property>
        <!--<property name="geocode demfilename">path_to_your_dem</property>-->
        <!--property name="geocode list">['merged/phsig.cor', 'merged/filt_topophase.unw', 'merged/los.rdr', 'merged/topophase.flat', 'merged/filt_topophase.flat','merged/topophase.cor','merged/filt_topophase.unw.conncomp']</property>-->
    </component>
    </topsApp>'''.format(config['swath'], roi, config['do_ESD'],config['do_interferogram'],config['ampcor_skip'],config['ampcor_search'],config['ampcor_window'],
                        config['do_ionosphere_correct'],config['ionosphere']['consider_burst'],config['ionosphere']['start'],config['ionosphere']['end'],config['ionosphere']['height'],
                        config['ionosphere']['apply_polyn_fit'],config['ionosphere']['max_window_phase'],config['ionosphere']['min_window_phase'],config['ionosphere']['max_window_azi'],config['ionosphere']['min_window_azi'], dem_path)
    with open(name + ".xml", "w") as fid:
        fid.write(xml)


def create_demApp(config, roi, dem_path, name='demApp'):
    # Setup topsApp.xml
    xml = '''<?xml version="1.0" encoding="UTF-8"?>
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
        <property name="range looks">19</property>
        <property name="azimuth looks">7</property>
        <property name="do ESD">{2}</property>
        <property name="do interferogram">True</property>
        <property name="do unwrap">True</property>
    	<property name="unwrapper name">snaphu_mcf</property>
        <property name="ESD coherence threshold">0.75</property>	    
        <property name="ampcor skip width">{3}</property>
        <property name="ampcor skip height">{3}</property>
        <property name="ampcor search window width">{4}</property>
        <property name="ampcor search window height">{4}</property>
        <property name="ampcor window width">{5}</property>
        <property name="ampcor window height">{5}</property>
        <property name="do ionosphere correction">{6}</property>
        <property name="apply ionosphere correction">{6}</property>
        <property name="consider burst properties in ionosphere computation">{7}</property>
        <property name="start ionosphere step">{8}</property>
        <property name="end ionosphere step">{9}</property>
        <property name="height of ionosphere layer in km">{10}</property>
        <property name="apply polynomial fit before filtering ionosphere phase">{11}</property>
        <property name="maximum window size for filtering ionosphere phase">{12}</property>
        <property name="minimum window size for filtering ionosphere phase">{13}</property>
        <property name="maximum window size for filtering ionosphere azimuth shift">{14}</property>
        <property name="minimum window size for filtering ionosphere azimuth shift">{15}</property>
        <property name="demfilename">{16}</property>
        <!--<property name="geocode demfilename">path_to_your_dem</property>-->
        <property name="geocode list">['merged/phsig.cor', 'merged/filt_topophase.unw', 'merged/range.rdr', 'merged/look_angle.rdr', 'merged/bp.rdr', 'merged/bn.rdr', 'merged/los.rdr', 'merged/z.rdr', 'merged/incid.rdr', 'merged/topophase.flat', 'merged/filt_topophase.flat','merged/topophase.cor', 'merged/elevation.rdr','merged/filt_topophase.unw.conncomp']</property>
    </component>
    </topsApp>'''.format(config['swath'], roi, config['do_ESD'],config['ampcor_skip'],config['ampcor_search'],config['ampcor_window'],
                        config['do_ionosphere_correct'],config['ionosphere']['consider_burst'],config['ionosphere']['start'],config['ionosphere']['end'],config['ionosphere']['height'],
                        config['ionosphere']['apply_polyn_fit'],config['ionosphere']['max_window_phase'],config['ionosphere']['min_window_phase'],config['ionosphere']['max_window_azi'],config['ionosphere']['min_window_azi'], dem_path)
    with open(name + ".xml", "w") as fid:
        fid.write(xml)
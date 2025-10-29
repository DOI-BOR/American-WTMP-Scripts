import os, sys
from com.rma.model import Project
from hec.heclib.util import HecTime

# print current path
print("Current paths: ", sys.path)

# create list of unwanted folders in sys.path
search_list = ["SacTrn", "Sacramento", "American"]

# initialize and search for unwanted paths
matching_paths = []
for p in sys.path:
    if any(phrase in p for phrase in search_list):
        matching_paths.append(p)

# print paths containing unwanted phrases
print("Paths to be removed:")
for path in matching_paths:
    print(path)

# remove matching paths from sys.path
for path in matching_paths:
    if path in sys.path:
        sys.path.remove(path)

# append path
sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))

import DSS_Tools
reload(DSS_Tools)

import DMS_preprocess
reload(DMS_preprocess)

import Forecast_preprocess as fpp
reload(fpp)

import create_balance_flow_jython as cbfj
reload(cbfj)



def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')

    rtw = computeOptions.getRunTimeWindow()
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()
    startyear_str = starttime_str[5:9]
    year = int(startyear_str)

    run_dir = computeOptions.getRunDirectory()
    model_name = 'ResSim'  # base Folsom forecast model
    if "nobypass" in run_dir.lower():
        model_name += " NoBypass"  
    shared_dir = os.path.join(fpp.study_dir_from_run_dir(run_dir),'shared')
    forecast_dss = os.path.join(shared_dir,'WTMP_American_forecast.dss')

    # forecast met data potentially has some weird units
    DMS_preprocess.fix_DMS_types_units(forecast_dss)

    # predict daily elevation, in case of start on arbitrary day
    fpp.write_forecast_elevations(currentAlternative, rtw, forecast_dss, shared_dir)

    # flow balance: subtract muni/pump flow from total release to get flow through dam
    fpp.subtract_muni_pump(forecast_dss)

    # get target temp location from dummy dss file somehow?

    # load Tair and Folsom FaveFlow
    start_doy = HecTime(starttime_str).dayOfYear() #.value()
    endtime_doy = HecTime(endtime_str).dayOfYear() #.value()

    # I don't think we need this split out; everything is dynamic including river outlets?
    folsom_outlets = fpp.write_qot_7outlets_flows(forecast_dss, starttime_str, endtime_str)

    # make equilibrium temp for ResSim target temp calcs
    currentAlternative.addComputeMessage("Computing equilibrium temperature, this may take a while...")
    # eq_temp(rtw,at,cl,ws,sr,td,eq_temp_out)
    fpp.eq_temp(rtw,
            [forecast_dss,"/MR AM.-NATOMA LAKE/FAIR OAKS/TEMP-AIR//1HOUR/251.40.53.1.1/"],
            [forecast_dss,"/MR AM.-LOWER AMERICAN R./MATHER AFB/%-CLOUD COVER-FRAC//1Hour/252.50.129.1.1/"],
            [forecast_dss,"/MR AM.-NATOMA LAKE/FAIR OAKS/SPEED-WIND//1HOUR/251.40.133.1.1/"],
            [forecast_dss,"/MR AM.-FOLSOM LAKE/THEORETICAL FOLSOM SOLAR/IRRAD-SOLAR//1HOUR/250.15.135.1.1/"],
            [forecast_dss,"/MR AM.-NATOMA LAKE/FAIR OAKS/TEMP-DEWPOINT//1HOUR/251.40.51.1.1/"],
            [forecast_dss,"/MR AM.-NATOMA LAKE/FAIR OAKS/Temp-Equil//1Hour/amer_bc_script/"]
           )

    # write an hourly forecast elevation based on starting elevation and flows
    DSS_Tools.resample_dss_ts(forecast_dss,'//FOLSOM/PUMPING (FP)//1Day/AMER_BC_SCRIPT/',rtw,forecast_dss,'1HOUR')
    DSS_Tools.resample_dss_ts(forecast_dss,'//Folsom-NF-in/FLOW-IN//1Day/AMER_BC_SCRIPT/',rtw,forecast_dss,'1HOUR')
    DSS_Tools.resample_dss_ts(forecast_dss,'//Folsom-SF-in/FLOW-IN//1Day/AMER_BC_SCRIPT/',rtw,forecast_dss,'1HOUR')
    DSS_Tools.resample_dss_ts(forecast_dss,'/AMERICAN RIVER/FOLSOM LAKE/FLOW-ACC-DEP//1Day/AMER_BC_SCRIPT/',rtw,forecast_dss,'1HOUR')
    inflow_records = ['//Folsom-NF-in/FLOW-IN//1Hour/AMER_BC_SCRIPT/',
                      '//Folsom-SF-in/FLOW-IN//1Hour/AMER_BC_SCRIPT/',
                      '/AMERICAN RIVER/FOLSOM LAKE/FLOW-ACC-DEP//1Hour/AMER_BC_SCRIPT/']  # this actually evap, but negative already, so it goes as inflow
    outflow_records = ['//FOLSOM/FLOW-RELEASE//1Hour/AMER_BC_SCRIPT/']
    
    # make constant records for linking
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=2.0, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='2cfs',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=2.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='2cfs',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=10.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='TENS',fpart='TENS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='evap', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    print('Making ELEV-INITIAL-GATE records...')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=401.0, what='elev', 
                        dss_type='PER-AVER', period='1DAY',cpart='ELEV-INITIAL-GATE',fpart='FULL-HEIGHT')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=362.0, what='elev', 
                        dss_type='PER-AVER', period='1DAY',cpart='ELEV-INITIAL-GATE',fpart='1-OUT')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=336.0, what='elev', 
                        dss_type='PER-AVER', period='1DAY',cpart='ELEV-INITIAL-GATE',fpart='2-OUT')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=307.0, what='elev', 
                        dss_type='PER-AVER', period='1DAY',cpart='ELEV-INITIAL-GATE',fpart='ALL-OUT')


    # if this is a non-iterative simulation, put the correct ensemble number into the Folsom.in file 
    if 'nobypass' in model_name.lower():
        print('Setting: Do not use Folsom River Bypass in forecast')
        fpp.remove_folsom_lower_river_use(forecast_dss,"/FOLSOM/LOWER RIVER OUTLET USEAGE/COUNT//1Day/AMER_BC_SCRIPT/")

    return True

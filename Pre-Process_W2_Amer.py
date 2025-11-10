import os, sys
from com.rma.model import Project
from hec.heclib.util import HecTime

# print current path
print("Current paths: ", sys.path)

# create list of unwanted folders in sys.path
search_list = ["SacTrn", "Sacramento"]

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
    model_name = 'W2 Folsom'  # base Folsom forecast model (iterative w/ bypass)
    if "FixedATSP" in run_dir:
        model_name += " FixedATSP"
    if "NoBypass" in run_dir:
        model_name += " NoBypass"
    model_dir = fpp.model_dir_from_run_dir(run_dir,'Folsom',model_name)
    targt_temp_npt_filepath = os.path.join(model_dir,'TargetSchedulesA.npt') # overwrite what's there    
    shared_dir = os.path.join(fpp.study_dir_from_run_dir(run_dir),'shared')
    forecast_dss = os.path.join(shared_dir,'WTMP_American_forecast.dss')
    schedule_csv = os.path.join(model_dir,'SchedulesA.csv')

    # forecast met data potentially has some weird units
    DMS_preprocess.fix_DMS_types_units(forecast_dss)

    # convert Folsom storage to monthly elevation
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_Folsom.csv'), 'Folsom')
    fpp.storage_to_elev('Folsom',elev_stor_area,forecast_dss,'//FOLSOM/STORAGE//1Month/AMER_BC_SCRIPT/',conic=False)

    # invent Natoma elevation record from folsom storage rec
    fpp.invent_elevation('Natoma',forecast_dss,'//FOLSOM/STORAGE//1Month/AMER_BC_SCRIPT/',123.0)

    # predict daily elevation, in case of start on arbitrary day
    fpp.write_forecast_elevations(currentAlternative, rtw, forecast_dss, shared_dir)

    # flow balance: subtract muni/pump flow from total release to get flow through dam
    fpp.subtract_muni_pump(forecast_dss)
    
    # flow balance: split evap for W2
    fpp.split_folsom_evap(forecast_dss,'/AMERICAN RIVER/FOLSOM LAKE/FLOW-ACC-DEP//1Day/AMER_BC_SCRIPT/')

    # divide Nimbus dam flow into 3, for the 3 gated spillways in W2 Natoma model
    fpp.split_nimbus_outflow(forecast_dss,'/AMERICAN RIVER/LAKE NATOMA/FLOW-NIMBUS ACTUAL//1Day/AMER_BC_SCRIPT/')

    # get target temp location from dummy dss file somehow?

    # load Tair and Folsom FaveFlow
    start_doy = HecTime(starttime_str).dayOfYear() #.value()
    endtime_doy = HecTime(endtime_str).dayOfYear() #.value()

    doys,FaveFlow,Tair = fpp.load_tt_data(forecast_dss, starttime_str, endtime_str) # day-of-year,CMS,C

    target_temp_write = fpp.write_target_temp_npt(year,1,doys,Tair,FaveFlow,schedule_csv,targt_temp_npt_filepath,lagWatt=True)

    folsom_outlets = fpp.write_qot_7outlets_flows(forecast_dss, starttime_str, endtime_str)

    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='TinyFlow',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='TinyFlow',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=10.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='TENS',fpart='TENS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')

    # update W2 interative control files for restart date
    # this functionality is already in the w2 plugin - did you check to see if folsom iterative compute option is checked??
    fpp.update_W2_Folsom_iterative_restart_date(rtw,model_dir)

    # if this is a non-iterative simulation, put the correct ensemble number into the 'Set Guess' file 
    # -- edit, 2025-01-27, B. Saenz.  We want to set this correct ensemble for all runs, as it tells also
    # tells the iterative version of Folsom W2 to start with a parituclar schedule, so save time
    #if 'FixedATSP' in model_name:

    w2run_base,_ = os.path.split(run_dir)
    fpp.update_W2_Folsom_iterative_schedule_number(w2run_base,model_dir)

    if target_temp_write and folsom_outlets:
        return True

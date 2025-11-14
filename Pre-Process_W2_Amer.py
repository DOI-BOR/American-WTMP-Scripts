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

shutter_links = ['init_elev_shutter_1','init_elev_shutter_2','init_elev_shutter_3']

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

    # get target temp location from dss
    TT_loc = fpp.get_downstream_loc(forecast_dss)

    # load Tair and Folsom FaveFlow
    start_doy = HecTime(starttime_str).dayOfYear() #.value()
    endtime_doy = HecTime(endtime_str).dayOfYear() #.value()

    doys,FaveFlow,Tair = fpp.load_tt_data(forecast_dss, starttime_str, endtime_str) # day-of-year,CMS,C

    target_temp_write = fpp.write_target_temp_npt(year,TT_loc,doys,Tair,FaveFlow,schedule_csv,targt_temp_npt_filepath,lagWatt=True)

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
    print('Making ELEV-INITIAL-GATE records...')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=401.0, what='elev', 
                        dss_type='PER-AVER', period='1DAY',cpart='ELEV-INITIAL-GATE',fpart='FULL-HEIGHT')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=362.0, what='elev', 
                        dss_type='PER-AVER', period='1DAY',cpart='ELEV-INITIAL-GATE',fpart='1-OUT')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=336.0, what='elev', 
                        dss_type='PER-AVER', period='1DAY',cpart='ELEV-INITIAL-GATE',fpart='2-OUT')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=307.0, what='elev', 
                        dss_type='PER-AVER', period='1DAY',cpart='ELEV-INITIAL-GATE',fpart='ALL-OUT')

    # find initial shutter elevations
    locations_obj = currentAlternative.getInputDataLocations()
    shutter_objs = DSS_Tools.organizeLocations(currentAlternative, locations_obj, shutter_links, return_dss_paths=False)
    w2_elevs = fpp.get_initial_shutter_positions(rtw, currentAlternative, computeOptions, shutter_objs)
    print('Initial W2 shutter elevs:',w2_elevs)

    # update W2 interative control files for restart date - wait this functionality is already in the w2 plugin - did you check to see if folsom iterative compute option is checked??
    # also update the inital gate elevations in the w2_con file
    fpp.update_W2_Folsom_iterative_restart_date_and_shutters(rtw,model_dir,w2_elevs)

    # if this is a non-iterative simulation, put the correct ensemble number into the 'Set Guess' file 
    #if 'FixedATSP' in model_name:
    # -- edit, 2025-01-27, B. Saenz.  We want to set the  correct ensemble for all runs, as it also
    # tells the iterative version of Folsom W2 to start with a parituclar schedule, to save time if you have a good guess
    w2run_base,_ = os.path.split(run_dir)
    fpp.update_W2_Folsom_iterative_schedule_number(w2run_base,model_dir)

    # The W2 selective file must be modified depending on the start gate to get the auto shutter to work correctly
    fpp.modify_w2_selective_start_date(rtw,os.path.join(model_dir,'w2_selective.npt'))

    if target_temp_write and folsom_outlets:
        return True

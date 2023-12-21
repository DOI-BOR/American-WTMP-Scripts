import os, sys
from com.rma.model import Project
from hec.heclib.util import HecTime

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
    model_dir = fpp.model_dir_from_run_dir(run_dir,'Folsom','W2 Folsom')
    targt_temp_npt_filepath = os.path.join(model_dir,'TargetSchedulesA.npt') # overwrite what's there    
    shared_dir = os.path.join(fpp.study_dir_from_run_dir(run_dir),'shared')
    forecast_dss = os.path.join(shared_dir,'WTMP_American_forecast.dss')
    schedule_csv = os.path.join(model_dir,'SchedulesA.csv')

    # forecast met data potentially has some weird units
    DMS_preprocess.fix_DMS_types_units(forecast_dss)

	# convert Folsom storage to elevation
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_Folsom.csv'), 'Folsom')
    fpp.storage_to_elev('Folsom',elev_stor_area,forecast_dss,'//FOLSOM/STORAGE//1Month/AMER_BC_SCRIPT/',conic=False)

    # invent Natoma elevation record from folsom storage rec
    fpp.invent_elevation('Natoma',forecast_dss,'//FOLSOM/STORAGE//1Month/AMER_BC_SCRIPT/',123.0)

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
    starting_elevation = DSS_Tools.first_value(forecast_dss,'//Folsom/ELEV//1Month/AMER_BC_SCRIPT/')
    print('starting_elevation ',starting_elevation)
    cbfj.predict_elevation(currentAlternative, rtw, 'Folsom', inflow_records, outflow_records, starting_elevation,
                         elev_stor_area, forecast_dss, '//Folsom/ELEV-FORECAST//1HOUR/AMER_BC_SCRIPT/', forecast_dss, shared_dir,
                         use_conic=False, alt_period=None, alt_period_string=None)

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
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='elev', 
                        dss_type='PER-AVER', period='1DAY',cpart='ELEV-INITIAL-GATE',fpart='FULL-HEIGHT')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='evap', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')

    return True

import sys,os,csv
print(sys.path)
from com.rma.io import DssFileManagerImpl
from hec.heclib.dss import HecDss
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE
import math
import datetime as dt
import hec.hecmath.TimeSeriesMath as tsmath
from hec.heclib.util import HecTime

import flowweightaverage
reload(flowweightaverage)

import DSS_Tools
reload(DSS_Tools)

import Forecast_preprocess as fpp
reload(fpp)

import tz_offset
reload(tz_offset)

##
#
# computeAlternative function is called when the ScriptingAlternative is computed.
# Arguments:
#   currentAlternative - the ScriptingAlternative. hec2.wat.plugin.java.impl.scripting.model.ScriptPluginAlt
#   computeOptions     - the compute options.  hec.wat.model.ComputeOptions
#
# return True if the script was successful, False if not.
# no explicit return will be treated as a successful return
#
##

# location.getName(),getParameter()

W2FolsomFlowTempLocs = [
  ['Qwo-105-str1','Two-105-str1'],
  ['Qwo-105-str2','Two-105-str2'],
  ['Qwo-105-str3','Two-105-str3'],
  ['Qwo-105-str4','Two-105-str4'],
  ['Qwo-105-str6','Two-105-str6'], # five is not used
  ['Qwo-105-str7','Two-105-str7'],
]

W2FolsomFlowLocs = [ft[0] for ft in W2FolsomFlowTempLocs]

OutputLocs = [
    'W2_Natoma_InflowQ',
    'W2_Natoma_InflowT',
    'W2_Folsom_Schedule_Temp',
    'W2_Folsom_Schedule',
    'W2_Folsom_Downstream_Temp'
]



def snap_folsom_elev(elev_in):
    '''Plotting routine are expecting exact values to discriminate shutter positions. Snap to known positions'''
    if elev_in >= 395.:
        return 401.0
    elif elev_in >= 355.:
        return 362.0
    elif elev_in >= 346.:
        return 349.0
    elif elev_in >= 320.:
        return 336.0
    elif elev_in >= 300.:
        return 307.0
    else:
        return elev_in

                
def read_str_csv(csv_file_path):
    ''' Reads W2 Folsom output str CSV file, getting shutter elevations and date

    OK, here is the deal.  In W2 forcast mode when bypass is allowed, if bypass is
    needed, W2 essentially turns penstock 1 into the bypass by lowering the intake
    elevation to 64.01 ft, whle at the same time shunting and re-caldulating powerplant
    flow so that is goes as additional flow through the remaining two powerhouse 
    penstocks. This means that we need to parse penstock 1 elevations, and identify when it
    reaches ecatly 64.01, and add that flow to the bypass flow. Normally, any bypass flow, which
    must be specified, occurs in column 13 (starting with column 0).
    
    '''
    strReader = csv.reader(open(csv_file_path), delimiter=',', quotechar='|')
    jday = []
    sElev1 = []
    sElev2 = []
    sElev3 = []
    bypassFlow = []
    revisedPenstock1Flow = []
    
    for i,row in enumerate(strReader):
        if i>1: # 2 header rows
            jday.append(float(row[0]))
            e1 = float(row[15])

            if e1 == 64.01:
                bypassFlow.append(float(row[8])+float(row[13]))
                revisedPenstock1Flow.append(0.0)
            else:
                bypassFlow.append(float(row[13]))
                revisedPenstock1Flow.append(float(row[8]))
            
            sElev1.append(snap_folsom_elev(e1*3.28084)) # m to ft
            sElev2.append(snap_folsom_elev(float(row[16])*3.28084)) # m to ft
            sElev3.append(snap_folsom_elev(float(row[17])*3.28084)) # m to ft            

    return jday,sElev1,sElev2,sElev3,bypassFlow,revisedPenstock1Flow

def read_storage_csv(csv_file_path):
    ''' Reads W2 Folsom output str CSV file, getting shutter elevations and date
    '''
    m3_to_acft = 0.000810714 # we are going to work in ac-ft since I don't know if the other codes/plotting can do other units
    strReader = csv.reader(open(csv_file_path), delimiter=' ', skipinitialspace=True)
    jday = []
    storage = []
    storage_lt_52 = []
    storage_lt_60 = []
    for i,row in enumerate(strReader):
        if i>1: # 2 header rows
            jday.append(float(row[0]))
            storage.append(float(row[1])*m3_to_acft)
            storage_lt_52.append(float(row[2])*m3_to_acft)
            storage_lt_60.append(float(row[3])*m3_to_acft)
    return jday,storage,storage_lt_52,storage_lt_60
    

def write_shutter_elevations_to_output_dss(str_csv,vol_csv,dss_file,output_tsc):
    '''Take str CSV shutter elevation output from W2 Folsom (interval 7.5 minutes?), 
    which is not ouput at the same interval as W2 output data (2 hours last time I checked),
    and merge it to the output time interval for writing to dss.
    output_tsc is an example, and should have the correct W2 Folsom fpart'''

    jday,sElev1,sElev2,sElev3,bypassFlow,revisedPenstock1Flow = read_str_csv(str_csv)

    jday_output = DSS_Tools.jday_from_tsc(output_tsc)

    parts = str(output_tsc.fullName).split('/')
    epart = parts[5]
    W2_fpart = parts[6]

    dss_fm = HecDss.open(dss_file)
    output_tsc.units= 'ft'

    jday_output_offset = tz_offset.days # there is a timezone shift on the DSS records currently
    
    for i,elev in enumerate([sElev1,sElev2,sElev3]):

        merged_data = merge_data_nearest_jday(jday_output,jday,elev,jday_output_offset)
        elev_name = 'W2_Folsom_Forecast_Shutter_%i'%(i+1)
        output_tsc.fullName = '/'.join(['','',elev_name,'ELEV','',epart,W2_fpart,''])
        output_tsc.values = merged_data
        dss_fm.put(output_tsc)       

    # write bypass flow also
    output_tsc.units = 'cms'
    
    output_tsc.values = merge_data_nearest_jday(jday_output,jday,bypassFlow,jday_output_offset)
    bflow_name = 'W2_Folsom_Forecast_BypassFlow'
    output_tsc.fullName = '/'.join(['','',bflow_name,'FLOW','',epart,W2_fpart,''])
    dss_fm.put(output_tsc)
    
    output_tsc.values = merge_data_nearest_jday(jday_output,jday,revisedPenstock1Flow,jday_output_offset)
    bflow_name = 'W2_Folsom_Forecast_RevisedPenstock1Flow'
    output_tsc.fullName = '/'.join(['','',bflow_name,'FLOW','',epart,W2_fpart,''])    
    dss_fm.put(output_tsc) 


    # write storage also
    jday,storage,storage_lt_52,storage_lt_60 = read_storage_csv(vol_csv)
    output_tsc.units= 'ac-ft'
    for name,stor in zip(['W2_Folsom_Storage','W2_Folsom_Storage_lt_52F','W2_Folsom_Storage_lt_60F'],
                         [storage,storage_lt_52,storage_lt_60]):
        output_tsc.values = merge_data_nearest_jday(jday_output,jday,stor,jday_output_offset)        
        output_tsc.fullName = '/'.join(['','',name,'STOR','',epart,W2_fpart,''])
        dss_fm.put(output_tsc)        
        
    dss_fm.close()
    

def merge_data_nearest_jday(jday1, jday2, data2, jd1_offset=0):
    """
    Merges the data from datasets (jday2, data2) to jday1 based on the nearest jday.

    Args:
        jday1: List of jday values for the first dataset.
        jday2: List of jday values for the second dataset.
        data2: List of corresponding data values for the second dataset.
        jd1_offset: offset in days. Sometimes timezones are a terrible thing

    Returns:
        data2_jday1: merged data2 of len(jday1)
    """

    data2_jday1 = []
    for i, jday1_val in enumerate(jday1):
        min_diff = 1.0e9  # Initialize with positive infinity
        nearest_jday2_val = None
        corresponding_T2_val = None

        for j, jday2_val in enumerate(jday2):
            diff = abs(jday1_val+jd1_offset - jday2_val)
            if diff < min_diff:
                min_diff = diff
                nearest_jday2_val = jday2_val
                corresponding_data2_val = data2[j]

        data2_jday1.append(corresponding_data2_val)

    return data2_jday1

    

def rectify_tsc_dates_to_model_year(tsc,model_year):
    '''If you mess up these dates, the DSS write fails silently, and may mess up future writes while the file is open!
    '''
    ystr = str(model_year)

    new_hec_times = []
    for j in range(tsc.numberValues):
        # Assuming hectime can be converted to Java Date or has method to get the equivalent
        date_str = tsc.getHecTime(j).dateAndTime()  # NOT 05Jan2010 0000, actually '5 January 2010, 24:00'
        #print(date_str)
        #date_str = date_str[:5]+ystr+date_str[9:]
        date_str = date_str[:-11] + ystr + date_str[-7:]
        print(date_str)
        new_hec_times.append(HecTime(date_str).value())

    tsc.times = new_hec_times
    tsc.startTime = new_hec_times[0]

    return tsc

def write_constant_1day_ts(dssFm,rec,rtw,constant_value):
    starttime_str = rtw.getStartTimeString()
    #endtime_str = timewindow.getEndTimeString()

    st = HecTime(starttime_str).value()
    times = [st+60*i for i in range(24)]
    values = [constant_value for i in range(24)]
    
    tsc = TimeSeriesContainer()
    tsc.startTime = st
    tsc.times = times
    tsc.values = values
    tsc.units = '#'
    tsc.type = 'INST-VAL'
    tsc.interval = 60
    tsc.numberValues = len(values)
    tsc.fullName = rec
    rec_parts = rec.split('/')
    tsc.location = rec_parts[2]
    tsc.parameter = rec_parts[3]
    tsc.version = rec_parts[6]

    dssFm.put(tsc)

def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
 
    locations_obj = currentAlternative.getInputDataLocations()
    locations = DSS_Tools.organizeLocationsPaired(currentAlternative, locations_obj, W2FolsomFlowTempLocs, return_dss_paths=True)
    currentAlternative.addComputeMessage('Found DSS paths:')
    for location in locations:
        for path in location:
            currentAlternative.addComputeMessage(str(path))
    
    currentAlternative.addComputeMessage('\n')
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    dss_file = computeOptions.getDssFilename()
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()
    startyear_str = starttime_str[5:9]
    year = int(startyear_str)

    scripting_run_dir = computeOptions.getRunDirectory()
    run_dir,_ = os.path.split(scripting_run_dir)
    model_name = 'W2 Folsom'  # base Folsom forecast model (iterative w/ bypass)
    if "FixedATSP" in run_dir:
        model_name += " FixedATSP"
    if "NoBypass" in run_dir:
        model_name += " NoBypass"
    #model_dir = fpp.model_dir_from_run_dir(scripting_run_dir,'Folsom',model_name)
    w2_run_dir = os.path.join(run_dir,'cequal-w2','Folsom',model_name)
    str_csv = os.path.join(w2_run_dir,'str_br1.csv')
    vol_file = os.path.join(w2_run_dir,'VOLUME_WB1.OPT')
    
    # get (and potentially fix) scripted output locations
    outputlocations_obj = currentAlternative.getOutputDataLocations()
    outputlocations = DSS_Tools.organizeLocations(currentAlternative, outputlocations_obj, OutputLocs)
    outputpaths = []
    for ol in outputlocations:
        opath = currentAlternative.createOutputTimeSeries(ol)
        tspath = str(opath).split('/')
        fpart = tspath[6]
    #    new_fpart = fpart.lower().split(':scripting-')[0] #fpart will have scripting in it, but we want w2 version
    #    new_fpart += ':cequalw2-' #replace scripting with W2
    #    new_fpart += '-'.join(computeOptions.getSimulationName().split('-')[:-1]) #sim name will have the sim group, so snip that
        if '|' in fpart:
            # remove everything before | including |
            tspath[6] = fpart[fpart.find('|')+1:]
        #tspath[6] = 'W2_Folsom'
        outputpaths.append('/'.join(tspath))

    # get W2 Folsom F-part, because we want some of these output to be plottable under the W2 f-part
    W2_fpart = locations[0][0].split('/')[6]

    # Add Folsom outflows into Natoma
    # ------------------------------------------------------------------------------------
    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpaths[0]))
    flow_locations = [loc[0] for loc in locations]  
    DSS_Tools.add_flows(currentAlternative, rtw, flow_locations, dss_file, outputpaths[0], dss_file)

    # copy to W2 alt for plotting
    DSS_Tools.copy_dss_ts(outputpaths[0],new_fpart=W2_fpart,dss_file_path=dss_file)    
    
    # Flow-weight average Folsom outflow temps for linking
    # ------------------------------------------------------------------------------------
    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpaths[1]))
    cfs_limit = 1.0 #float
    tsc_natoma_inflow_temp = flowweightaverage.FWA2(currentAlternative, dss_file, rtw, locations, outputpaths[1], 
                                                    cfs_limit, 10.0, return_tsc=True)
    DSS_Tools.copy_dss_ts(outputpaths[1],new_fpart=W2_fpart,dss_file_path=dss_file)    

    # Write out target temps and schedule number
    # ------------------------------------------------------------------------------------
    run_dir = computeOptions.getRunDirectory() # runs/model_alt_name/scripting
    model_name = 'W2 Folsom'  # base Folsom forecast model (iterative w/ bypass)
    if "FixedATSP" in run_dir:
        model_name += " FixedATSP"
    if "NoBypass" in run_dir:
        model_name += " NoBypass"
    model_dir = fpp.model_dir_from_run_dir(run_dir,'Folsom',model_name)
    run_base_dir,_ = os.path.split(run_dir)
    model_run_dir_Folsom = os.path.join(run_base_dir,'CeQual-W2','Folsom',model_name)
    shared_dir = os.path.join(fpp.study_dir_from_run_dir(run_dir),'shared')
    forecast_dss = os.path.join(shared_dir,'WTMP_American_forecast.dss')
    
    # read TEMP_LOG.OPT and find last schedule used (matching ensemble number if Schedules are loaded in WTMP as target temps)
    # ----------------------------------------------------------------
    nSchedule = None
    with open(os.path.join(model_run_dir_Folsom,'TEMP_LOG.OPT'),'r') as fp:
        for line in fp.readlines():
            # the last schedule read in the file is the last (potentially iterative) schedule used
            if line.startswith("OPEN FILE:TargetSchedulesA.npt"):
                # extract schedule number from line that looks like this: OPEN FILE:TargetSchedulesA.npt   122.000   39    7
                # where 39 is the schedule  number
                tokens = line.split() # split on whitespace
                nSchedule = int(tokens[3])

    # write (copy) schedule target temps to output DSS
    # ----------------------------------------------------------------
    schedule_rec = "/WTMP_American/Folsom Iterative Target/TEMP-WATER//1Day/C:0000"+ "%02i"%nSchedule +"|ScheduleA/"
    dssIn = HecDss.open(os.path.join(shared_dir,"ScheduleA_TT_daily.dss"))
    tsc_sched = dssIn.get(schedule_rec,True)
    dssIn.close()

    # need to line up schedule with run time window ... hmm use doy of start date and find nearest date? Read the TT in output dss for dates?
    tsc_sched = rectify_tsc_dates_to_model_year(tsc_sched,year)
    # convert to C just so we don't go crazy using DSS
    if tsc_sched.units == 'F':
        tc = []
        for tf in tsc_sched.values:
            tc.append((tf-32.)*5.0/9.0)
        tsc_sched.values = tc
        tsc_sched.units = 'C'    
    dssOut = HecDss.open(dss_file)

    # write schedule number to output DSS
    # ----------------------------------------------------------------
    sched_parts = tsc_sched.fullName.split('/')
    out_parts = outputpaths[2].split('/')
    out_parts[5] = sched_parts[5] # get period from original record
    tsc_sched.fullName = '/'.join(out_parts)
    #tsc_sched.numberValues = len(tsc_sched.values)
    #tsc_sched.setStoreAsDoubles(True)
    print('writing W2 Schedule Target Temps: ')
    print('    '+dss_file)
    print('    '+tsc_sched.fullName)
    
    #print(tsc_sched.values)
    #print(len(tsc_sched.values))
    dssOut.put(tsc_sched)  # TODO:  Why does this never show up in output file???
    out_parts[6] = W2_fpart
    tsc_sched.fullName = '/'.join(out_parts)
    dssOut.put(tsc_sched) # write out copy under W2_fpart
    
    # use above tsc to write constant nSchedule to dss, I guess
    # ----------------------------------------------------------------
    print('writing W2 Schedule Number: '+outputpaths[3])
    print('    '+dss_file)
    print('    '+outputpaths[3])
    write_constant_1day_ts(dssOut,outputpaths[3],rtw,nSchedule)
    out_parts = outputpaths[3].split('/')
    out_parts[6] = W2_fpart
    write_constant_1day_ts(dssOut,'/'.join(out_parts),rtw,nSchedule)  # write out copy under W2_fpart

    # back calc dowstream temp using W2 regressions
    # ----------------------------------------------------------------
    currentAlternative.addComputeMessage("Back-calculating downstream temperature using W2 regressions...")
    doys,FaveFlow,Tair = fpp.load_tt_data(forecast_dss, starttime_str, endtime_str) # day-of-year,CMS,C    
    #TModel_tsm = DSS_Tools.dss_read_ts_safe(forecast_dss,dssRec,starttime_str,endtime_str)
    #TModel_tsc = currentAlternative.loadTimeSeries(locations[0])

    TModel_tsm = tsmath(tsc_natoma_inflow_temp)
    TModel_tsc_daily = DSS_Tools.standardize_interval(TModel_tsm,'1day').getData()
    dtt,DownstreamTempWatt = fpp.calc_downstream_temp_W2(year,1,doys,Tair,TModel_tsc_daily.values,FaveFlow,lagWatt=True)
    dtt,DownstreamTempHazel = fpp.calc_downstream_temp_W2(year,3,doys,Tair,TModel_tsc_daily.values,FaveFlow,lagWatt=True)

    out_parts = outputpaths[4].split('/')
    out_parts[5] = '1Day'
    out_parts[2] = 'W2_DownstreamRegressionWatt'
    TModel_tsc_daily.fullName = '/'.join(out_parts)

    print('writing W2 Downstream Regression calc: '+TModel_tsc_daily.fullName)
    #print(DownstreamTemp)
    TModel_tsc_daily.values = DownstreamTempWatt
    #print('    '+dss_file)
    #print('    '+TModel_tsc_daily.fullName)
    dssOut.put(TModel_tsc_daily)
    script_fpart = out_parts[6]
    out_parts[6] = W2_fpart
    TModel_tsc_daily.fullName = '/'.join(out_parts)
    dssOut.put(TModel_tsc_daily) # write out copy under W2_fpart

    out_parts[6] = script_fpart
    out_parts[2] = 'W2_DownstreamRegressionHazel'
    TModel_tsc_daily.fullName = '/'.join(out_parts)
    print('writing W2 Downstream Regression calc: '+TModel_tsc_daily.fullName)
    TModel_tsc_daily.values = DownstreamTempHazel
    dssOut.put(TModel_tsc_daily)
    out_parts[6] = W2_fpart
    TModel_tsc_daily.fullName = '/'.join(out_parts)
    dssOut.put(TModel_tsc_daily) # write out copy under W2_fpart

    # we need shuttler elevation for creating buzz plots
    # ----------------------------------------------------------------
    tsc = dssOut.get(str(locations[0][0]),True)
    print('DEBUG:','got data!!!')
    dssOut.close()    
    write_shutter_elevations_to_output_dss(str_csv,vol_file,dss_file,tsc)

    return True

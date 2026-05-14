import os
from com.rma.io import DssFileManagerImpl
from hec.heclib.dss import HecDss
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from hec.heclib.util import HecTime
import hec.hecmath.TimeSeriesMath as tsmath

import DSS_Tools
reload(DSS_Tools)

import Forecast_preprocess as fpp
reload(fpp)

# NOTE: As of 2025-02 this script/model alternative is unused. For W2_Folsom - only models, this script have been
# replaced by theoutput link sciprt for folsom for consistency. -ben saenz


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
def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )

    rtw = computeOptions.getRunTimeWindow()
    dss_file = computeOptions.getDssFilename()
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()
    startyear_str = starttime_str[5:9]
    year = int(startyear_str)

    ### Determine the alternative name ###
    # This step is necessary to set the model up correctly based on American configurations
    # Get the current folder
    study_root = computeOptions.getRunDirectory()  # returns watershed directory
    
    # Split the path and strip off the scripts, alternative-group, and run folder to obtain the WAT directory
    study_main = study_root.split(os.sep)
    study_main = os.sep.join(study_main[:-3])
    
    # Create the directory to the WAT simulation group
    study_wat = study_main + os.sep + 'wat' + os.sep + 'simGroups'
    
    # Get the contents of the directory
    all_groups = os.listdir(study_wat)
    
    # Filter to only forecast groups
    forecast_groups = [x for x in all_groups if ".fsimgrp" in x]
    
    # Remove backup files not caught by the previous filter
    forecast_groups = [x for x in forecast_groups if '.bak' not in x]
    
    # Strip off the file extensions from each name
    forecast_groups = [x.split(".fsimgrp")[0] for x in forecast_groups]
    
    # Sort by the list decreasing
    forecast_groups.sort(reverse=True)
    
    # Get the simulation name which is in alternative-simulation group format separated by a dash
    sim_name = computeOptions.getSimulationName()
    
    # Loop and find the what forecast group is within the simulation name
    for group_name in forecast_groups:
        # Check if the group name is in the simulation name
        if group_name in sim_name:
            # Valid group name component has been found. Break to continue processing 
            break
    
    # Remove the group from the simulation name
    alt_name = sim_name.split(group_name)[0]
    
    # Remove the dash from the alternative name
    alt_name = alt_name[:-1]
    
    # Construct the new simulation name that is the disk name
    alt_name_underscore = alt_name.replace(' ', '_')
    sim_name_underscore = alt_name_underscore + '-' + group_name

    ### Adjust the model name based on the alternative name ###
    # Get the current run directory
    run_dir = computeOptions.getRunDirectory()

    # Start with the base model name. This is an assumed convention. 
    model_name = "W2 Folsom"  # base Folsom forecast model (iterative w/ bypass)
    
    # Check the alternative if the model needs adjustment for the fixed ATSP mode
    if "FixedATSP" in alt_name:
        # FixedATSP is found in the alternative name. Append it to the run directory
        model_name += " FixedATSP"
    
    # Check the alternative if the model needs adjustment for the no bypass mode
    if "NoBypass" in alt_name:
        # NoBypass is found in the alternative name. Append it to the run directory
        model_name += " NoBypass"
    
    ### Setup the paths ###
    model_dir = os.path.join(study_main, 'runs', sim_name_underscore, 'CeQual-W2', 'Folsom', model_name)
    targt_temp_npt_filepath = os.path.join(model_dir,'TargetSchedulesA.npt') # overwrite what's there    
    shared_dir = os.path.join(study_main,'shared')
    forecast_dss = os.path.join(shared_dir,'WTMP_American_forecast.dss')
    schedule_csv = os.path.join(model_dir,'SchedulesA.csv')

    locations = currentAlternative.getInputDataLocations() # should be only one, W2 Folsom outflow temp
    locations_path = str(currentAlternative.loadTimeSeries(locations[0]))

    print('locations_path:')
    print(locations_path)

    # get W2 Folsom F-part from input locations    
    W2_fpart = locations_path.split('/')[6]

    # get scripted output locations
    outputlocations = currentAlternative.getOutputDataLocations()
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
        
    # read TEMP_LOG.OPT and find last schedule used (matching ensemble number if Schedules are loaded in WTMP as target temps)
    nSchedule = None
    with open(os.path.join(model_dir,'TEMP_LOG.OPT'),'r') as fp:
        for line in fp.readlines():
            # the last schedule read in the file is the last (potentially iterative) schedule used
            if line.startswith("OPEN FILE:TargetSchedulesA.npt"):
                # extract schedule number from line that looks like this: OPEN FILE:TargetSchedulesA.npt   122.000   39    7
                # where 39 is the schedule  number
                tokens = line.split() # split on whitespace
                nSchedule = int(tokens[3])

    # write (copy) schedule target temps to output DSS
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
    sched_parts = tsc_sched.fullName.split('/')
    out_parts = outputpaths[0].split('/')
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
    
    # use above tsc to write constant nSchedule to dss, I guess
    print('writing W2 Schedule Number: '+outputpaths[1])
    print('    '+dss_file)
    print('    '+outputpaths[1])
    write_constant_1day_ts(dssOut,outputpaths[1],rtw,nSchedule)
    

    # back calc dowstream temp using W2 regressions
    currentAlternative.addComputeMessage("Back-calculating downstream temperature using W2 regressions...")
    doys,FaveFlow,Tair = fpp.load_tt_data(forecast_dss, starttime_str, endtime_str) # day-of-year,CMS,C    
    #TModel_tsm = DSS_Tools.dss_read_ts_safe(forecast_dss,dssRec,starttime_str,endtime_str)
    TModel_tsc = currentAlternative.loadTimeSeries(locations[0])
    TModel_tsm = tsmath(TModel_tsc)
    TModel_tsc_daily = DSS_Tools.standardize_interval(TModel_tsm,'1day').getData()
    dtt,DownstreamTemp = fpp.calc_downstream_temp_W2(year,1,doys,Tair,TModel_tsc_daily.values,FaveFlow,lagWatt=False)

    out_parts = outputpaths[2].split('/')
    out_parts[5] = '1Day'
    TModel_tsc_daily.fullName = '/'.join(out_parts)
    print('writing W2 Downstream Regression calc: '+TModel_tsc_daily.fullName)
    #print(DownstreamTemp)
    TModel_tsc_daily.values = DownstreamTemp
    print('    '+dss_file)
    print('    '+TModel_tsc_daily.fullName)
    dssOut.put(TModel_tsc_daily)
    dssOut.close()
    
    return True


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
    

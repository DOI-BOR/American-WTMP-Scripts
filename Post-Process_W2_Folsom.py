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

# NOTE: As of 2025-02 this script/model alternative is unused.
#
# Historical Context
# ------------------
# This script originally provided post-processing for standalone
# W2_Folsom forecast alternatives.
#
# It has since been replaced by the Folsom output-link script so that
# all Folsom workflows use a consistent post-processing pathway.
#
# Maintenance Guidance
# --------------------
# Before modifying this script, verify whether any legacy forecast
# configurations still reference it. The primary production workflow
# may no longer execute this code.
#
# Original note:
# "For W2_Folsom-only models, this script has been replaced by the
# output link script for Folsom for consistency."

##
#
# computeAlternative function is called when the ScriptingAlternative
# is computed by HEC-WAT.
#
# Arguments:
#
#   currentAlternative
#       ScriptingAlternative object.
#       hec2.wat.plugin.java.impl.scripting.model.ScriptPluginAlt
#
#   computeOptions
#       Compute configuration object.
#       hec.wat.model.ComputeOptions
#
# Return:
#
#   True
#       Successful execution.
#
#   False
#       Failure.
#
# Notes:
#
#   No explicit return value is treated as success by HEC-WAT.
#
# Purpose
# -------
#
# This script:
#
#   1. Identifies the final W2 schedule used during model execution.
#   2. Writes schedule target temperatures to DSS.
#   3. Writes active schedule number to DSS.
#   4. Generates downstream temperature regression outputs.
#
# These products are primarily forecast diagnostics rather than
# direct model outputs.
#
##


def computeAlternative(currentAlternative, computeOptions):
    
    #######################################################################
    # Script initialization
    #######################################################################

    # Record execution in the HEC-WAT compute log.
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )

    #######################################################################
    # Runtime configuration
    #
    # These values define the forecast period and output DSS file.
    #######################################################################

    # Forecast runtime window.
    rtw = computeOptions.getRunTimeWindow()
    
    # DSS output file for the current simulation.
    dss_file = computeOptions.getDssFilename()
    
    # Runtime start and end timestamps.
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()
    
    # Extract forecast year.
    #
    # Schedule libraries often contain historical years and must later
    # be aligned with the current forecast year.
    startyear_str = starttime_str[5:9]
    year = int(startyear_str)

    #######################################################################
    # Determine alternative name and simulation group
    #
    # Forecast directory structures are derived from simulation names.
    #
    # HEC-WAT simulation names typically follow:
    #
    #     <AlternativeName>-<ForecastGroup>
    #
    # This section reconstructs both components.
    #######################################################################
    
    # Run directory provided by HEC-WAT.
    study_root = computeOptions.getRunDirectory()  # returns watershed directory
    
    # Remove run-specific folders and obtain study root
    study_main = study_root.split(os.sep)
    study_main = os.sep.join(study_main[:-3])
    
    # Simulation group definition directory.
    study_wat = study_main + os.sep + 'wat' + os.sep + 'simGroups'
    
    # Read available simulation groups.
    all_groups = os.listdir(study_wat)
    
    # Retain forecast groups only.
    forecast_groups = [x for x in all_groups if ".fsimgrp" in x]
    
    # Ignore backup files.
    forecast_groups = [x for x in forecast_groups if '.bak' not in x]
    
    # Remove filename extension.
    forecast_groups = [x.split(".fsimgrp")[0] for x in forecast_groups]
    
    # Reverse sorting favors newer forecast groups.
    forecast_groups.sort(reverse=True)
    
    # Simulation name supplied by HEC-WAT.
    sim_name = computeOptions.getSimulationName()
    
    #######################################################################
    # Identify forecast group
    #
    # Assumes forecast group name appears within simulation name.
    #######################################################################
    
    for group_name in forecast_groups:
        
        # Check if the group name is in the simulation name
        if group_name in sim_name:
            
            # Valid group name component has been found. Break to continue processing 
            break
    
    # Remove simulation group suffix.
    alt_name = sim_name.split(group_name)[0]
    
    # Remove trailing separator
    alt_name = alt_name[:-1]
    
    # Convert to directory naming convention.
    alt_name_underscore = alt_name.replace(' ', '_')
    
    # Construct run-directory simulation name.
    sim_name_underscore = alt_name_underscore + '-' + group_name

    #######################################################################
    # Determine CE-QUAL-W2 model variant
    #
    # Alternative naming conventions control model selection.
    #######################################################################

    # Current run directory.
    run_dir = computeOptions.getRunDirectory()

    # Default model.
    model_name = "W2 Folsom"  # base Folsom forecast model (iterative w/ bypass)
    
    # Fixed ATSP variant.
    if "FixedATSP" in alt_name:
        
        # FixedATSP is found in the alternative name. Append it to the run directory
        model_name += " FixedATSP"
    
    # Check the alternative if the model needs adjustment for the no bypass mode
    if "NoBypass" in alt_name:
        
        # NoBypass is found in the alternative name. Append it to the run directory
        model_name += " NoBypass"
    
    #######################################################################
    # Construct model file locations
    #
    # Expected directory structure:
    #
    # runs/
    #   <simulation>/
    #     CeQual-W2/
    #       Folsom/
    #         <model_name>/
    #######################################################################
    
    model_dir = os.path.join(study_main, 'runs', sim_name_underscore, 'CeQual-W2', 'Folsom', model_name)
    
    # Historical target schedule file.
    #
    # Note:
    # This variable is not used later in the script but documents the
    # expected model schedule location.
    targt_temp_npt_filepath = os.path.join(model_dir,'TargetSchedulesA.npt') # overwrite what's there    
    
    # Shared forecast resources.
    shared_dir = os.path.join(study_main,'shared')
    
    # Forecast DSS containing meteorological inputs.
    forecast_dss = os.path.join(shared_dir,'WTMP_American_forecast.dss')
    
    # CSV schedule export generated by model workflows.
    schedule_csv = os.path.join(model_dir,'SchedulesA.csv')
    
    # Load the input DSS locations for the current alternative.
    locations = currentAlternative.getInputDataLocations() # should be only one, W2 Folsom outflow temp
    
    # Convert the first input time series location to a string path.
    locations_path = str(currentAlternative.loadTimeSeries(locations[0]))

    print('locations_path:')
    print(locations_path)

    # Extract the F-part from the W2 Folsom input DSS path.
    # DSS pathname parts are slash-delimited; index 6 is the F-part.    
    W2_fpart = locations_path.split('/')[6]

    # Get the scripted output DSS locations for the current alternative.
    outputlocations = currentAlternative.getOutputDataLocations()
    outputpaths = []
    
    # Create output time series records and normalize their F-parts.
    for ol in outputlocations:
        opath = currentAlternative.createOutputTimeSeries(ol)
        tspath = str(opath).split('/')
        fpart = tspath[6]
   
        # If the F-part contains a pipe, keep only the portion after it.
        if '|' in fpart:
            
            # remove everything before | including |
            tspath[6] = fpart[fpart.find('|')+1:]
        
        # Reassemble the cleaned pathname and save it.
        outputpaths.append('/'.join(tspath))
        
    # Read TEMP_LOG.OPT to determine the last schedule number used.
    # This matches the final iterative schedule loaded by the model.
    nSchedule = None
    with open(os.path.join(model_dir,'TEMP_LOG.OPT'),'r') as fp:
        for line in fp.readlines():
            
            # The last matching line in the file is the last schedule applied.
            if line.startswith("OPEN FILE:TargetSchedulesA.npt"):
                # extract schedule number from line that looks like this: OPEN FILE:TargetSchedulesA.npt   122.000   39    7
                # where 39 is the schedule  number
                # Split the log line on whitespace and parse the schedule number.
                tokens = line.split() # split on whitespace
                nSchedule = int(tokens[3])

    # Build the DSS record path for the selected target temperature schedule.
    schedule_rec = "/WTMP_American/Folsom Iterative Target/TEMP-WATER//1Day/C:0000"+ "%02i"%nSchedule +"|ScheduleA/"
    
    # Open the schedule DSS file and retrieve the matching time series.
    dssIn = HecDss.open(os.path.join(shared_dir,"ScheduleA_TT_daily.dss"))
    tsc_sched = dssIn.get(schedule_rec,True)
    dssIn.close()

    # Align the schedule dates with the model run year before writing.
    # This avoids silent DSS write failures due to date mismatches.
    tsc_sched = rectify_tsc_dates_to_model_year(tsc_sched,year)
    
    # Convert schedule temperatures from Fahrenheit to Celsius if needed.
    if tsc_sched.units == 'F':
        tc = []
        for tf in tsc_sched.values:
            tc.append((tf-32.)*5.0/9.0)
        tsc_sched.values = tc
        tsc_sched.units = 'C' 

    # Open the model output DSS file for writing.
    dssOut = HecDss.open(dss_file)

    # Update the schedule record pathname to match the output location.
    # Preserve the original period part from the schedule record.
    sched_parts = tsc_sched.fullName.split('/')
    out_parts = outputpaths[0].split('/')
    out_parts[5] = sched_parts[5] # get period from original record
    tsc_sched.fullName = '/'.join(out_parts)
    
    print('writing W2 Schedule Target Temps: ')
    print('    '+dss_file)
    print('    '+tsc_sched.fullName)
    
    # Write the schedule target temperature series to the output DSS.
    dssOut.put(tsc_sched)  # TODO:  Why does this never show up in output file???
    
    # Write a constant time series carrying the selected schedule number.
    # Reuse the run time window and second output path.
    print('writing W2 Schedule Number: '+outputpaths[1])
    print('    '+dss_file)
    print('    '+outputpaths[1])
    write_constant_1day_ts(dssOut,outputpaths[1],rtw,nSchedule)
    

    # Back-calculate downstream temperature using W2 regression methods.
    currentAlternative.addComputeMessage("Back-calculating downstream temperature using W2 regressions...")
    
    # Load forecast day-of-year, average flow, and air temperature inputs.
    doys,FaveFlow,Tair = fpp.load_tt_data(forecast_dss, starttime_str, endtime_str) # day-of-year,CMS,C    
    
    # Load the model output temperature time series from the input location.
    TModel_tsc = currentAlternative.loadTimeSeries(locations[0])
    TModel_tsm = tsmath(TModel_tsc)
    
    # Standardize the model temperature series to a daily interval.
    TModel_tsc_daily = DSS_Tools.standardize_interval(TModel_tsm,'1day').getData()
    
    # Compute downstream temperatures using the W2 regression routine.
    dtt,DownstreamTemp = fpp.calc_downstream_temp_W2(year,1,doys,Tair,TModel_tsc_daily.values,FaveFlow,lagWatt=False)

    # Prepare the downstream temperature output pathname with a 1Day interval.
    out_parts = outputpaths[2].split('/')
    out_parts[5] = '1Day'
    
    # Replace the daily model temperatures with the computed downstream values.
    TModel_tsc_daily.fullName = '/'.join(out_parts)
    print('writing W2 Downstream Regression calc: '+TModel_tsc_daily.fullName)
    
    # Replace the daily model temperatures with the computed downstream values.
    TModel_tsc_daily.values = DownstreamTemp
    print('    '+dss_file)
    print('    '+TModel_tsc_daily.fullName)
    
    # Write the downstream regression results and close the DSS file.
    dssOut.put(TModel_tsc_daily)
    dssOut.close()
    
    # Indicate successful completion of the workflow.
    return True


def rectify_tsc_dates_to_model_year(tsc,model_year):
    '''
    If you mess up these dates, the DSS write fails silently, and may mess up future writes 
    while the file is open!
    '''
    
    # Convert the model year to a string for date replacement.
    ystr = str(model_year)

    # Collect updated HEC times matching the model year.
    new_hec_times = []
    for j in range(tsc.numberValues):
        # Read each timestamp as a formatted date string from the time series.
        # Assuming hectime can be converted to Java Date or has method to get the equivalent
        date_str = tsc.getHecTime(j).dateAndTime()  # NOT 05Jan2010 0000, actually '5 January 2010, 24:00'
        
        # Replace the year portion while preserving month, day, and time.
        date_str = date_str[:-11] + ystr + date_str[-7:]
        print(date_str)
        new_hec_times.append(HecTime(date_str).value())
    
    # Overwrite the time series timestamps and update the start time.
    tsc.times = new_hec_times
    tsc.startTime = new_hec_times[0]
    
    # Return the modified time series container.
    return tsc

def write_constant_1day_ts(dssFm,rec,rtw,constant_value):
    
    # Get the run window start time string from the runtime window object.
    starttime_str = rtw.getStartTimeString()

    # Convert the start time to HEC time and create 24 hourly timestamps.
    st = HecTime(starttime_str).value()
    times = [st+60*i for i in range(24)]
    values = [constant_value for i in range(24)]
    
    # Create a DSS time series container for the constant output series.
    tsc = TimeSeriesContainer()
    tsc.startTime = st
    tsc.times = times
    tsc.values = values
    tsc.units = '#'
    tsc.type = 'INST-VAL'
    tsc.interval = 60
    tsc.numberValues = len(values)
    tsc.fullName = rec
    
    # Populate pathname metadata fields from the DSS record string.
    rec_parts = rec.split('/')
    tsc.location = rec_parts[2]
    tsc.parameter = rec_parts[3]
    tsc.version = rec_parts[6]

    # Write the constant time series to the DSS file manager.
    dssFm.put(tsc)
    

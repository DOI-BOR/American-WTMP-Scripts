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
#
# Arguments:
#   currentAlternative - the ScriptingAlternative.
#                        hec2.wat.plugin.java.impl.scripting.model.ScriptPluginAlt
#   computeOptions     - the compute options.
#                        hec.wat.model.ComputeOptions
#
# Return:
#   True  = script executed successfully
#   False = script failed
#
# Note:
#   If no explicit return value is provided, the script is treated as successful.
#
##

# Common location methods:
# location.getName()
# location.getParameter()

# Mapping between flow locations and corresponding temperature locations.
# These locations are used when processing W2 Folsom model outputs.

W2FolsomFlowTempLocs = [
  ['Qwo-105-str1','Two-105-str1'],
  ['Qwo-105-str2','Two-105-str2'],
  ['Qwo-105-str3','Two-105-str3'],
  ['Qwo-105-str4','Two-105-str4'],
  ['Qwo-105-str6','Two-105-str6'], # Penstock 5 is not used
  ['Qwo-105-str7','Two-105-str7'],
]

# Extract only the flow location names.
W2FolsomFlowLocs = [ft[0] for ft in W2FolsomFlowTempLocs]

# DSS output locations written by this workflow.
OutputLocs = [
    'W2_Natoma_InflowQ',
    'W2_Natoma_InflowT',
    'W2_FOLSOM_SCHEDULE_TEMP_FINAL',
    'W2_FOLSOM_SCHEDULE_FINAL',
    'W2_Folsom_Downstream_Temp'
]



def snap_folsom_elev(elev_in):
    """
    Maps a Folsom Dam intake elevation (ft) to the nearest known operational
    shutter intake elevation, removing small floating-point variations introduced
    by the CE-QUAL-W2 model output.

    Downstream plotting and analysis tools identify shutter configurations by
    comparing elevations to fixed reference values. Without snapping, minor
    floating-point drift in W2 output can cause misidentification.

    Inputs:
      elev_in -- float; Folsom intake elevation in feet (converted from W2 meters output)

    Output:
      Returns a float representing the snapped operational intake elevation (ft):
        >= 395 ft -> 401.0 ft  (all shutters installed)
        >= 355 ft -> 362.0 ft  (one shutter removed)
        >= 346 ft -> 349.0 ft  (deganged-middle shutter)
        >= 320 ft -> 336.0 ft  (intermediate configuration)
        >= 300 ft -> 307.0 ft  (lower / penstock configuration)
        <  300 ft -> elev_in   (passthrough; no defined snap position)
    """
    
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
    """
    Reads a CE-QUAL-W2 Folsom structure output file (str_br1.csv) and returns
    Julian day, snapped shutter elevations for three penstocks, bypass flow,
    and a revised Penstock 1 flow that accounts for the W2 bypass mode.

    In W2 bypass mode, Penstock 1 is repurposed as a bypass outlet by lowering
    its intake elevation to exactly 64.01 m. When this occurs:
      - The flow reported in column 8 (Penstock 1) is added to column 13 (bypass)
        to produce the total bypass flow.
      - The revised Penstock 1 flow is set to zero (it is no longer a penstock).

    Shutter elevations are converted from meters to feet and snapped to known
    operational values using snap_folsom_elev. The first two rows are header rows
    and are skipped.

    Inputs:
      csv_file_path -- full path to the W2 structure output CSV file (str_br1.csv)

    Output:
      Returns a tuple (jday, sElev1, sElev2, sElev3, bypassFlow, revisedPenstock1Flow)
      where each element is a list of floats with one entry per output time step:
        jday                 -- Julian day of each output row
        sElev1               -- snapped Penstock 1 intake elevation (ft)
        sElev2               -- snapped Penstock 2 intake elevation (ft)
        sElev3               -- snapped Penstock 3 intake elevation (ft)
        bypassFlow           -- total bypass flow (CMS; includes Penstock 1 flow in bypass mode)
        revisedPenstock1Flow -- Penstock 1 flow (CMS; zero when bypass mode is active)
    """
    
    
    strReader = csv.reader(open(csv_file_path), delimiter=',', quotechar='|')
    jday = []
    sElev1 = []
    sElev2 = []
    sElev3 = []
    
    bypassFlow = []
    revisedPenstock1Flow = []
    
    for i,row in enumerate(strReader):
        
        # Skip two header rows
        if i>1: 
            jday.append(float(row[0]))
            
            # Penstock 1 intake elevation (meters)
            e1 = float(row[15])

            # Detect W2 bypass mode
            if e1 == 64.01:
                bypassFlow.append(float(row[8])+float(row[13]))
                revisedPenstock1Flow.append(0.0)
            else:
                bypassFlow.append(float(row[13]))
                revisedPenstock1Flow.append(float(row[8]))
            
            # Convert meters to feet and snap to known elevations
            sElev1.append(snap_folsom_elev(e1*3.28084)) # m to ft
            sElev2.append(snap_folsom_elev(float(row[16])*3.28084)) # m to ft
            sElev3.append(snap_folsom_elev(float(row[17])*3.28084)) # m to ft            

    return jday,sElev1,sElev2,sElev3,bypassFlow,revisedPenstock1Flow

def read_storage_csv(csv_file_path):
    """
    Reads a CE-QUAL-W2 reservoir volume/storage output file (VOLUME_WB1.OPT) and
    returns Julian day, total storage, and two temperature-stratified storage
    volumes, all converted from cubic meters to acre-feet.

    The file is space-delimited with two header rows that are skipped. Column
    assignments are: [0] Julian day, [1] total storage, [2] storage below 52 deg F,
    [3] storage below 60 degF.

    Inputs:
      csv_file_path - full path to the W2 volume output file (VOLUME_WB1.OPT)

    Output:
      Returns a tuple (jday, storage, storage_lt_52, storage_lt_60) where each
      element is a list of floats with one entry per output time step (acre-ft):
        jday          - Julian day of each output row
        storage       - total reservoir storage (acre-ft)
        storage_lt_52 - storage volume below 52 degF (acre-ft)
        storage_lt_60 - storage volume below 60 degF (acre-ft)
    """
    
    # factor to convert to ac-ft 
    m3_to_acft = 0.000810714 
    
    strReader = csv.reader(open(csv_file_path), delimiter=' ', skipinitialspace=True)
    
    jday = []
    storage = []
    storage_lt_52 = []
    storage_lt_60 = []
    
    for i,row in enumerate(strReader):
        
        # Skip two header rows
        if i>1: 
            jday.append(float(row[0]))
            storage.append(float(row[1])*m3_to_acft)
            storage_lt_52.append(float(row[2])*m3_to_acft)
            storage_lt_60.append(float(row[3])*m3_to_acft)
    return jday,storage,storage_lt_52,storage_lt_60
    

def write_shutter_elevations_to_output_dss(str_csv,vol_csv,dss_file,output_tsc):
    '''
    Reads W2 shutter elevation and storage output files and writes
    processed results to DSS.

    Because shutter elevation output is typically written at a different
    time interval than other W2 model outputs, this routine merges the
    shutter data onto the DSS output time grid before writing.

    Outputs include:
      - Shutter elevations (Penstocks 1-3) (ELEV, ft)
      - Bypass flow (FLOW, CMS)
      - Revised Penstock 1 flow (FLOW, CMS)
      - Reservoir storage (STOR, acre-ft)
      - Storage below 52F (STOR, acre-ft)
      - Storage below 60F (STOR, acre-ft)

    output_tsc is used as a template container for DSS writes and
    provides the appropriate W2 forecast F-part.
    '''

    # Read structure output file
    jday,sElev1,sElev2,sElev3,bypassFlow,revisedPenstock1Flow = read_str_csv(str_csv)

    # Retrieve DSS output time grid
    jday_output = DSS_Tools.jday_from_tsc(output_tsc)

    # Parse DSS pathname components
    parts = str(output_tsc.fullName).split('/')
    epart = parts[5]
    W2_fpart = parts[6]

    # Open the DSS file for output
    dss_fm = HecDss.open(dss_file)
    
    # Write elevations in feet
    output_tsc.units= 'ft'

     # Optional timezone adjustment
    jday_output_offset = 0
    
    # Output each of penstock 2, penstock 3, and the leakage terms
    for i,elev in enumerate([sElev1,sElev2,sElev3]):

        merged_data = merge_data_nearest_jday(jday_output,jday,elev,jday_output_offset)
        
        elev_name = 'W2_Folsom_Forecast_Shutter_%i'%(i+1)
        
        output_tsc.fullName = '/'.join(['','',elev_name,'ELEV','',epart,W2_fpart,''])
        output_tsc.values = merged_data
        dss_fm.put(output_tsc)       

    # Reset the units
    output_tsc.units = 'cms'
    
    # Output the bypass flows
    output_tsc.values = merge_data_nearest_jday(jday_output,jday,bypassFlow,jday_output_offset)
    bflow_name = 'W2_Folsom_Forecast_BypassFlow'
    output_tsc.fullName = '/'.join(['','',bflow_name,'FLOW','',epart,W2_fpart,''])
    dss_fm.put(output_tsc)
    
    # Revised Penstock 1 flow
    output_tsc.values = merge_data_nearest_jday(jday_output,jday,revisedPenstock1Flow,jday_output_offset)
    bflow_name = 'W2_Folsom_Forecast_RevisedPenstock1Flow'
    output_tsc.fullName = '/'.join(['','',bflow_name,'FLOW','',epart,W2_fpart,''])    
    dss_fm.put(output_tsc) 
    
    # Write storage 
    jday,storage,storage_lt_52,storage_lt_60 = read_storage_csv(vol_csv)
    output_tsc.units= 'ac-ft'
    
    for name,stor in zip(['W2_Folsom_Storage','W2_Folsom_Storage_lt_52F','W2_Folsom_Storage_lt_60F'],
                         [storage,storage_lt_52,storage_lt_60]):
        
        output_tsc.values = merge_data_nearest_jday(jday_output,jday,stor,jday_output_offset)        
        output_tsc.fullName = '/'.join(['','',name,'STOR','',epart,W2_fpart,''])
        dss_fm.put(output_tsc)        
    
    # Close the DSS file
    dss_fm.close()
    

def merge_data_nearest_jday(jday1, jday2, data2, jd1_offset=0):
    """
    Maps values from a source Julian-day time grid (jday2/data2) onto a target
    time grid (jday1) using a nearest-forward lookup.

    For each Julian day in jday1, the function finds the first Julian day in jday2
    that is greater than or equal to it and assigns the corresponding data2 value.
    This is a forward-fill nearest-neighbor approach, not an interpolation.

    When no jday2 value satisfies the condition for a given jday1 entry (i.e., jday1
    extends beyond the end of jday2), the corresponding output value will be None.

    Inputs:
      jday1      -- list of floats; target time grid (Julian days) to map data onto
      jday2      -- list of floats; source time grid (Julian days) of the input data
      data2      -- list of floats; source values corresponding to jday2; same length as jday2
      jd1_offset -- float; optional offset added to jday1 values before lookup,
                    used for timezone corrections (default 0)

    Output:
      Returns a list with the same length as jday1, containing values from data2
      aligned to the jday1 time grid. Entries where no matching jday2 value was
      found will be None.
    """

    # Create a list to hold the ouput data
    data2_jday1 = []
    
    # Loop on the dates of the list to be filled to find the correct interpolation values
    for i, jday1_val in enumerate(jday1):
        
        # Initialize placeholder values for the nearest date and value
        nearest_jday2_val = None
        corresponding_data2_val = None
        
        # Loop on the dates of the seccond list to find the first date greater or equal to the target date
        for j, jday2_val in enumerate(jday2):
            
            # If the target date is less than or equal to the current value, update the value
            if jday1_val <= jday2_val:
                # Update the value and date
                nearest_jday2_val = jday2_val
                corresponding_data2_val = data2[j]
                # Break from the loop with the value
                break

        # Append the value into the output list
        data2_jday1.append(corresponding_data2_val)
    
    # Return to the calling function
    return data2_jday1

    

def rectify_tsc_dates_to_model_year(tsc,model_year):
    """
    Replaces the year component of all timestamps in a TimeSeriesContainer
    with a specified model year, preserving month, day, and time-of-day.

    This is needed when reading schedule temperature records from a library DSS
    file that may contain timestamps from a different year than the current
    forecast. DSS writes can fail silently when timestamps are malformed, so
    this function ensures all timestamps are valid for the forecast year before
    the container is written.

    Inputs:
      tsc        -- HEC TimeSeriesContainer object whose timestamps are to be updated
      model_year -- integer; the four-digit calendar year to substitute into all
                    timestamps

    Output:
      Returns the modified TimeSeriesContainer with all timestamps updated to
      model_year and startTime set to the first updated timestamp.
    """
    '''
    Replaces the year component of all timestamps in a
    TimeSeriesContainer.

    Important:
      DSS writes can fail silently if timestamps are malformed.
      This helper updates all timestamps to the model year while
      preserving month/day/time information.
    '''
    
    ystr = str(model_year)

    new_hec_times = []
    for j in range(tsc.numberValues):
        
        # Assuming hectime can be converted to Java Date or has method to get the equivalent
        date_str = tsc.getHecTime(j).dateAndTime()  # NOT 05Jan2010 0000, actually '5 January 2010, 24:00'
        date_str = date_str[:-11] + ystr + date_str[-7:]
        print(date_str)
        new_hec_times.append(HecTime(date_str).value())

    tsc.times = new_hec_times
    tsc.startTime = new_hec_times[0]

    return tsc

def write_constant_1day_ts(dssFm,rec,rtw,constant_value):
    """
    Writes a 24-hour constant-value time series (one value per hour, hourly interval)
    to a DSS file using an open HecDss file handle.

    Used to expose scalar values such as the selected ATSP schedule number as a
    DSS time series so that reporting tools can read them in a standard format.
    The time series starts at the run time window start time and spans 24 hourly
    time steps. Units are set to '#' (dimensionless count) and type to INST-VAL.

    Inputs:
      dssFm          -- open HecDss file handle to write the record into
      rec            -- DSS pathname string for the output record
      rtw            -- WAT run time window object providing the start time string
      constant_value -- float or int; the constant value to fill all 24 time steps

    Output:
      No return value. Writes one 24-step hourly INST-VAL DSS record to dssFm
      at the pathname specified by rec.
    """
    
    starttime_str = rtw.getStartTimeString()

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

def nSchedule_from_AutoRunTempLog(model_run_dir_Folsom):
    """
    Reads AutoRunTempLog.opt and returns the last valid ATSP
    schedule used during the compliance season (May 1 to Nov 30).

    The AutoRunTempLog.opt file records the schedule loaded at each time step during
    iterative W2 simulations. The function scans all rows and returns the schedule
    number from the last row whose Julian day is at or before 333 (approximately
    Nov 30), ignoring the reload that occurs after the season ends on day 334.

    This is the preferred method for determining the active schedule. Use this
    function instead of nSchedule_from_TEMP_LOG.

    Inputs:
      model_run_dir_Folsom -- full path to the CE-QUAL-W2 Folsom model run directory
                              containing AutoRunTempLog.opt

    Output:
      Returns an integer representing the last valid ATSP schedule number used during
      the compliance season, or None if no qualifying row is found.
    """
    
    nSchedule = None
    
    with open(os.path.join(model_run_dir_Folsom,'AutoRunTempLog.opt'),'r') as fp:
        for i,line in enumerate(fp.readlines()):
            if i>0:
                
                # the last schedule read in the file is the last (potentially iterative) schedule used May1-Nov30
                tokens = line.split(',') 
                
                # Schedule seems to revert to original on jday 334
                if float(tokens[0].strip()) <= 333.0:  
                    nSchedule = int(tokens[5].strip())
    return nSchedule


def nSchedule_from_TEMP_LOG(model_run_dir_Folsom):
    """
    Reads the CE-QUAL-W2 TEMP_LOG.OPT file and returns the last ATSP schedule
    number loaded during the compliance season (May 1 to Nov 30).

    DEPRECATED: This function is retained for reference only. TEMP_LOG.OPT is
    not always reliable for this purpose. Use nSchedule_from_AutoRunTempLog instead.

    The function scans all lines beginning with 'OPEN FILE:TargetSchedulesA.npt'
    and extracts the schedule number and day-of-year from the whitespace-delimited
    tokens. Only rows with a load day-of-year below 333 are considered, to avoid
    the post-season schedule reload.

    Inputs:
      model_run_dir_Folsom -- full path to the CE-QUAL-W2 Folsom model run directory
                              containing TEMP_LOG.OPT

    Output:
      Returns an integer representing the last valid ATSP schedule number found
      before Julian day 333, or None if no qualifying line is found.
    """
    
    nSchedule = None
    with open(os.path.join(model_run_dir_Folsom,'TEMP_LOG.OPT'),'r') as fp:
        for line in fp.readlines():
            
            # the last schedule read in the file is the last (potentially iterative) schedule used
            if line.startswith("OPEN FILE:TargetSchedulesA.npt"):
                print(line)
                
                # extract schedule number from line that looks like this: OPEN FILE:TargetSchedulesA.npt   122.000   39    7
                # where 39 is the schedule  number, 122.0 is the day-of-year
                tokens = line.split() # split on whitespace

                # sometimes the last schedule read can simply reload the original scheule used for iterative simulations
                # (because the different scheudules are technicall not defined after 1Dec). We don't want that one, we
                # want the last one loaded that was used for the schedule period (May through Nov).
                if float(tokens[2]) < 333.0:  
                    
                    # is load day-of-year < ~30 Nov
                    nSchedule = int(tokens[3])
    return nSchedule


def computeAlternative(currentAlternative, computeOptions):
    """
    Entry point for the WAT scripting alternative compute for the CE-QUAL-W2
    Folsom Lake post-processing workflow.

    Orchestrates ten sequential operations:
      1. Discovers and organizes W2 outlet flow/temperature input DSS locations.
      2. Determines the simulation year, alternative name, and W2 model variant
         (base iterative, FixedATSP, or NoBypass) from the simulation name.
      3. Locates the CE-QUAL-W2 Folsom model run directory and output files.
      4. Sums all Folsom outlet flows to produce the Natoma inflow record.
      5. Computes the flow-weighted average Folsom release temperature (FWA2).
      6. Reads the final ATSP schedule number from AutoRunTempLog.opt.
      7. Exports the schedule target temperatures (converted to Celsius if needed)
         and schedule number to the forecast DSS file.
      8. Computes downstream regression temperatures at Watt Avenue and Hazel
         Avenue using calc_downstream_temp_W2.
      9. Writes shutter elevation diagnostics from the W2 structure output file.
     10. Writes reservoir storage diagnostics from the W2 volume output file.

    All DSS records are written twice: once under the scripting alternative F-part
    and once under the original W2 F-part, for compatibility with legacy plotting
    tools that cannot locate scripting-alternative records.

    Inputs:
      currentAlternative -- WAT scripting alternative object providing input/output
                            data locations, compute messages, and DSS path creation
      computeOptions     -- WAT compute options object providing the DSS filename,
                            run time window, run directory, and simulation name

    Output:
      Returns True on successful completion.
      Writes the following DSS records to the forecast DSS file (under both
      scripting and W2 F-parts):
        - W2_Natoma_InflowQ              (summed Folsom outlet flows)
        - W2_Natoma_InflowT              (flow-weighted release temperature)
        - W2_FOLSOM_SCHEDULE_TEMP_FINAL  (ATSP schedule target temperatures)
        - W2_FOLSOM_SCHEDULE_FINAL       (ATSP schedule number, constant 24-step series)
        - W2_DownstreamRegressionWatt    (Watt Avenue regression temperature)
        - W2_DownstreamRegressionHazel   (Hazel Avenue regression temperature)
        - W2_Folsom_Forecast_Shutter_1/2/3 (penstock intake elevations, ft)
        - W2_Folsom_Forecast_BypassFlow  (total bypass flow, CMS)
        - W2_Folsom_Forecast_RevisedPenstock1Flow (Penstock 1 flow, CMS)
        - W2_Folsom_Storage / _lt_52F / _lt_60F   (storage volumes, acre-ft)
    """
    
    # Log the start of script execution in the HEC-WAT compute messages.
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
 
    # Retrieve all input locations configured for this scripting alternative.
    locations_obj = currentAlternative.getInputDataLocations()
    
    # Organize locations into paired flow/temperature groups and return DSS paths.
    locations = DSS_Tools.organizeLocationsPaired(currentAlternative, locations_obj, W2FolsomFlowTempLocs, return_dss_paths=True)
    
    # Echo the discovered DSS paths to the compute log for troubleshooting.
    currentAlternative.addComputeMessage('Found DSS paths:')
    
    # Loop through each paired location set.
    for location in locations:
        
        # Display each DSS pathname that was found.
        for path in location:
            currentAlternative.addComputeMessage(str(path))
    
    # Add a blank line to improve readability of the compute log.
    currentAlternative.addComputeMessage('\n')
    
    # Output DSS file associated with this forecast run.
    dss_file = computeOptions.getDssFilename()
    
    # Runtime window controls the simulation start/end dates.
    rtw = computeOptions.getRunTimeWindow()
    
     # Retrieve DSS filename again for local convenience.
    dss_file = computeOptions.getDssFilename()
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()
    
    # Parse the model year from the HEC time string.
    startyear_str = starttime_str[5:9]
    year = int(startyear_str)

    ###########################################################################
    # Determine forecast alternative name and simulation group
    #
    # This section reconstructs the forecast alternative name from the
    # simulation name and available simulation groups on disk. The result
    # is used later to identify the proper CE-QUAL-W2 model directory.
    ###########################################################################
    
    # Root run directory provided by HEC-WAT.
    study_root = computeOptions.getRunDirectory()  # returns watershed directory
    
    # Split the path and strip off the scripts, alternative-group, and run folder to obtain the WAT directory
    study_main = study_root.split(os.sep)
    study_main = os.sep.join(study_main[:-3])
    
    # Directory containing simulation group definitions.
    study_wat = study_main + os.sep + 'wat' + os.sep + 'simGroups'
    
    # Read all files in the simulation group directory.
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

    ###########################################################################
    # Determine CE-QUAL-W2 model variant
    #
    # Multiple model variants exist. Alternative naming conventions drive
    # selection of the proper model directory.
    ###########################################################################
    
    # Current run directory (not currently used elsewhere).
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
    
     ###########################################################################
    # Construct CE-QUAL-W2 output file locations
    ###########################################################################
    
    # Primary W2 model output directory.
    w2_run_dir = os.path.join(study_main, 'runs', sim_name_underscore, 'CeQual-W2', 'Folsom', model_name)
    
    # Structure output file containing intake elevations and flows.
    str_csv = os.path.join(w2_run_dir,'str_br1.csv')
    
    # Reservoir storage output file.
    vol_file = os.path.join(w2_run_dir,'VOLUME_WB1.OPT')
    
    ###########################################################################
    # Build DSS output paths
    ###########################################################################
    
    # Obtain output locations configured in the scripting alternative.
    outputlocations_obj = currentAlternative.getOutputDataLocations()
    
    # Organize locations according to expected output names.
    outputlocations = DSS_Tools.organizeLocations(currentAlternative, outputlocations_obj, OutputLocs)
    
    # Final DSS output paths.
    outputpaths = []
    
    # Create DSS pathnames for each configured output.
    for ol in outputlocations:
        
        opath = currentAlternative.createOutputTimeSeries(ol)
        
        # Split DSS pathname into components.
        tspath = str(opath).split('/')
        
        # Retrieve F-part
        fpart = tspath[6]
    
        # Remove scripting metadata prefixes if present.
        if '|' in fpart:
            # remove everything before | including |
            tspath[6] = fpart[fpart.find('|')+1:]
        
        # Store normalized pathname.
        outputpaths.append('/'.join(tspath))

    # Retrieve original W2 F-part for plotting compatibility.
    W2_fpart = locations[0][0].split('/')[6]

    ###########################################################################
    # Generate Natoma inflow from Folsom outflows
    ###########################################################################
    
    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpaths[0]))
    
    # Flow DSS records from all Folsom outlets.
    flow_locations = [loc[0] for loc in locations]  
    
    # Sum outlet flows into Natoma inflow.
    DSS_Tools.add_flows(currentAlternative, rtw, flow_locations, dss_file, outputpaths[0], dss_file)

    # copy to W2 alt for plotting
    DSS_Tools.copy_dss_ts(outputpaths[0],new_fpart=W2_fpart,dss_file_path=dss_file)    
    
    ###########################################################################
    # Compute flow-weighted outlet temperature
    ###########################################################################
    
    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpaths[1]))
    
    # Ignore extremely small flows during weighting.
    cfs_limit = 1.0 #float
    
    # Generate representative Natoma inflow temperature.
    tsc_natoma_inflow_temp = flowweightaverage.FWA2(currentAlternative, dss_file, rtw, locations, outputpaths[1], 
                                                    cfs_limit, 10.0, return_tsc=True)
    
     # Duplicate under W2 F-part.
    DSS_Tools.copy_dss_ts(outputpaths[1],new_fpart=W2_fpart,dss_file_path=dss_file)    

    ###########################################################################
    # Schedule processing
    ###########################################################################
    
    # W2 model directory used for schedule lookup.
    model_run_dir_Folsom = w2_run_dir
    
    # Shared forecast resources directory.
    shared_dir = os.path.join(study_main,'shared')
    
    # Forecast DSS containing meteorological inputs.
    forecast_dss = os.path.join(shared_dir,'WTMP_American_forecast.dss')
    
    # read AutoRunTempLog.opt and find last schedule used
    # ----------------------------------------------------------------
    nSchedule = nSchedule_from_AutoRunTempLog(model_run_dir_Folsom)

    # Build schedule DSS pathname dynamically.
    schedule_rec = "/WTMP_American/Folsom Iterative Target/TEMP-WATER//1Day/C:0000"+ "%02i"%nSchedule +"|ScheduleA/"
    
    # Open schedule library.
    dssIn = HecDss.open(os.path.join(shared_dir,"ScheduleA_TT_daily.dss"))
    
    # Read selected schedule.
    tsc_sched = dssIn.get(schedule_rec,True)
    
    # Close input DSS file.
    dssIn.close()

    # Adjust schedule dates to current forecast year.
    tsc_sched = rectify_tsc_dates_to_model_year(tsc_sched,year)
    
    # Convert Fahrenheit schedules to Celsius.
    if tsc_sched.units == 'F':
        
        tc = []
        
        for tf in tsc_sched.values:
            tc.append((tf-32.)*5.0/9.0)
        
        tsc_sched.values = tc
        tsc_sched.units = 'C'    
    
    # Open forecast output DSS.
    dssOut = HecDss.open(dss_file)


    ###########################################################################
    # Export schedule target temperatures
    #
    # These records represent the target release temperatures associated
    # with the final schedule selected by the W2 model.
    ###########################################################################
    
    # write schedule temps to output DSS
    sched_parts = tsc_sched.fullName.split('/')
    
    # Use the configured scripting output record as a template.
    out_parts = outputpaths[2].split('/')
    
    # Preserve the original schedule interval (typically 1Day).
    out_parts[5] = sched_parts[5] # get period from original record
    
    # Build the final output pathname.
    tsc_sched.fullName = '/'.join(out_parts)
    
    # Diagnostic logging for DSS troubleshooting.
    print('writing W2 Schedule Target Temps: ')
    print('    '+dss_file)
    print('    '+tsc_sched.fullName)
    
    # Write schedule temperatures to scripting-alternative F-part.
    #
    # NOTE:
    # Historical comment indicates this record may occasionally fail
    # to appear in DSS. If schedule outputs disappear in the future,
    # begin troubleshooting here.
    
    dssOut.put(tsc_sched)  # TODO:  Why does this never show up in output file???
    
    # Duplicate record under original W2 F-part.
    out_parts[6] = W2_fpart
    tsc_sched.fullName = '/'.join(out_parts)
    dssOut.put(tsc_sched) # write out copy under W2_fpart
    
    ###########################################################################
    # Export schedule number
    #
    # Some reporting tools require the active schedule number rather than
    # the full temperature schedule. A constant daily DSS record is used
    # to expose the selected schedule.
    ###########################################################################
    
    print('writing W2 Schedule Number: '+outputpaths[3])
    print('    '+dss_file)
    print('    '+outputpaths[3])
    
    # Write schedule number under scripting-alternative F-part.
    write_constant_1day_ts(dssOut,outputpaths[3],rtw,nSchedule)
    
    # Build equivalent record under W2 F-part.
    out_parts = outputpaths[3].split('/')
    out_parts[6] = W2_fpart
    
    # Duplicate schedule number for W2 plotting tools.
    write_constant_1day_ts(dssOut,'/'.join(out_parts),rtw,nSchedule)  # write out copy under W2_fpart

    ###########################################################################
    # Downstream temperature regression calculations
    #
    # W2 directly simulates release temperatures but not all downstream
    # forecast locations. Empirical regression equations are used to
    # estimate temperatures at key downstream locations.
    #
    # Current outputs:
    #   - Watt Avenue
    #   - Hazel Avenue
    ###########################################################################
    
    currentAlternative.addComputeMessage("Back-calculating downstream temperature using W2 regressions...")
    
    # Load forecast meteorological and flow inputs.
    #
    # Returns:
    #   doys      = day-of-year values
    #   FaveFlow  = average flow (CMS)
    #   Tair      = air temperature (C)
    
    doys,FaveFlow,Tair = fpp.load_tt_data(forecast_dss, starttime_str, endtime_str) # day-of-year,CMS,C    

    # Convert output temperature series into a tsmath object so it can
    # be resampled to daily values.
    TModel_tsm = tsmath(tsc_natoma_inflow_temp)
    
    # Standardize interval for regression calculations.
    TModel_tsc_daily = DSS_Tools.standardize_interval(TModel_tsm,'1day').getData()
    
    # Compute Watt Avenue regression temperatures.
    dtt,DownstreamTempWatt = fpp.calc_downstream_temp_W2(year,1,doys,Tair,TModel_tsc_daily.values,FaveFlow,lagWatt=True)
    
    # Compute Hazel Avenue regression temperatures.
    dtt,DownstreamTempHazel = fpp.calc_downstream_temp_W2(year,3,doys,Tair,TModel_tsc_daily.values,FaveFlow,lagWatt=True)

    ###########################################################################
    # Export Watt Avenue regression temperatures
    ###########################################################################
    
    # Use configured downstream temperature record as template.
    out_parts = outputpaths[4].split('/')
    
    # Regression outputs are daily.
    out_parts[5] = '1Day'
    
    # Replace location name with Watt regression identifier.
    out_parts[2] = 'W2_DownstreamRegressionWatt'
    
    # Assign DSS pathname.
    TModel_tsc_daily.fullName = '/'.join(out_parts)

    # Log output location.
    print('writing W2 Downstream Regression calc: '+TModel_tsc_daily.fullName)
    
    # Replace temperature values with regression output.
    TModel_tsc_daily.values = DownstreamTempWatt
    
    # Write scripting-alternative version.
    dssOut.put(TModel_tsc_daily)
    
    # Save scripting F-part before modification.
    script_fpart = out_parts[6]
    
    # Duplicate under W2 F-part.
    out_parts[6] = W2_fpart
    TModel_tsc_daily.fullName = '/'.join(out_parts)
    
    dssOut.put(TModel_tsc_daily) # write out copy under W2_fpart

    ###########################################################################
    # Export Hazel Avenue regression temperatures
    ###########################################################################
    
    # Restore scripting F-part.
    out_parts[6] = script_fpart
    
    # Replace location with Hazel regression identifier.
    out_parts[2] = 'W2_DownstreamRegressionHazel'
    
    # Build DSS pathname.
    TModel_tsc_daily.fullName = '/'.join(out_parts)
    print('writing W2 Downstream Regression calc: '+TModel_tsc_daily.fullName)
    
     # Replace values with Hazel regression output.
    TModel_tsc_daily.values = DownstreamTempHazel
    
    # Write scripting-alternative version.
    dssOut.put(TModel_tsc_daily)
    
    # Duplicate under W2 F-part.
    out_parts[6] = W2_fpart
    TModel_tsc_daily.fullName = '/'.join(out_parts)
    
    dssOut.put(TModel_tsc_daily) # write out copy under W2_fpart

    ###########################################################################
    # Export shutter elevations and storage diagnostics
    #
    # These outputs support:
    #   - Buzz plots
    #   - Intake elevation diagnostics
    #   - Bypass flow verification
    #   - Reservoir storage analysis
    #
    # Data are extracted directly from CE-QUAL-W2 output files rather
    # than existing DSS records.
    ###########################################################################
    
    # Read a representative DSS record so timing information can be
    # reused when constructing shutter elevation outputs.
    tsc = dssOut.get(str(locations[0][0]),True)
    
    # Simple diagnostic confirming the record was retrieved.
    print('DEBUG:','got data!!!')
    
    # Close output DSS before auxiliary write routine reopens it.
    dssOut.close()   

    # Process W2 structure and storage files and write diagnostics.
    write_shutter_elevations_to_output_dss(str_csv,vol_file,dss_file,tsc)


    ###########################################################################
    # Successful completion
    #
    # At this point the script has:
    #
    #   1. Generated Natoma inflow.
    #   2. Computed flow-weighted release temperature.
    #   3. Exported ATSP schedule temperatures.
    #   4. Exported schedule number.
    #   5. Generated downstream temperature regressions.
    #   6. Exported shutter elevations.
    #   7. Exported storage diagnostics.
    #
    # Returning True signals success to the HEC-WAT scripting framework.
    ###########################################################################
    
    return True

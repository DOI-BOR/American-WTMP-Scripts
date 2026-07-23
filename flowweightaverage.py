from hec.heclib.dss import HecDss
import hec.hecmath.TimeSeriesMath as tsmath
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE
import math,sys,datetime
import DSS_Tools
reload(DSS_Tools)

import tz_offset
reload(tz_offset)

def organizeLocations(currentAlternative, locations):
    """
    Organizes a flat list of WAT data location objects into a list of
    [flow_path, temperature_path] pairs by resolving and correcting
    the DSS pathname for each location.

    Locations are expected to be ordered as alternating flow/temperature
    pairs (index 0 = flow, index 1 = temperature, index 2 = flow, etc.).
    If the total number of locations is odd, the function logs an error
    and calls sys.exit(1).

    DSS pathnames are resolved using currentAlternative.loadTimeSeries()
    and corrected via DSS_Tools.fixInputLocationFpart() to match the
    current alternative's input F-part context.

    Inputs:
      currentAlternative -- WAT scripting alternative object for logging
                            and resolving linked time-series paths
      locations          -- list of WAT DataLocation objects in alternating
                            flow/temperature order; length must be even

    Output:
      Returns a list of 2-element lists, each of the form
      [flow_dss_path, temp_dss_path], where each element is a corrected
      DSS pathname string. Length equals len(locations) // 2.
      Calls sys.exit(1) if the number of locations is odd.
    """
    
    # Initialize list that will store flow/temperature DSS path pairs
    locations_list = []
    
    # Verify that locations are provided as flow/temp pairs
    if len(locations) % 2 != 0:
        currentAlternative.addComputeMessage("Uneven amount of Flow/Temp pairings. Check inputs.")
        sys.exit(1)
    
    # Process each location in sequence
    for li, location in enumerate(locations):
        
        # Load and standardize the DSS path for the current location
        tspath =str(currentAlternative.loadTimeSeries(location))
        tspath = DSS_Tools.fixInputLocationFpart(currentAlternative, tspath)
        
        # Even indices represent the start of a flow/temp pair
        if li % 2 == 0: 
            current_pair = [tspath]
            
        else:
            
            # Odd indices complete the pair and add it to the output list
            current_pair.append(tspath)
            locations_list.append(current_pair)
            
    # Return organized list of DSS path pairs        
    return locations_list


def flow_in_cfs(units,flows):
    """
    Converts a list of flow values to cubic feet per second (CFS), applying a
    unit conversion from cubic meters per second (CMS) if necessary.

    Inputs:
      units -- string; the units label of the source flow record (e.g., 'cfs', 'cms');
               comparison is case-insensitive
      flows -- list of floats; flow values in the units specified by the units argument

    Output:
      Returns a list of floats containing flow values in CFS.
      Calls sys.exit(-1) if the units string is not recognized ('cfs' or 'cms').
    """
    
    # No conversion needed when flow data is already in cfs
    if units.lower()=='cfs':
        return flows
        
    # Convert cms values to cfs    
    elif units.lower()=='cms':
        values_converted = []
        
        # Apply conversion factor to each flow value
        for f in flows:
            values_converted.append(f * 35.314666213)
            
        return values_converted
        
    else:
        # Abort if flow units are not recognized
        print('FWA2: flow units not known:',units)
        sys.exit(-1)

def temperature_in_C(units,temps):
    """
    Converts a list of temperature values to degrees Celsius (C), applying a
    Fahrenheit-to-Celsius conversion if necessary.

    Inputs:
      units -- string; the units label of the source temperature record
               (e.g., 'C', 'deg C', 'F', 'deg F'); comparison is case-insensitive
      temps -- list of floats; temperature values in the units specified by the
               units argument

    Output:
      Returns a list of floats containing temperature values in degrees Celsius.
      Calls sys.exit(-1) if the units string is not recognized.
    """
    
    # No conversion needed when temperatures are already Celsius
    if units.lower()=='c' or units.lower()=='deg c':
        return temps
        
    # Convert Fahrenheit temperatures to Celsius    
    elif units.lower()=='f' or units.lower()=='deg f':
        values_converted = []
        # Apply Fahrenheit-to-Celsius conversion
        for t in temps:
            values_converted.append((t - 32.0)*5.0/9.0)
            
        return values_converted
        
    else:
        # Abort if temperature units are not recognized
        print('FWA2: temperature units not known:',units)
        sys.exit(-1)

def FWA2(currentAlt, dssFile, timewindow, DSSPaths_list, outputname, cfs_limit=None, bad_data_fill_tempC=None, 
         last_override=False,return_tsc=False):
    """
    Computes a flow-weighted average temperature (FWAT) from multiple DSS
    flow/temperature record pairs and writes the result to a DSS record.

    For each time step, valid flow/temperature pairs are identified by applying
    the following filters:
      - Neither flow nor temperature is NaN
      - Flow is greater than cfs_limit (default 0.0) and less than 9e6 CFS
      - Temperature is between 0.0 and 80.0 deg C

    The flow-weighted average is computed as:
        FWAT[i] = sum(flow[i] * temp[i]) / sum(flow[i])

    When no valid pairs exist at a time step, the fill value (bad_data_fill_tempC
    or UNDEFINED_DOUBLE) is used. Optionally, the last DSS pair's temperature can
    override the weighted average when its flow and temperature pass the same filters.

    This is an improved replacement for FWA, which was producing unexpected results.

    Inputs:
      currentAlt          -- WAT scripting alternative object for logging
      dssFile             -- full path to the DSS file for reading and writing
      timewindow          -- WAT run time window object providing start/end time strings
      DSSPaths_list       -- list of 2-element lists, each [flow_dss_path, temp_dss_path];
                             flow units converted from CMS to CFS automatically;
                             temperature units converted from F to C automatically
      outputname          -- DSS pathname string for the output FWAT record
      cfs_limit           -- float; minimum flow threshold in CFS below which a pair is
                             excluded from the weighted average (default None = 0.0)
      bad_data_fill_tempC -- float; temperature value used when no valid pairs exist at
                             a time step (default None = UNDEFINED_DOUBLE)
      last_override       -- bool; if True, the last DSS pair's temperature overrides the
                             weighted average when it passes the flow/temperature filters
                             (default False)
      return_tsc          -- bool; if True, returns the populated TimeSeriesContainer
                             instead of 0 (default False)

    Output:
      Writes the FWAT time series to dssFile at outputname.
      Returns the populated TimeSeriesContainer if return_tsc is True.
      Returns 0 on successful completion if return_tsc is False.
      Calls sys.exit(-1) if flow and temperature record lengths do not match.
    """
    
    # Flow-weighted average temperature calculation with improved validation
    '''
    Made a new flow-weighted average temperature function; other one was producing weirdness 
    '''
    
    # Extract DSS time window strings
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    
    # Log the processing period
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    
    # Open DSS file for reading and writing
    dssFm = HecDss.open(dssFile)

    # Initialize accumulators for weighted-average calculation
    flow_total = []
    flowtemp_total = []
    n_pairs = []

    # Configure optional flow threshold and fill value
    flow_limit = 0.0 if cfs_limit is None else cfs_limit
    fill_value = UNDEFINED_DOUBLE if bad_data_fill_tempC is None else bad_data_fill_tempC
    
    # Process each flow/temperature DSS path pair
    for dspi, dsspaths in enumerate(DSSPaths_list):
        
        # Extract flow and temperature DSS paths
        flow_dss_path = dsspaths[0]
        temp_dss_path = dsspaths[1]
        
        # Log and read flow data
        currentAlt.addComputeMessage(str(flow_dss_path))
        print('FWA2 Reading:',flow_dss_path)
        
        tsc_flow = dssFm.read(flow_dss_path, starttime_str, endtime_str, False).getData()
        
        # Convert flow values to cfs if necessary
        flows = flow_in_cfs(tsc_flow.units,tsc_flow.values)
        
        # Read associated temperature data
        print('FWA2 Reading:',temp_dss_path)
        
        # Track whether this record can be used for override logic
        last_rec_valid = False  # test to see if we can use override values
        
        try:
            
            # Read and standardize temperature values
            tsc_temp = dssFm.read(temp_dss_path, starttime_str, endtime_str, False).getData()
            temps = temperature_in_C(tsc_temp.units,tsc_temp.values)
            
            # Diagnostic output of first records
            print('tscf',tsc_flow.values[0])
            print(flows[0])
            print('tsct',tsc_temp.values[0])
            print(temps[0])
    
            # Use metadata from the first valid temperature record
            if dspi==0:
                nrecs = len(flows)
                temp_type = tsc_temp.type
    
            # Ensure flow and temperature records have matching lengths
            if len(flows) != nrecs or len(temps) != nrecs:
                currentAlt.addComputeMessage("FWA2: record lengths do not match!")
                print("FWA2: record lengths do not match!",nrecs,len(flows),len(temps))
                sys.exit(-1)
    
            # Process each timestep in the record
            for i in range(nrecs):
                
                # Initialize accumulators during first DSS pair
                if dspi==0:
                    n_pairs.append(0)  
                    flow_total.append(0.0)
                    flowtemp_total.append(0.0)
                    
                if not math.isnan(flows[i]) and not math.isnan(temps[i]):
                    
                    # Reject very small or unrealistic flow values
                    if flows[i] > flow_limit and flows[i] < 9.0e6: 
                        
                        # Reject unrealistic temperature values
                        if temps[i] >= 0.0 and temps[i] <= 80.0:
                            # passed the data checks
                            
                            # Passed all validation checks
                            n_pairs[i] += 1
                            flow_total[i] += flows[i]
                            flowtemp_total[i] += flows[i]*temps[i]
    
            # Mark record as successfully processed
            last_rec_valid = True
            
        except:
            
            # Skip records that cannot be read or processed
            currentAlt.addComputeMessage('FWA2: data not addeded for record: '+temp_dss_path)
            last_rec_valid = False
    
    # Store final flow-weighted average temperature values    
    fwat = []
    
    # Diagnostic output showing record count
    print('nrecs:',nrecs)
    
    # Calculate weighted-average temperature for each timestep
    for i in range(nrecs):
        
        # Compute flow-weighted average when valid pairs exist
        if n_pairs[i] > 0:
            fwat.append(flowtemp_total[i]/flow_total[i])        
        else:
            
            # Use configured fill value when no valid data are available
            fwat.append(fill_value)
            
        # Optionally override weighted average with last valid record    
        if last_override and last_rec_valid:
            
            # Apply same validation criteria used during accumulation
            if flows[i] > flow_limit and flows[i] < 9.0e6: 
                
                # Ensure override temperature is within valid range
                if temps[i] >= 0.0 and temps[i] <= 80.0:
                    fwat[i] = temps[i]
        
    # Reuse the last temperature container for DSS output
    tsc_temp.type = temp_type
    
    # Assign output DSS pathname
    tsc_temp.fullName = outputname
    
    # Replace values with calculated flow-weighted temperatures
    tsc_temp.values = fwat
    
    # Write results back to the DSS file
    dssFm.write(tsc_temp)
    
    # Close DSS file and release resources
    dssFm.close()
    
    # Optionally return the populated time-series container
    if return_tsc:
        return tsc_temp
        
    else:
        
        # Return success code when container is not requested
        return 0

def FWA(currentAlt, dssFile, timewindow, DSSPaths_list, outputname, cfs_limit=None):
    """
    Computes a flow-weighted average temperature (FWAT) from multiple DSS
    flow/temperature record pairs and writes the result to a DSS record.

    For each DSS pair, flow*temperature products are computed at each time step.
    Flows below cfs_limit are set to zero before weighting. The final FWAT is:
        FWAT[i] = sum(flow[i] * temp[i]) / sum(flow[i])

    NaN flow-temperature products from any pair are skipped during summation.
    When total flow at a time step is zero, UNDEFINED_DOUBLE is written.

    NOTE: This function was producing unexpected results and has been superseded
    by FWA2. It is retained for reference and backward compatibility.

    Inputs:
      currentAlt    -- WAT scripting alternative object for logging
      dssFile       -- full path to the DSS file for reading and writing
      timewindow    -- WAT run time window object providing start/end time strings
      DSSPaths_list -- list of 2-element lists, each [flow_dss_path, temp_dss_path];
                       CMS flows are converted to CFS automatically
      outputname    -- DSS pathname string for the output FWAT record
      cfs_limit     -- float; flows below this threshold are zeroed out before weighting
                       (default None = no threshold applied)

    Output:
      Writes the FWAT time series to dssFile at outputname.
      Returns 0 on successful completion.
      Logs the number of values written to the compute message log.
    """
    
    # Extract DSS time window strings for reading data
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    
    # Log the requested processing period
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    
    # Open DSS file and initialize storage dictionary
    dssFm = HecDss.open(dssFile)
    dss_data = {}
    
    # Read each flow/temperature DSS path pair
    for dspi, dsspaths in enumerate(DSSPaths_list):
        
        # Retrieve flow DSS pathname
        flow_dss_path = dsspaths[0]
        currentAlt.addComputeMessage(str(flow_dss_path))
        
        # Read flow time series for the specified time window
        flowTS = dssFm.read(flow_dss_path, starttime_str, endtime_str, False)
        
        
        # Extract underlying time-series container
        flowTS = flowTS.getData()
        
        # Save time stamps for later output construction
        hecstarttimes = flowTS.times
        
        # Determine flow units and convert if necessary
        flow_units = flowTS.units
        
        if flow_units.lower() == 'cms':
            
            # Convert cubic meters per second to cubic feet per second
            currentAlt.addComputeMessage('Converting cms to cfs')
            flowvals = []
            
            for flow in flowTS.values:
                flowvals.append(flow * 35.314666213)
                
            # Initialize storage for this DSS pair    
            dss_data[dspi] = {'flow': flowvals} #start the dict
        
        else:
            
            # Store flow values directly when already in cfs
            dss_data[dspi] = {'flow': flowTS.values} #start the dict

        # Optionally remove flows below threshold
        if cfs_limit != None:
            
            for fi, flow in enumerate(dss_data[dspi]['flow']):
                
                if flow < cfs_limit:
                    
#                   # Replace low-flow values with zero
                    dss_data[dspi]['flow'][fi] = 0
                
        # Read corresponding temperature DSS record
        temp_dss_path = dsspaths[1]
        currentAlt.addComputeMessage(str(temp_dss_path))
       
        TempTS = dssFm.read(temp_dss_path, starttime_str, endtime_str, False)
        
        # Extract temperature data container
        TempTS = TempTS.getData()
        
        # Save metadata needed for output record creation
        tempunits = TempTS.units
        temptype = TempTS.type
        
        # Store temperature values for later calculations
        dss_data[dspi]['temp'] = TempTS.values


    # Compute flow x temperature products for each DSS pair   
    for dspi in dss_data.keys():
        
        flowtemps = []
        offset = 0
        
        # Calculate weighted temperature contribution at each timestep
        for i, flow in enumerate(dss_data[dspi]['flow']):
            
            temp = dss_data[dspi]['temp'][i]
            flowtemps.append(flow*temp)
        
        # Store flow-temperature products        
        dss_data[dspi]['flowtemp'] = flowtemps
         
    # Compute total flow across all DSS paths (sum across multiple datasets) 
    total_flows = []
    
    # Use the first DSS path as a reference index base
    dspi = dss_data.keys()[0]
    
    # Loop through each time step
    for i, flow in enumerate(dss_data[dspi]['flow']):
        
        # start with base dataset flow
        temptotalflow = flow
        
        # Sum flow from all other datasets at same index
        for key in dss_data.keys():
            if key != dspi:
                temptotalflow += dss_data[key]['flow'][i]
        total_flows.append(temptotalflow)
    
    # Compute total flow-weighted temperature numerator (flow * temp summed)
    # Handles NaN values carefully
    total_flowtemp = []
    
    dspi = dss_data.keys()[0]
    
    for i, flowtemp in enumerate(dss_data[dspi]['flowtemp']):
        temptotalflowtemp = flowtemp
        
        for key in dss_data.keys():
            if key != dspi:   
                
                # skip invalid temperature values
                if not math.isnan(dss_data[key]['flowtemp'][i]):
                    
                    # if current accumulator is NaN, replace instead of add
                    if math.isnan(temptotalflowtemp):
                        temptotalflowtemp = dss_data[key]['flowtemp'][i]
                    else:
                        temptotalflowtemp += dss_data[key]['flowtemp'][i]
        total_flowtemp.append(temptotalflowtemp)
    
    # Compute flow-weighted average temperature
    # FW_Avg = sum(flow*temp) / sum(flow)
    FW_Avg_vals = []
    
    for i, flow in enumerate(total_flows):
        flowtemp = total_flowtemp[i]
        
         # avoid divide-by-zero
        if flow == 0:
            FW_Avg_vals.append(UNDEFINED_DOUBLE)
            
        else:
            FW_Avg_vals.append(flowtemp / flow)
    
    # Write results to DSS TimeSeriesContainer
    tsc = TimeSeriesContainer()
    tsc.times = hecstarttimes
    tsc.fullName = outputname
    tsc.values = FW_Avg_vals
    tsc.startTime = hecstarttimes[0]
    tsc.units = tempunits
    tsc.type = temptype
    tsc.endTime = hecstarttimes[-1]
    tsc.numberValues = len(FW_Avg_vals)
    tsc.startHecTime = timewindow.getStartTime()
    tsc.endHecTime = timewindow.getEndTime()
    
    # write output time series to DSS file
    dssFm.write(tsc)
    dssFm.close()
    
    # log output size
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(FW_Avg_vals)))
    return 0


# DAILY FLOW-WEIGHTED AVERAGE TEMPERATURE FUNCTION
def FWA_Daily(currentAlt, dssFile, timewindow, DSSPaths_list, outputname, cfs_limit=None,delay_days=0):
    """
    Computes a daily flow-weighted average temperature (FWAT) from multiple DSS
    sub-daily flow/temperature record pairs and writes a daily output record to DSS.

    For each DSS pair, sub-daily flow and temperature values are grouped into 24-hour
    windows anchored at the HEC '24' hour boundary. Within each daily window, the
    flow-weighted temperature contribution is computed as:
        FWTemp_daily = sum(flow[i] * temp[i]) / sum(flow[i])

    Daily contributions from all DSS pairs are then combined:
        FWAT_daily[d] = sum(FWTemp[d] * DailyFlow[d]) / sum(DailyFlow[d])

    When the total daily flow sum is less than 24.0 CFS, UNDEFINED_DOUBLE is written
    for that day to avoid meaningless results from very low-flow periods.

    An optional start delay (delay_days) advances the read window start to account
    for routing travel time.

    NOTE: This function has known issues and has been superseded by FWA2_Daily.
    It is retained for reference.

    Inputs:
      currentAlt    -- WAT scripting alternative object for logging
      dssFile       -- full path to the DSS file for reading and writing
      timewindow    -- WAT run time window object providing start/end time strings
      DSSPaths_list -- list of 2-element lists, each [flow_dss_path, temp_dss_path];
                       CMS flows are converted to CFS automatically
      outputname    -- DSS pathname string for the daily output FWAT record
      cfs_limit     -- float; flows below this threshold are zeroed out before
                       weighting (default None = no threshold applied)
      delay_days    -- int; number of days to advance the start of the read window
                       (default 0)

    Output:
      Writes the daily FWAT time series (deg C) to dssFile at outputname.
      Returns 0 on successful completion.
      Logs the number of values written to the compute message log.
    """
    
    # Convert start time string and optionally apply delay
    starttime_str = timewindow.getStartTimeString()
    
    # Apply day delay if specified
    if delay_days > 0:
        dt_start = DSS_Tools.hec_str_time_to_dt(starttime_str) + datetime.timedelta(days=delay_days)
        starttime_str = dt_start.strftime('%d%b%Y %H%M')      
    
    endtime_str = timewindow.getEndTimeString()
    
    # Log computation time window
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    
    # Open DSS file for reading
    dssFm = HecDss.open(dssFile)
    
    # Dictionary to store flow and temperature datasets per path
    dss_data = {}
    
    # Loop through each DSS path pair (flow, temp)
    for dspi, dsspaths in enumerate(DSSPaths_list):
        
        flow_dss_path = dsspaths[0]
        currentAlt.addComputeMessage(str(flow_dss_path))
        
        # Read flow time series
        flowTS = dssFm.read(flow_dss_path, starttime_str, endtime_str, False)
        
        # Convert timestamps to readable form (for daily grouping)
        readabledates = []
        
        for i in range(len(flowTS.getData().times)):
            hecdate = flowTS.getData().getHecTime(i).dateAndTime(4) #delete later
            readabledates.append(hecdate)
        
        # Extract raw flow time series data
        flowTS = flowTS.getData()
        hecstarttimes = flowTS.times
        
        # Identify "day start" indices (based on HEC time format ending in 24)
        daystart_idx = [ti for ti, t in enumerate(readabledates) if t.split(':')[0][-2:] == '24']
        
        #  Extract DSS integer times marking daily boundaries
        daily_times = [hecstarttimes[t] for t in daystart_idx]
        daily_times_readable = [readabledates[t] for t in daystart_idx]
        
        # FLOW PROCESSING Section
        flow_units = flowTS.units
        
        # Convert CMS to CFS if needed
        if flow_units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')
            flowvals = []
            
            for flow in flowTS.values:
                flowvals.append(flow * 35.314666213)
                
             # Store converted flow values
            dss_data[dspi] = {'flow': flowvals} 
            
        else:
            
             # Store raw flow values if already in correct units
            dss_data[dspi] = {'flow': flowTS.values} 
        
        # Apply optional flow cutoff (filtering out low flows)
        if cfs_limit != None:
            for fi, flow in enumerate(dss_data[dspi]['flow']):
                if flow < cfs_limit:
#                   print('Flow of {0} removed for being under limit: {1}'.format(flow, cfs_limit)) #noisy..
                    dss_data[dspi]['flow'][fi] = 0

        # TEMPERATURE PROCESSING Section
        temp_dss_path = dsspaths[1]
        currentAlt.addComputeMessage(str(temp_dss_path))
        
        # Read temperature series
        TempTS = dssFm.read(temp_dss_path, starttime_str, endtime_str, False)
        TempTS = TempTS.getData()
        
        tempunits = TempTS.units
        dss_data[dspi]['temp'] = TempTS.values
        
        # Compute pointwise flow x temperature products
        flowtemps = []
        total_flows = []
        offset = 0
        
        for i, flow in enumerate(dss_data[dspi]['flow']): #thats everyday bro            
            temp = dss_data[dspi]['temp'][i]
            
            # debug output for first few values
            if i < 10:
                print(i,'Flow-temp pair:',flow,temp)
            
            flowtemps.append(flow*temp)
        
        # DAILY AGGREGATION of Flow-Weighted Temps
        FWTemps_daily = []
        flowsums_daily = []
        
        for ti in daystart_idx:
            
           # approximate 24-hour block; ensure non-negative start index
           dei = ti - 23
           if dei < 0:
                dei = 0
           
           # Sum flows over daily window
           daily_flowsums = sum(dss_data[dspi]['flow'][dei:ti+1])
           flowsums_daily.append(daily_flowsums)
           
           flowtemp_div_flows = []
           
           # Normalize flow-temp products by daily flow sum
           for flowtemp in flowtemps[dei:ti]:
                
                if daily_flowsums == 0:
                    flowtemp_div_flows.append(MISSING_DOUBLE)
                
                else:
                    flowtemp_div_flows.append(flowtemp/daily_flowsums)
           
           # Store daily weighted temperature
           FWTemps_daily.append(sum(flowtemp_div_flows))
           
        # Save intermediate results for dataset
        dss_data[dspi]['FWTemps_daily'] = FWTemps_daily
        dss_data[dspi]['DailyFlowsums'] = flowsums_daily

    # SUM ALL DAILY FLOWS ACROSS DATASETS
    total_Daily_flowsums = []
    dspi = dss_data.keys()[0]
    dailyflowsums = dss_data[dspi]['DailyFlowsums']
    
    for i, dfs in enumerate(dailyflowsums):
        temptotalflow = dfs
        for key in dss_data.keys():
            if key != dspi:
                temptotalflow += dss_data[key]['DailyFlowsums'][i]
        total_Daily_flowsums.append(temptotalflow)

    # COMBINE FLOW-WEIGHTED TEMPERATURE TERMS
    total_FWA_FlowTemps = []
    
    # store intermediate results per dataset
    dspi = dss_data.keys()[0]
    FWTemps = dss_data[dspi]['FWTemps_daily']
    
    for i, fwtd in enumerate(FWTemps):
        
        # Base dataset contribution
        dflow_sum = dss_data[dspi]['DailyFlowsums'][i]
        FW_Flowtemps = fwtd * dflow_sum
        
        temp_Fwflowtemps_summed = FW_Flowtemps           
        
        # Add contributions from all other datasets
        for key in dss_data.keys():
            if key != dspi:
                dflow_sum = dss_data[key]['FWTemps_daily'][i]
                dFWTemp = dss_data[key]['DailyFlowsums'][i]
                temp_Fwflowtemps_summed += dflow_sum * dFWTemp
        
        total_FWA_FlowTemps.append(temp_Fwflowtemps_summed)
      
    # FINAL FLOW-WEIGHTED AVERAGE TEMPERATURE      
    FW_Avg_vals = []
    
    print('len total_FWA_FlowTemps:',len(total_FWA_FlowTemps))
    print('len total_Daily_flowsums:',len(total_Daily_flowsums))
    
    for i, tfwaft in enumerate(total_FWA_FlowTemps):
        
        total_Daily_flowsum = total_Daily_flowsums[i]
        
        # Debug first few values
        if i < 10:
            print(i,'total_Daily_flowsum:',total_Daily_flowsum)
            print(i,'tfwaft:',tfwaft)
        
        # Avoid invalid or too-small totals
        if total_Daily_flowsum < 24.0:
            FW_Avg_vals.append(UNDEFINED_DOUBLE)
        
        else:
            FW_Avg_vals.append(tfwaft/total_Daily_flowsum)
    
    # WRITE OUTPUT DAILY TIME SERIES to DSS
    tsc = TimeSeriesContainer()
    tsc.times = daily_times
    tsc.fullName = outputname
    tsc.values = FW_Avg_vals
    tsc.startTime = daily_times[0]
    tsc.units = tempunits
    tsc.endTime = daily_times[-1]
    tsc.numberValues = len(FW_Avg_vals)

    # Write time series and close file
    dssFm.write(tsc)
    dssFm.close()
    
    # Log output size
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(FW_Avg_vals)))
    
    return 0

# SIMPLE TEMPERATURE CONVERSION UTILITY
def F_to_C(t,is_in_F):
    """
    Converts a single temperature value from Fahrenheit to Celsius if needed.

    Inputs:
      t       -- float; temperature value to convert
      is_in_F -- bool; if True, applies the Fahrenheit-to-Celsius conversion;
                  if False, returns the value unchanged

    Output:
      Returns a float representing the temperature in degrees Celsius.
    """
    
    # Convert Fahrenheit to Celsius if needed
    if is_in_F:
        return (t-32.)*5./9.
    
    else:
        return t

def FWA2_Daily(currentAlt, dssFile, timewindow, DSSPaths_list, outputname, cfs_limit=-1,delay_days=0,delay_hours=0):
    """
    Computes a daily flow-weighted average temperature (FWAT) from multiple DSS
    sub-daily flow/temperature record pairs and writes a 1-day output record to DSS.

    Sub-daily flow and temperature values are grouped into calendar days using
    timezone-corrected Python datetime objects (via tz_offset.timedelta). Within
    each daily bin, valid flow/temperature pairs are accumulated:
      - Flow must be greater than cfs_limit
      - Temperature must be between 0.0 and 120.0 deg C (or F before conversion)

    The daily FWAT is computed as:
        FWAT[d] = sum(flow[i] * temp_C[i]) / sum(flow[i])   for all i in day d

    When total daily flow is zero, UNDEFINED_DOUBLE is written for that day.

    This function supersedes FWA_Daily, which had known aggregation issues.

    Inputs:
      currentAlt    -- WAT scripting alternative object for logging
      dssFile       -- full path to the DSS file for reading and writing
      timewindow    -- WAT run time window object providing start/end time strings
      DSSPaths_list -- list of 2-element lists, each [flow_dss_path, temp_dss_path];
                       CMS flows are converted to CFS via a scalar multiplier;
                       Fahrenheit temperatures are converted to Celsius via F_to_C
      outputname    -- DSS pathname string for the daily output FWAT record
      cfs_limit     -- float; minimum flow threshold in CFS; time steps with flow
                       at or below this value are excluded from the weighted average
                       (default -1 = include all positive flows)
      delay_days    -- int; number of days to advance the read window start
                       (default 0)
      delay_hours   -- int; additional hours to advance the read window start
                       alongside delay_days (default 0)

    Output:
      Writes the daily FWAT time series (deg C, derived from the 1-day resampled
      template of the first temperature record) to dssFile at outputname.
      Returns 0 on successful completion.
    """
    
    # Get start and end time strings from time window
    starttime_str = timewindow.getStartTimeString()
    
    # Apply optional delay (days + hours) to start time
    if delay_days > 0:
        dt_start = DSS_Tools.hec_str_time_to_dt(starttime_str) + datetime.timedelta(days=delay_days,hours=delay_hours)
        starttime_str = dt_start.strftime('%d%b%Y %H%M')      
    
    # Compute end time string for reading DSS records
    endtime_str = timewindow.getEndTimeString()
    
    # Log computation range and open DSS file
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    
    # Open the DSS file for reading
    dssFm = HecDss.open(dssFile)

    # Loop through each DSS flow/temperature pair
    for dspi, dsspaths in enumerate(DSSPaths_list):
        
        flow_dss_path = dsspaths[0]
        currentAlt.addComputeMessage(str(flow_dss_path))        
        
        # Read flow time series and convert to data object
        tsc_flow = dssFm.read(flow_dss_path, starttime_str, endtime_str, False).getData()
        
        # Store datetime array once (only for first dataset)
        if dspi==0:
            dtt = DSS_Tools.hectime_to_datetime(tsc_flow)
        
        # Convert flow units from CMS to CFS if needed
        to_cfs = 1.0
        if tsc_flow.units.lower() == 'cms':
            to_cfs = 35.314666213

        # Read temperature time series
        temp_dss_path = dsspaths[1]
        currentAlt.addComputeMessage(str(temp_dss_path))
        
        tsc_temp = dssFm.read(temp_dss_path, starttime_str, endtime_str, False).getData()
        
        # Initialize daily arrays only once (first dataset)
        if dspi==0:    
            
            # Convert to daily time series template
            out_tsc = tsmath(tsc_temp).transformTimeSeries('1DAY',"","AVE").getData()
            
            # Number of daily bins
            nDaily = len(out_tsc.values)
            
            # Initialize arrays for daily accumulated flow, flow*temp, and final FWA
            flow = [0.0 for i in range(nDaily)]
            flowtemp = [0.0 for i in range(nDaily)]
            fwa = [UNDEFINED_DOUBLE for i in range(nDaily)]
         
        # Temperature unit conversion flag         
        to_C = False
        if tsc_temp.units.lower() == 'f' or tsc_temp.units.lower() == 'degf':
            to_C = True

        # Iterate over sub-daily values and accumulate into daily bins
        di = 0
        
        # initialize first timestamp with timezone offset
        dt0 = dtt[0] + tz_offset.timedelta
        last_day = dt0.day
        for i,(f,t) in enumerate(zip(tsc_flow.values,tsc_temp.values)):
            
            # Convert timestamp with timezone correction
            dttz = dtt[i] + tz_offset.timedelta
            
             # Debug print for first values
            if i < 50:
                print(i,di,last_day,dttz.strftime('%d%b%Y %H%M'))
            
             # Move to next day when day changes
            if dttz.day != last_day:
                last_day = dttz.day
                di += 1
                
                # Stop if we exceed expected number of daily bins
                if di >= nDaily:
                    break
                    
            # Convert flow to CFS      
            f1 = f*to_cfs
            
            # Debug output for early iterations
            if i < 50:
                print(i,di,last_day,dttz.strftime('%d%b%Y %H%M'),f1,F_to_C(t,to_C))
            
            # Apply filters (flow threshold + valid temperature range)
            if f1 > cfs_limit and t >= 0.0 and t < 120.:  
                
                flow[di] += f1
                flowtemp[di] += f1*F_to_C(t,to_C)

    # Compute final flow-weighted average temperature per day
    for i in range(nDaily):  
        if i < 4:
            print(i,flowtemp[i],flow[i])
        # avoid divide-by-zero
        if flow[i] > 0.0:
            fwa[i] = flowtemp[i]/flow[i]
    
    # Write output DSS time series
    out_tsc.fullName = outputname
    out_tsc.units = 'C'
    out_tsc.values = fwa
    
    dssFm.write(out_tsc)
    dssFm.close()
        
    return 0
           




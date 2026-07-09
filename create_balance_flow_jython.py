

import math
from hec.heclib.dss import HecDss
from hec.hecmath import HecMathException
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
from hec.heclib.util import HecTime
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
import hec.hecmath.TimeSeriesMath as tsmath
from rma.util.RMAConst import MISSING_DOUBLE
import math
import sys
import datetime as dt
import os


from hec.io import DSSIdentifier
from hec.heclib.util import HecTime
from com.rma.io import DssFileManagerImpl
from java.util import TimeZone

def linear_interpolation(x_values, y_values, x):
    
    # Validate that x and y arrays have matching lengths.
    # At least two points are required for interpolation.
    if len(x_values) != len(y_values) or len(x_values) < 2:
        raise ValueError("Input lists must have the same length and contain at least 2 data points.")

    # Search for the interval containing the interpolation point.
    for i in range(1, len(x_values)):
        if x <= x_values[i]:
            
            # Define lower and upper bounding points.
            x0, y0 = x_values[i - 1], y_values[i - 1]
            x1, y1 = x_values[i], y_values[i]

            # Compute interpolated value using a linear relationship.
            y = y0 + (y1 - y0) * (x - x0) / (x1 - x0)

            # Return interpolated result immediately.
            return y

    # Requested x value falls outside the lookup table range.
    raise ValueError("Interpolation point is outside the range of provided data.")

def read_elev_storage_area_file(file_name, res_name):
    # Input file contains elevation, storage, and area information.
    # Units are [ft, acre-ft, acre].
    elevstorarea = {} 
    
    # Store each parameter separately for lookup operations.
    elev = []
    stor = []
    area = []
    
    # Import os module for diagnostics.
    import os
    print('cwd: ' + os.getcwd())
    
    # Natoma file contains only elevation and area values.
    if res_name.lower() == 'natoma':
        with open(file_name, 'r') as fn:
            for line in fn:
                
                # Split CSV line into individual fields.
                sline = line.strip().split(',')
                
                # Store elevation and area values.
                elev.append(float(sline[0]))
                area.append(float(sline[1]))
    else:
        
        # Most reservoirs contain elevation, storage, and area.
        with open(file_name, 'r') as fn:
            for line in fn:
                
                # Parse CSV record.
                sline = line.strip().split(',')
                
                # Populate elevation-storage-area arrays.
                elev.append(float(sline[0]))
                stor.append(float(sline[1]))
                area.append(float(sline[2]))
                
    # Store arrays in a dictionary for downstream use.
    elevstorarea['elev'] = elev
    elevstorarea['stor'] = stor
    elevstorarea['area'] = area
    
    # Return structured lookup data.
    return elevstorarea

def build_conic_storage_array(elev, area, firstStorageValue=0.0):
    '''Find storage of slabs between measurement points on the elevation area curve,
    using a conic estimation.  Adapted from storage.java from HEC ResSim, 2022-06-17'''
    # calculate storage at each elevation using conic formula
    n_measures = len(elev)
    
    # Initialize storage array with starting value.
    storage = []
    storage.append(firstStorageValue)
    
    # Compute cumulative storage between elevation layers.
    for i in range(1, n_measures):
        
        # Layer thickness between elevations.
        h = elev[i] - elev[i-1]
        
        # Apply conic volume equation and accumulate storage.
        storage.append(h/3. * (area[i-1] + area[i] + math.sqrt(area[i-1] * area[i])) + storage[i-1])
    
    # Return cumulative storage lookup table.
    return storage


def conic_storage_interp(interpElev, elev, area, conicStorage, idx):
    '''Find storage between measurement points on the elevation area curve,
    using interpolation between conic layers.  Adapted from storage.java from
    HEC ResSim, 2022-06-17'''
    
    # Compute normalized position within layer.
    h = (interpElev - elev[idx]) / (elev[idx+1] - elev[idx])
    
    # Geometric mean area used by conic interpolation.
    geomMean = math.sqrt(area[idx] * area[idx+1])
    
    # Estimate area at target elevation.
    areaInterp = area[idx] + 2.*(geomMean - area[idx])*h + (area[idx] + area[idx+1] - 2.*geomMean)*h*h
    
    # Compute interpolated storage volume.
    storageInterp = (interpElev - elev[idx])/3. * (area[idx] + areaInterp + math.sqrt(area[idx] * areaInterp)) + conicStorage[idx]
    
    # Return estimated storage.
    return storageInterp


def get_elev_layer_idx(elev, obs_elev, elev_stor_area):
    # find lower bounding index of where elevation lands in elev-stor-area table
    # Returned index represents lower bounding layer.

    idx = UNDEFINED_DOUBLE
    min_val = None
    
    # Search all elevation values.
    for i in range(len(elev)):
        
        # Calculate absolute distance from observation.
        valchk = abs(elev[i]-obs_elev) 
        
        # Handle NaN values.
        if math.isnan(valchk):
            min_val = valchk
            idx = i
            
        # Initialize minimum distance.
        elif min_val == None:
            min_val = valchk
            idx = i
            
        # Update nearest elevation index.
        elif valchk < min_val:
            min_val = valchk
            idx = i
            
    # Adjust index if nearest value lies above observation.
    if idx != UNDEFINED_DOUBLE:
        if elev_stor_area['elev'][idx] > obs_elev: 
            idx -= 1
    else:
        # Return invalid index when no match found.
        idx = -1
    return idx

def get_balance_period(balance_period):
    # Convert DSS interval strings into hours.
    
    # Example: "1Hour" -> 1.0
    if 'hour' in balance_period.lower():
        return float(balance_period.lower().replace('hour', ''))
     
    # Example: "1Day" -> 24.0     
    elif 'day' in balance_period.lower():
        return float(balance_period.lower().replace('day', '')) * 24
        
    # Example: "15Min" -> 0.25
    elif 'min' in balance_period.lower():
        return float(balance_period.lower().replace('min', '')) / 60

def check_dss_intervals(records, balance_period, currentAlt):
    # Verify DSS pathname intervals match computation interval.
    for r in records:
        
        # Stop execution if interval mismatch is detected.
        if balance_period.lower() not in r.lower():
            currentAlt.addComputeMessage('DSS record {0} not matching time interval {1}'.format(r, balance_period))
            sys.exit(-1)


def read_ts_rec_w_optional_fname(dssFm, pathname, starttime_str, endtime_str):
    '''pathname may contain the dss filepath additionally before the dss ts path, separated by '::'
       If so, use that dss file.'''
       
    # Check whether an alternate DSS file was specified.  
    if '::' in pathname:
        
        # Log DSS record being opened.
        print('Splitting and reading:',pathname)
        
        # Separate DSS filename from DSS pathname.
        alt_dss_file,pathname_clean = pathname.split('::')
        
        # Open alternate DSS file.
        dssFmRec = HecDss.open(alt_dss_file)
        
        # Read requested time series over specified time window.
        tsc = dssFmRec.read(pathname_clean, starttime_str, endtime_str, False).getData()
        
        # Close alternate DSS file after reading.
        dssFmRec.close()
    else:
        
        # Read directly from default DSS file.
        tsc = dssFm.read(pathname, starttime_str, endtime_str, False).getData()
        
    # Return populated TimeSeriesContainer.
    return tsc


def read_inflows_outflows(currentAlt, dss_file, inflow_records, outflow_records, starttime_str, endtime_str,
                          starttime_hectime, endtime_hectime):

    # Open DSS file containing reservoir flow records.
    dssFm = HecDss.open(dss_file)

    # Initialize aggregate inflow and outflow arrays.
    inflows = []
    outflows = []
    
    # Store DSS timestamps from first valid record.
    times = []

    # Read inflows
    print('Reading inflows')
    
    # Loop through all inflow DSS records.
    for j, inflow_record in enumerate(inflow_records): #for each of the dss paths in inflow_records
        
        # Store pathname for readability.
        pathname = inflow_record
        
        # Log current inflow being processed.
        print('\nreading: ' + str(pathname))
        try:
            
            # Print requested time window.
            print(starttime_str, endtime_str)
            
            # Print source DSS file.
            print(dss_file)
            
            # Read time series from DSS.
            tsc = read_ts_rec_w_optional_fname(dssFm, pathname, starttime_str, endtime_str)
            
            # Extract values and metadata.
            values = tsc.values
            hectimes = tsc.times
            units = tsc.units

            # Trim leading values outside requested time window.
            if hectimes[0] < starttime_hectime: #if startdate is before the timewindow..
                print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], starttime_hectime))
                
                # Calculate number of records to skip.
                st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
                
                # Trim arrays to requested start time.
                values = values[st_offset:]
                hectimes = hectimes[st_offset:]
                
            # Trim trailing values beyond requested time window.
            if hectimes[-1] > endtime_hectime:
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                
                # Calculate number of excess records.
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
                
                # Remove excess values.
                values = values[:(len(hectimes) - st_offset)]
                hectimes = hectimes[:(len(hectimes) - st_offset)]

        except HecMathException:
            
            # Log DSS read failure and terminate.
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))
            sys.exit(-1)

        # Convert metric flow units to cfs if needed.
        if units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')
            print('Converting inflow to cms to cfs')
            convvals = []
            
            # Apply m3/s-to-ft3/s conversion factor.
            for flow in values:
                convvals.append(flow * 35.314666213)
            values = convvals

        # First inflow initializes accumulator.
        if len(inflows) == 0:
            inflows = values
            
            # Save common DSS timestamps.
            times = hectimes #TODO: check how this handles missing values
        else:
            
            # Add current inflow to cumulative inflow series.
            for vi, v in enumerate(values):
                inflows[vi] += v

    
    
    # Read outflows
    print('Reading outflow records')
    
    # Loop through all outflow DSS records.
    for j, outflow_record in enumerate(outflow_records):  # for each of the dss paths in inflow_records    
        
        # Store pathname for readability.
        pathname = outflow_record
        
        # Log current record being processed.
        currentAlt.addComputeMessage('reading' + str(pathname))
        
        # Read outflow time series.
        tsc = read_ts_rec_w_optional_fname(dssFm, pathname, starttime_str, endtime_str)
        try:
            
            # Extract values and metadata.
            values = tsc.values
            hectimes = tsc.times
            units = tsc.units
            
            # Trim records before simulation window.
            if hectimes[0] < starttime_hectime: #if startdate is before the timewindow..
                print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], endtime_hectime))
                st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
                values = values[st_offset:]
                hectimes = hectimes[st_offset:]
                
            # Trim records after simulation window.    
            if hectimes[-1] > endtime_hectime:
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
                values = values[:(len(hectimes) - st_offset)]
                hectimes = hectimes[:(len(hectimes) - st_offset)]
          
        except HecMathException:
            
            # Log DSS read error and terminate.
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))
            sys.exit(-1)

        # Convert metric flow units if required.
        if units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')
            print('Converting outflow cms to cfs')
            convvals = []
            
            # Convert each flow value.
            for flow in values:
                convvals.append(flow * 35.314666213)
            values = convvals

        # First record initializes outflow accumulator.
        if len(outflows) == 0:
            outflows = values
        else:
            # Sum all outflow components together.
            for vi, v in enumerate(values):
                outflows[vi] += v

    # Close DSS file after all reads are complete.
    dssFm.close()

    # Compute inflow minus outflow for each timestep.
    inflow_outflow = []
    for i in range(len(inflows[1:])):
        inflow_outflow.append(inflows[i+1] - outflows[i+1])
   
    # Result is a period-average flow in CFS.
    currentAlt.addComputeMessage("Len inflow_outflow:"+str(len(inflow_outflow)))
    currentAlt.addComputeMessage("Len times:"+str(len(times)))

    # Return timestamps and net reservoir flow.
    return times,inflow_outflow


def predict_elevation(currentAlt, starttime_str, endtime_str, res_name, inflow_records, outflow_records, starting_elevation,
                         elev_stor_area, dss_file, output_dss_record_name, output_dss_file, shared_dir,
                         use_conic=False, alt_period=None, alt_period_string=None, balance_period_str='1Hour'):
    '''From inflows/outflows, predict hourly elevation, useful for lookback/starting elevation for forecasts starting
    on arbitrary dates during forecast period
    '''
    
    # Convert DSS interval string (e.g., "1Hour") into numeric hours.
    balance_period = get_balance_period(balance_period_str) # convert to (float) hours
    
    # Ensure inflow and outflow records match simulation timestep.
    check_dss_intervals(inflow_records, balance_period_str, currentAlt)
    check_dss_intervals(outflow_records, balance_period_str, currentAlt)
    
    # Convert flow volume from CFS to acre-feet per timestep.    
    cfs_2_acreft = balance_period * 3600. / 43559.9
    
    # Precompute inverse conversion for convenience (not directly used later).
    acreft_2_cfs = 1. / cfs_2_acreft

    # Convert time strings to HEC internal time representation.
    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    
    # Read net inflow minus outflow time series for the reservoir.
    times,inflow_outflow = read_inflows_outflows(currentAlt, dss_file, inflow_records, outflow_records, 
                                                 starttime_str, endtime_str, starttime_hectime, endtime_hectime)

    # Log data length consistency for debugging.
    currentAlt.addComputeMessage("Len inflow_outflow:"+str(len(inflow_outflow)))
    currentAlt.addComputeMessage("Len times:"+str(len(times)))
    
    # Initialize storage using elevation-storage curve.
    storage = linear_interpolation(elev_stor_area['elev'], elev_stor_area['stor'], starting_elevation)
    
    # Store initial condition as first element in storage series.
    storage = [storage,]
    
    # Initialize elevation prediction output array.
    elev_predicted = []
    
    #################################################################### 
    # Forward mass-balance simulation 
    #################################################################### 
    
    # Step through each timestep of net inflow.
    for i in range(len(inflow_outflow)):
        
        # Update storage using net flow (CFS to acre-feet conversion).
        storage.append( storage[-1] + inflow_outflow[i]*cfs_2_acreft )
        
        # Convert updated storage back to elevation using inverse lookup.
        elev_predicted.append( linear_interpolation(elev_stor_area['stor'], elev_stor_area['elev'], storage[-1]) )

    #################################################################### 
    # Write elevation output to DSS 
    ####################################################################

    # Open output DSS file for writing results.
    dssFm_out = HecDss.open(output_dss_file)
    
    # Compute timestep spacing from input time series.
    steptime = times[1]-times[0]
    
    # Create DSS time series container for elevation results.
    tsc = TimeSeriesContainer()
    
    # Set starting time for output series.
    tsc.startTime = times[0] 
    
    # Define DSS interval in minutes.
    tsc.interval = int(balance_period)*60
    
    # Assign DSS pathname for elevation output.
    tsc.fullName = output_dss_record_name
    
    # Include initial elevation plus predicted values.
    tsc.values = [starting_elevation] + elev_predicted
    
    # Define output units and type for DSS compatibility.
    tsc.units = 'ft'
    tsc.type = 'INST-VAL'
    tsc.numberValues = len(tsc.values)
    
    # Write elevation time series to DSS.
    dssFm_out.write(tsc)

    #################################################################### 
    # Write storage output to DSS 
    ####################################################################

    # Modify DSS pathname to represent storage output.
    recparts = output_dss_record_name.split('/')
    recparts[3] = 'STORAGE-PREDICTED'
    
    # Reset time alignment for storage output.
    tsc.startTime = times[0]
    tsc.fullName = '/'.join(recparts)
    
    # Write storage values (initial + predicted storage).
    tsc.values = storage
    tsc.units = 'ac-ft'
    tsc.type = 'PER-CUM'
    
    # Store storage time series in DSS.
    dssFm_out.write(tsc)

    #################################################################### 
    # Optional resampling step 
    ####################################################################

    # If alternate output resolution is requested, transform series.
    if alt_period is not None:
        
        # Only resample if interval differs from base computation.
        if alt_period_string.lower() != balance_period_str.lower():
            
            # Read back written DSS record.
            tsm = dssFm_out.read(output_dss_record_name)
            
            # Convert to requested averaging interval.
            tsm_new_interval = tsm.transformTimeSeries(alt_period_string, "", "AVE")
            
            # Write resampled version back to DSS.
            dssFm_out.write(tsm_new_interval)

    # Close DSS file after writing outputs.
    dssFm_out.close()


def create_balance_flows(currentAlt, timewindow, res_name, inflow_records, outflow_records, stage_record, evap_record,
                         elev_stor_area, dss_file, output_dss_record_name, output_dss_file, shared_dir,
                         storage_dss_record_name='', evap_dss_record_name='',
                         balance_period_str="1HOUR", use_conic=False, write_evap=False, write_storage=False,
                         alt_period=None,alt_period_string=None, lookback_padding=1440):

    ''' 
    Compute reservoir mass-balance residuals using: 
    inflows, outflows, stage changes, and evaporation losses. 
    Outputs a DSS time series of balance flow residuals. 
    ''' 
    
    #################################################################### 
    # Validate input time series consistency 
    #################################################################### 
    
    # Ensure all inflow records align with model timestep.
    check_dss_intervals(inflow_records, balance_period_str, currentAlt)
    
    # Ensure all outflow records align with model timestep.
    check_dss_intervals(outflow_records, balance_period_str, currentAlt)
    
    # Ensure stage and evaporation series also match timestep.
    check_dss_intervals([stage_record, evap_record], balance_period_str, currentAlt)
    
    # Convert DSS interval string into numeric hours.
    balance_period = get_balance_period(balance_period_str) 
    
    # Log selected simulation period for debugging.
    print('balance_period ' + str(balance_period))
    
    #################################################################### 
    # Unit conversion factors 
    ####################################################################
    
    # Convert CFS over timestep to acre-feet.
    cfs_2_acreft = balance_period * 3600. / 43559.9
    
    # Convert acre-feet back to CFS-equivalent.
    acreft_2_cfs = 1. / cfs_2_acreft

    #################################################################### 
    # Extract simulation window 
    ####################################################################
    
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()

    # Convert simulation times to HEC integer time format.
    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    
    # Log simulation window.
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))

    #################################################################### 
    # Compute inflow minus outflow time series 
    ####################################################################
    
    
    times,inflow_outflow = read_inflows_outflows(currentAlt, dss_file, inflow_records, outflow_records, 
                                                 starttime_str, endtime_str, starttime_hectime, endtime_hectime)
    
    # Log data consistency.
    print('len times:',len(times))
    print('len inflow_outflow:',len(inflow_outflow))

    #################################################################### 
    # Read stage (reservoir elevation) 
    ####################################################################

    dssFm = HecDss.open(dss_file)

    # Read stage
    print('Reading stage')
    tsc = read_ts_rec_w_optional_fname(dssFm, stage_record, starttime_str, endtime_str)
    try:
        stage = tsc.values
        hectimes = tsc.times
        
        # Trim stage data before simulation window.
        if hectimes[0] < starttime_hectime: #if startdate is before the timewindow..
            print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], endtime_hectime))
            st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
            stage = stage[st_offset:]
            hectimes = hectimes[st_offset:]
            
        # Trim stage data after simulation window.
        if hectimes[-1] > endtime_hectime:
            print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
            st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
            stage = stage[:(len(hectimes) - st_offset)]
            hectimes = hectimes[:(len(hectimes) - st_offset)]
        print('Number Stage Values: {0}'.format(len(stage)))

        # Convert meters to feet if needed.
        if tsc.units.lower() == 'm':
            currentAlt.addComputeMessage('Converting stage m to ft')
            print('Converting stage cms to cfs')
            convvals = []
            for elev in stage:
                convvals.append(elev * 3.280839895)
            stage = convvals
        
    except HecMathException:
        currentAlt.addComputeMessage('ERROR reading' + str(stage_record))
        sys.exit(-1)

    #################################################################### 
    # Read evaporation time series 
    ####################################################################
    
    print('Reading evap')
    tsc = read_ts_rec_w_optional_fname(dssFm, evap_record, starttime_str, endtime_str)
    try:
        evap = tsc.values
        hectimes = tsc.times
        
        # Trim evap before start time.
        if hectimes[0] < starttime_hectime: #if startdate is before the timewindow..
            print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], endtime_hectime))
            st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
            evap = evap[st_offset:]
            hectimes = hectimes[st_offset:]
            
        # Trim evap after end time
        if hectimes[-1] > endtime_hectime:
            print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
            st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
            evap = evap[:(len(hectimes) - st_offset)]
            hectimes = hectimes[:(len(hectimes) - st_offset)]
        print('Number Evap Values: {0}'.format(len(evap)))
    except HecMathException:
        currentAlt.addComputeMessage('ERROR reading' + str(evap_record))
        sys.exit(-1)

    #################################################################### 
    # Build storage-area relationship (conic method) 
    ####################################################################
    
    conic_storage = build_conic_storage_array(elev_stor_area['elev'], elev_stor_area['area'])

    #################################################################### 
    # Initialize balance calculation arrays 
    ####################################################################
    
    n = len(stage) - 1
    flow_resid = []
    flow_evap = []    
    storage_record = []

    #################################################################### 
    # Main mass balance loop 
    ####################################################################
    
    for k in range(n):
        
        # Reservoir elevation at start and end of timestep.
        stage_start = stage[k]
        stage_end = stage[k+1]
        
        ################################################################ 
        # Compute storage from elevation 
        ################################################################

        if use_conic:
            
            # Identify layer index for interpolation.
            idx1 = get_elev_layer_idx(elev_stor_area['elev'], stage_start, elev_stor_area)
            
            # Convert elevation to storage using conic method.
            storage_start = conic_storage_interp(stage_start, elev_stor_area['elev'], elev_stor_area['area'], conic_storage, idx1)
            idx2 = get_elev_layer_idx(elev_stor_area['elev'], stage_end, elev_stor_area)
            storage_end = conic_storage_interp(stage_end, elev_stor_area['elev'], elev_stor_area['area'], conic_storage, idx2)
        else:
            
            # Use simple linear storage interpolation.
            storage_start = linear_interpolation(elev_stor_area['elev'], elev_stor_area['stor'], stage_start)
            storage_end = linear_interpolation(elev_stor_area['elev'], elev_stor_area['stor'], stage_end)

        ################################################################ 
        # Compute storage change and flux components 
        ################################################################

        # Change in storage over timestep (acre-ft).
        delta_stor_from_stage = storage_end - storage_start  # in acre-ft
        
        # Convert storage change to equivalent flow (CFS).
        delta_stor_flow = delta_stor_from_stage * acreft_2_cfs # in cfs
        
        # Net inflow minus outflow from DSS data.
        inflow_minus_outflow = inflow_outflow[k]  # in cfs
        
        # Average reservoir surface area for evaporation.
        area_avg = 0.5 * (linear_interpolation(elev_stor_area['elev'], elev_stor_area['area'], stage_start) +
                          linear_interpolation(elev_stor_area['elev'], elev_stor_area['area'], stage_end))
        
        # Convert evaporation depth to volumetric flow loss.
        evap_flow_loss = (evap[k] * area_avg) * acreft_2_cfs  # in cfs

        ################################################################ 
        # Mass balance residual 
        ################################################################
        
        # Residual = observed storage change - modeled flux balance.
        resid = delta_stor_flow - (inflow_minus_outflow - evap_flow_loss)
        
        # Store outputs for diagnostics.
        flow_resid.append(resid)
        flow_evap.append(evap_flow_loss)
        storage_record.append(storage_start)

    #################################################################### 
    # Export diagnostic CSV (fallback debugging output) 
    ####################################################################
    
    if True:
        print('Writing to CSV')
        # dump to CSV if DSS is mis-behaving or if flows are -999. etc.
        with open(os.path.join(shared_dir, "{0}_balance_flow.csv".format(res_name)), 'w') as opf:
            opf.write('date, balance_flow [cfs]\n')
            for i in range(len(flow_resid)):
                new_line = ','.join([str(times[i]), str(flow_resid[i]), '\n'])
                opf.write(new_line)
        # pd.DataFrame({'date':pd.to_datetime([tstart + model_time_step*i for i in range(len(flow_resid))]),'balance_flow [cfs]':flow_resid}).to_csv("%s balance flow.csv"%res_name)

    dssFm_out = HecDss.open(output_dss_file)
    
    # Output record

    # sometimes ResSim does not include the start record in period average simulations, so if one flow or elevation data
    # record is missing, the calc can sometimes go way off.  Constrain to realistic values, set invalid to zero.
    # also, recs offset by timezone and/or daily records not at midnight can introduced bad values on the first or last days.
    # So, filter at least the first 24 hours, last 24 hours
    
    check_steps = 1
    
    # Replace NaNs and extreme values with zeros.
    bad_flow_bound = 1.e7
    
    # For hourly models, extend validation window.
    if balance_period_str.lower() == '1hour':
        check_steps = 24
    for i in range(check_steps):
        for idx in [i,-1-i]:
            if math.isnan(flow_resid[idx]) or flow_resid[idx] > bad_flow_bound or flow_resid[idx] < -bad_flow_bound:
                flow_resid[idx] = 0.0
                
    #################################################################### 
    # Write DSS balance flow output 
    ####################################################################

    steptime = times[1]-times[0]
    tsc = TimeSeriesContainer()
    
    # Shift start time backward to ensure alignment with ResSim.
    tsc.startTime = times[0] - steptime
    
    # Convert timestep to DSS minutes.
    tsc.interval = int(balance_period)*60
    
    # Assign DSS pathname for output.
    tsc.fullName = output_dss_record_name
    

    # copy back 1st balance flow record 2 steps, instead of writing from 1st valid balance calc.
    # otherwise, time-averaging the balanece flows later leaves off the 1st time step needed for a ResSim run
    # best we can do I guess to make ResSim computes work
    
    # Duplicate first value to avoid missing initial timestep issues.
    tsc.values = [flow_resid[0],flow_resid[0]] + flow_resid
    tsc.units = 'CFS'
    tsc.type = 'PER-AVER'
    tsc.numberValues = len(tsc.values)
    
    # Write balance residual series to DSS.
    dssFm_out.write(tsc)

    #################################################################### 
    # Optional resampling 
    ####################################################################

    if alt_period is not None:
        if alt_period_string.lower() != balance_period_str.lower():
            tsm = dssFm_out.read(output_dss_record_name)
            tsm_new_interval = tsm.transformTimeSeries(alt_period_string, "", "AVE")
            dssFm_out.write(tsm_new_interval)
     
    #################################################################### 
    # Optional evaporation output 
    ####################################################################
    
    if write_evap:
        tsc = TimeSeriesContainer()
        tsc.times = times
        tsc.fullName = evap_dss_record_name
        tsc.values = flow_evap
        tsc.startTime = times[1]
        tsc.units = 'CFS'
        tsc.type = 'PER-AVER'
        tsc.endTime = times[-1]
        tsc.numberValues = len(flow_resid)
        tsc.startHecTime = timewindow.getStartTime()
        tsc.endHecTime = timewindow.getEndTime()
        dssFm_out.write(tsc)
   
    #################################################################### 
    # Optional storage output 
    ####################################################################
    
    if write_storage:
        tsc = TimeSeriesContainer()
        tsc.times = times
        tsc.fullName = storage_dss_record_name
        tsc.values = storage_record
        tsc.startTime = times[1]
        tsc.units = "AC-FT"
        tsc.type = 'INST-VAL'  # is this right?
        tsc.endTime = times[-1]
        tsc.numberValues = len(flow_resid)
        tsc.startHecTime = timewindow.getStartTime()
        tsc.endHecTime = timewindow.getEndTime()
        dssFm_out.write(tsc)

    # Close DSS files.
    dssFm.close()
    dssFm_out.close()
    
    # Indicate successful completion of balance computation.
    return True

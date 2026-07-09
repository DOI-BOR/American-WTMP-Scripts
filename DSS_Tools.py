from hec.heclib.dss import HecDss
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE
from hec.hecmath import HecMathException
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
import hec.hecmath.TimeSeriesMath as tsmath
from com.rma.model import Project
import os,shutil,copy,sys,math
from java.util import Vector, Date

import datetime
from hec.heclib.util import HecTime


def copy_dss_ts(dss_rec,new_fpart=None,new_dss_rec=None,
                dss_file_path=None,dss_file_handle=None):

    # Validate that a DSS source has been supplied.
    # Either a file path or an already-open DSS handle is required.
    if dss_file_path is None and dss_file_handle is None:
        raise ValueError('copy_dss_rec_to_new_fpart: you must supply either a valid dss_file_path OR dss_file_handle')
        
    # Validate that a destination record definition exists.
    # Either a new F-part or a complete DSS pathname is required.    
    if new_fpart is None and new_dss_rec is None:
        raise ValueError('copy_dss_rec_to_new_fpart: you must supply either a new_fpart OR new_dss_rec')

    # Use existing DSS handle if provided.
    # Otherwise open the DSS file.
    if dss_file_handle is not None:
        dss_fm = dss_file_handle
        
    else:
         dss_fm = HecDss.open(dss_file_path)
         
    # Read the full DSS record into a TimeSeriesContainer.     
    tsc = dss_fm.get(dss_rec,True)

    # Use supplied DSS pathname if provided.
    if new_dss_rec is not None:
        dss_rec_out = new_dss_rec
        
    else:
        # Otherwise replace only the F-part of the pathname.
        rec_parts = tsc.fullName.split('/')
        # F-part is position 6 in DSS pathname structure.
        rec_parts[6] = new_fpart
        # Reassemble pathname.
        dss_rec_out = '/'.join(rec_parts)    

    # Assign new pathname to the time series.
    tsc.fullName = dss_rec_out
    
    # Write copied record back to DSS.
    dss_fm.put(tsc)

    # Close file only if this function opened it.
    if dss_file_handle is None:
        dss_fm.close()

def jday_from_tsc(tsc):
    
    # Convert DSS HEC times into Python datetime objects.
    dtt = hectime_to_datetime(tsc)
    
    # Return decimal day-of-year values for each timestamp.
    return [decimal_doy(dt) for dt in dtt]        


def decimal_doy(dt):
    
    # Extract integer day-of-year from the datetime object.
    doy = dt.timetuple().tm_yday
    
    # Compute fractional portion of the day from the time components.
    fractional_day = (dt.hour / 24.0) + (dt.minute / 1440.0) + (dt.second / 86400.0) + (dt.microsecond / 86400000000.0)
    
    # Return day-of-year including fractional day contribution.
    return doy + fractional_day


def organizeLocations(curAlt, location_objs, loc_names, return_dss_paths=False):
    
    # Initialize output list that will contain ordered locations or DSS paths
    locations_list = []
    
    # Print the total number of available location objects
    print('num_locs:',len(location_objs))
    
    # Process locations in the order specified by loc_names
    for name in loc_names:
        
        # Find the index of the requested location
        i_loc = findLocationOrder(curAlt,location_objs,name)
        
        # Diagnostic output showing the matched location index
        print('name:',name,'i_loc:',i_loc)
        
        # Return DSS path strings instead of location objects if requested
        if return_dss_paths:
            
            # Retrieve the corresponding location object
            lo1 = location_objs[i_loc]
            print(lo1)
            
            # Handle locations linked to a previous model differently
            if lo1.isLinkedToPreviousModel():
                
                # Load and adjust the linked time-series path
                tspath = str(curAlt.loadTimeSeries(lo1))
                tspath = fixInputLocationFpart(curAlt, tspath)
                
            else:
                
                # Use the location's DSS path directly
                tspath = str(lo1.getDssPath())
                
            # Display the final DSS path for debugging    
            print(tspath)
            
            # Add DSS path string to output list
            locations_list.append(tspath)
            
        else:
            
            # Add the location object itself to the output list
            locations_list.append(location_objs[i_loc])
            
    # Return ordered list of location objects or DSS paths        
    return locations_list

def organizeLocationsPaired(curAlt, location_objs, loc_names_paired, return_dss_paths=False):   
    
    # Apply organizeLocations to each group of paired location names
    return [organizeLocations(curAlt, location_objs, pn, return_dss_paths) for pn in loc_names_paired]


def findLocationOrder(curAlt,location_objs,name):
    
    # Search through available location objects for a matching name
    for i,loc in enumerate(location_objs):
        
        # Diagnostic output showing the current location being checked
        print(i,'Checking loc: ',loc.getName())
        
        # Return index when a matching location name is found
        if name == loc.getName():
            return i
    # If we make it here, the requested location was not found
    curAlt.addComputeMessage("Scripting - Location name not found: "+name)
    # Stop execution because the requested location is required
    sys.exit(1)

def first_value(dss_file,dss_rec,start_str=None,end_str=None):
    
    # Open the DSS file for reading
    dssFm = HecDss.open(dss_file)  

    # Read the entire record when no time window is specified    
    if start_str is None and end_str is None:
        tsc = dssFm.get(dss_rec,True)    
    else:
        # Read only data within the requested time range
        tsc = dssFm.read(dss_rec,start_str,end_str,False).getData()
        
    # Close the DSS file to release resources
    dssFm.close()
    
    # Return the first value in the retrieved dataset
    return tsc.values[0]


def first_value(dss_file,dss_rec,start_str=None,end_str=None):
    
    # Open the DSS file for reading
    dssFm = HecDss.open(dss_file)    

    # Read the entire record when no time window is specified    
    if start_str is None and end_str is None:
        tsc = dssFm.get(dss_rec,True)
        
    else:
        
        # Read only data within the requested time range
        tsc = dssFm.read(dss_rec,start_str,end_str,False).getData()
        
    # Close the DSS file to release resources    
    dssFm.close()
    
    # Return the first value in the retrieved dataset
    return tsc.values[0]

def standardize_interval(tsm, interval, makePerAver=True):
    
    # Retrieve the underlying time-series container
    tsc = tsm.getData()
    
    # Convert interval name to HEC-DSS interval in minutes
    if interval.lower()=='1hour':
        intint = 60
        
    elif interval.lower()=='1day':
        intint=1440
        
    elif interval.lower()=='1week':
        intint=10800
        
    else:
        
        # Abort if an unsupported interval is requested
        print('interval not supported:',interval)
        sys.exit(-1)

    # Check whether the time series already uses the target interval
    if tsc.interval != intint:
        
        # Optionally enforce PER-AVER type before resampling
        if makePerAver:
            tsm.setType('PER-AVER')
            
        # Transform the time series to the requested interval using averaging
        return tsm.transformTimeSeries(interval, "", "AVE")
        
    else:
        
        # No transformation needed if interval already matches
        return tsm

def get_sanitized_record_list(dss_file_path):
    '''The DSS library seems to return lists of paths with dates in them (getPathnameList()), and some of those
    dates don't even exist in the file or cannot be read and throw an error.  As of Jan 2024,
    this is an orregular problem, and the manual soluation is throwing away the DSS file but
    in many cases that is problematic. So, here we filter dates and check for duplicates.'''
    
    # Open the DSS file and retrieve all pathname records.
    dss = HecDss.open(dss_file_path)
    recs = dss.getPathnameList()
    dss.close()
    
    # Store unique date-independent pathnames.
    sanitized_recs = []
    
    # Remove D-part dates and eliminate duplicates.
    for r in recs:
        rec_tokens = r.split('/')
        
        # Clear the D-part to normalize records.
        rec_tokens[4] = ''  
        r_sanitized = '/'.join(rec_tokens)
        
        # Only keep unique sanitized paths.
        if not r_sanitized in sanitized_recs:
            sanitized_recs.append(r_sanitized)
            
    return sanitized_recs  


def dss_read_ts_safe(dssFilePath,dssRec,start_date=None,end_date=None,returnTSM=False,debug=False):

    # Open DSS file for reading.
    # This call will raise an exception if the file does not exist.
    dss = HecDss.open(dssFilePath,True)

    if start_date is None and end_date is None:
        
        # Read the entire record using get().
        # This is preferred because read() sometimes truncates data.
        tsc = dss.get(dssRec,True) 
        
        dss.close()
        
        # Optional debug logging.
        if debug:
            print('Reading DSS in script...')
            print('    file: '+dssFilePath)
            print('    record: '+dssRec)
            
        # Return either a TimeSeriesMath object
        # or the raw TimeSeriesContainer.
        if returnTSM:
            return tsmath(tsc)
            
        else:
            return tsc
            
    elif start_date is not None and end_date is not None:
        
        # Read only the requested time window.
        tsm = dss.read(dssRec,start_date,end_date,False) 
        
        dss.close()
        
        # Optional debug logging
        if debug:
            print('Reading DSS in script between '+start_date+' and ',+end_date)
            print('    file: '+dssFilePath)
            print('    record: '+dssRec)
            
        # Return requested object type.    
        if returnTSM:
            return tsm
            
        else:
            return tsm.getData()


def data_from_dss(dss_file,dss_rec,starttime_str, endtime_str):
    
    # Open DSS file for reading.
    dssFm = HecDss.open(dss_file)
    
    # Read either the full record or a specified window.
    if starttime_str is None and endtime_str is None:
        tsc = dssFm.get(dss_rec,True)
        
    else:
        tsc = dssFm.read(dss_rec, starttime_str, endtime_str, False).getData()
        
    # Close file after reading.
    dssFm.close()
    
    # Return only the value array.
    return tsc.values


def hectime_to_datetime(tsc):

    # Output list of Python datetime objects.
    dtt = []
    
    # Convert each DSS timestamp individually.
    for j in range(tsc.numberValues):
        
        # Assuming hectime can be converted to Java Date or has method to get the equivalent
        # Extract year and adjust for Java Date convention.
        year = tsc.getHecTime(j).year() - 1900
        month = tsc.getHecTime(j).month() - 1
        
        # Create equivalent Java Date object.
        java_date = Date(year, month, tsc.getHecTime(j).day(), tsc.getHecTime(j).hour(), tsc.getHecTime(j).minute())
        
        # Convert milliseconds since epoch to Python datetime.
        timestamp = (java_date.getTime() / 1000)
        dtt.append(datetime.datetime.fromtimestamp(timestamp))

    return dtt
    
def hectime_to_julian(time_string):

    # Convert the HEC Time to a Java time
    java_date = Date(time_string.year() - 1900, time_string.month() - 1, time_string.day(), time_string.hour(), time_string.minute())
    
    # Convert from the Java date into the python date
    python_date = datetime.datetime.fromtimestamp(java_date.getTime() / 1000)
    
    # Calculate the day of the year
    doy = decimal_doy(python_date)  
    
    # Return to the calling function
    return doy


def fixInputLocationFpart(currentAlternative, tspath):
    
    # Build replacement F-part prefix from current input F-part.
    new_fpart_start = ':'.join(currentAlternative.getInputFPart().split(':')[:-1])
    
    # Split DSS pathname into parts.
    tspath = tspath.split('/')
    fpart = tspath[6]
    
    # Preserve final suffix from original F-part.
    fpart_split = fpart.split(':')
    new_fpart = new_fpart_start + ':' + fpart_split[-1]
    
    # Replace F-part and rebuild pathname.
    tspath[6] = new_fpart
    tspath = '/'.join(tspath)
    
    return tspath

def appendAPart(current_path, ApartAppend):
    
    # Split DSS pathname into components.
    tspath = tspath.split('/')
    Apart = tspath[1]
    
    # If A-part is empty, use supplied value.
    if len(Apart) == 0:
        new_Apart = ApartAppend
        
    else:
        
        # Otherwise append with underscore separator.
        new_Apart = Apart + '_' + ApartAppend
    
    # Update pathname and rebuild string.
    tspath[1] = new_Apart
    tspath = '/'.join(tspath)
    
    return tspath

def getDataLocationDSSInfo(location, currentAlternative, computeOptions):
    
    # Check whether the location references a previous model.
    if location.isLinkedToPreviousModel():
        
        # Load DSS pathname from upstream model output.
        tspath = str(currentAlternative.loadTimeSeries(location))
        tspath = fixInputLocationFpart(currentAlternative, tspath)
        
        # Use DSS file associated with current compute options.
        dsspath = computeOptions.getDssFilename()
    
    else:
        
        # Retrieve linked DSS pathname directly.
        tspath = location.getLinkedToLocation().getDssPath()
        
        # Build full path to DSS file.
        rundir = Project.getCurrentProject().getProjectDirectory()
        dsspath = location.getLinkedToLocation().get_dssFile()
        dsspath = os.path.join(rundir, dsspath)
        
    return tspath, dsspath

def strip_templateID_and_rename_records(dssFilePath,currentAlt):

    # make copy of dss file
    shutil.copyfile(dssFilePath,dssFilePath+'.bak')

    # Open DSS file and retrieve all pathname records.
    dss = HecDss.open(dssFilePath)
    rec_names = dss.getPathnameList()
    new_rec_names = Vector()
    
    # Loop through each DSS pathname and update the F-part.
    for i,r in enumerate(rec_names):          
        parts = r.split('/')
        
        # Skip processing if the F-part does not contain a template ID.
        if not '-' in parts[-2]:
            return
            
        # Remove the first four characters from the F-part.
        parts[-2] = parts[-2][4:]
        new_rec_names.add('/'.join(parts))
        
        # Log the pathname transformation.
        currentAlt.addComputeMessage('Fixing path: '+r+' --> '+new_rec_names[-1])
    
    # Rename all records in a single DSS operation.    
    dss.renameRecords(rec_names, new_rec_names)
    
    # Close the DSS file.
    dss.close()

def add_DSS_Data(currentAlt, dssFile, timewindow, input_data, output_path):
    
    # Get simulation time window boundaries.
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    
    # Open DSS file for reading source records.
    dssFm = HecDss.open(dssFile)
    output_data = []
    
    # Read each input DSS record and accumulate values.
    for dsspath in input_data:
        
        print('reading', str(dsspath))
        ts = dssFm.read(dsspath, starttime_str, endtime_str, False)
        ts = ts.getData()
        
        # Extract values and timestamps.
        values = ts.values
        times = ts.times
        
        # Store metadata from the source record.
        units = ts.units
        
        # Capture DSS data type (INST-VAL, PER-AVER, etc.).
        dsstype = ts.type
        
        # Initialize output array or add values into running total.
        if len(output_data) == 0:
            
            output_data = values
            
        else:
            
            for vi, val in enumerate(values):
                output_data[vi] += val
    
    # Create output DSS time-series container.    
    tsc = TimeSeriesContainer()
    
    # Assign timestamps to output record
    tsc.times = times
    
    # Set DSS pathname for output record.
    tsc.fullName = output_path
    
    # Store aggregated values.
    tsc.values = output_data
    tsc.startTime = times[0]
    
    # Assign engineering units.
    tsc.units = units
    
    # Preserve DSS record type.
    tsc.type = dsstype
    tsc.endTime = times[-1]
    
    # Populate DSS container bookkeeping fields.
    tsc.numberValues = len(output_data)
    tsc.startHecTime = timewindow.getStartTime()
    tsc.endHecTime = timewindow.getEndTime()
    
    # Write output record and close file.
    dssFm.write(tsc)
    dssFm.close()
    
    # Log number of values written.
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(output_data)))
    return 0

def resample_dss_ts(inputDSSFile, inputRec, timewindow, outputDSSFile, newPeriod, 
                    pad_1mon=False,pad_value=None):
    '''Can upsample an even period DSS timeseries, e.g. go from 1DAY -> 1HOUR, or downsample.  However, hecmath likes to
    clip of days that don't have the complete 24 hour cycle.  So, we pad here, but there is a chance we ask for data not
    available. The read gives garbage data and doesn't complain.  
    TODO: figure out how to check for bounds for non-midnight start and end times.
    '''
    
    dssFm = HecDss.open(inputDSSFile)
    
    # Read only the specified time window when provided.
    if timewindow is not None:
        starttime_str = timewindow.getStartTimeString()
        endtime_str = timewindow.getEndTimeString()
        
        # Force full-day boundaries for resampling.
        starttime_str = starttime_str[:-4] + '0000'
        endtime_str = endtime_str[:-4] + '2400' # clipped days don't work in computes ... hope the downloaded DMS data is long enough to do this.
        
        # Optionally extend read period by one month on both ends.
        if pad_1mon:
            dt_start = hec_str_time_to_dt(starttime_str) - datetime.timedelta(days=31)
            starttime_str = dt_start.strftime('%d%b%Y %H%M')      
            dt_end = hec_str_time_to_dt(endtime_str) + datetime.timedelta(days=31)
            endtime_str = dt_end.strftime('%d%b%Y %H%M')
        
        print('Resampling',newPeriod, inputRec,starttime_str,endtime_str)
        tsm = dssFm.read(inputRec, starttime_str, endtime_str, False)
    
    else:
        
        # Read entire record when no time window is supplied.
        print('Resampling',newPeriod, inputRec)
        tsm = dssFm.read(inputRec)  # caution - 'read' sometimes doesn't get whole record?  Need to use get?

    # Optionally pad series with synthetic values at both ends.
    if pad_value is not None:
        tsc = tsm.getData()
        t_diff = tsc.times[1] - tsc.times[0]
        
        new_times = []
        new_values = []
        
        # Copy original data into temporary lists.
        for i in range(tsc.numberValues):
            new_times.append(tsc.times[i])
            new_values.append(tsc.values[i])
        
        # Add one value before and after the record.
        tsc.times = [new_times[0]-t_diff] + new_times + [new_times[-1]+t_diff]
        tsc.values = [pad_value] + new_values + [pad_value]
        
        # Update DSS container metadata.
        tsc.numberValues = tsc.numberValues + 2
        tsc.startTime = tsc.times[0]
        tsc.endTime = tsc.times[-1]
        
        # Perform time-series transformation.
        tsm_new = tsmath(tsc).transformTimeSeries(newPeriod,"","AVE")  
    else:
        
         # Resample directly without padding.
        tsm_new = tsm.transformTimeSeries(newPeriod,"","AVE")
    
    # Close input DSS file.
    dssFm.close()

    # Write transformed series to output DSS file.
    dssFmout = HecDss.open(outputDSSFile)
    dssFmout.write(tsm_new)
    dssFmout.close()


def airtemp_lapse(dss_file,dss_rec,lapse_in_C,dss_outfile,f_part):
    
    # Read temperature time series from DSS.
    dss = HecDss.open(dss_file)
    tsm = dss.read(dss_rec)
    lapse = lapse_in_C
    
    # Convert lapse value if data are stored in Fahrenheit.
    if 'f' in tsm.getUnits().lower():
        lapse = lapse*9.0/5.0+32.0
        
    # Apply temperature offset to entire series.    
    tsm = tsm.add(lapse)
    tsc = tsm.getData()
    dss.close()

    # Replace F-part of pathname for output record.
    pathparts = dss_rec.split('/')
    pathparts[-2] = f_part
    tsc.fullName = '/'.join(pathparts)
    
    # Write modified record to output DSS file.
    dss_out = HecDss.open(dss_outfile)
    dss_out.write(tsc)
    dss_out.close()

def min_ts(dss_file,dss_rec,min_value,dss_outfile,f_part):
    
    # Retrieve DSS time series.
    dss = HecDss.open(dss_file)
    tsc = dss.get(dss_rec,True)
    dss.close()

    # Enforce a minimum value threshold.
    for vi, v in enumerate(tsc.values):
        tsc.values[vi] = max(v, min_value)

    # Rebuild the DSS pathname with the updated F part to reflect the modification
    pathparts = dss_rec.split('/')
    pathparts[-2] = f_part
    tsc.fullName = '/'.join(pathparts)
    
    # Write modified series to output DSS file.
    dss_out = HecDss.open(dss_outfile)
    dss_out.write(tsc)
    dss_out.close()

def add_flows(currentAlt, timewindow, inflow_records, dss_file, output_dss_record_name, output_dss_file):
     
    # Unit conversion references.
    # cfs_2_acreft = balance_period * 3600. / 43559.9
    # acreft_2_cfs = 1. / cfs_2_acreft
    
    # Determine computation time window.
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    
    # Convert DSS times to HEC integer timestamps.
    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dss_file)

    inflows = []
    times = []

    # Read and accumulate all inflow records.
    print('Reading inflows')
    
    #for each of the dss paths in inflow_records
    for j, inflow_record in enumerate(inflow_records): 
        
        pathname = inflow_record
        currentAlt.addComputeMessage('reading' + str(pathname))
        print('\nreading' + str(pathname))
        
        try:
            print(starttime_str, endtime_str)
            
            # Support alternate DSS file syntax using "file::record".
            if '::' in inflow_record:
                dss_file_alt,inflow_rec_alt = inflow_record.split('::')
                dssFm_alt = HecDss.open(dss_file_alt)
                ts = dssFm_alt.read(inflow_rec_alt, starttime_str, endtime_str, False)
                dssFm_alt.close()
                print(dss_file_alt)
                
            else:
                print(dss_file)
                ts = dssFm.read(pathname, starttime_str, endtime_str, False)
                
            # Extract data and metadata from DSS record.    
            ts_data = ts.getData()
            values = ts_data.values
            hectimes = ts_data.times
            units = ts_data.units
            tstype = ts_data.type
            
            # Trim data before requested time window.
            if hectimes[0] < starttime_hectime: 
                
                #if startdate is before the timewindow..
                print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], starttime_hectime))
                
                st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
                values = values[st_offset:]
                hectimes = hectimes[st_offset:]
            
            # Trim data after requested time window.
            if hectimes[-1] > endtime_hectime:
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
                values = values[:(len(hectimes) - st_offset)]
                hectimes = hectimes[:(len(hectimes) - st_offset)]
                

        except HecMathException:
            
            # Abort computation if DSS read fails.
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))
            sys.exit(-1)

        # Convert CMS flows to CFS if necessary.
        if units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')
            convvals = []
            
            for flow in values:
                convvals.append(flow * 35.314666213)
            
            values = convvals

        # Initialize total flow series or accumulate into it.
        if len(inflows) == 0:
            inflows = values
            times = hectimes 
        
        else:
            
            for vi, v in enumerate(values):
                inflows[vi] += v

    # Build output DSS time-series container.
    tsc = TimeSeriesContainer()
    tsc.times = times
    tsc.fullName = output_dss_record_name
    tsc.values = inflows
    
    # Output flow metadata.
    tsc.units = 'CFS'
    tsc.type = tstype
    tsc.numberValues = len(inflows)

    # Write combined inflow record.
    dssFm_out = HecDss.open(output_dss_file)
    dssFm_out.write(tsc)

    # Close DSS resources.
    dssFm.close()
    dssFm_out.close()


def add_or_subtract_flows(currentAlt, timewindow, inflow_records, dss_file, operation,
                       output_dss_record_name, output_dss_file, multiplier=None,
                       delay_days=0,delay_hours=0):
    # Define a function to add or subtract multiple flow records based on flags
    
    # Retrieve start time as string from the time window
    starttime_str = timewindow.getStartTimeString()
    
    if delay_days > 0:
        
        # Apply optional delay to avoid mismatches when data starts later
        dt_start = hec_str_time_to_dt(starttime_str) + datetime.timedelta(days=delay_days,hours=delay_hours)
        starttime_str = dt_start.strftime('%d%b%Y %H%M')          
    
    # Retrieve end time for the time window
    endtime_str = timewindow.getEndTimeString()
    

    # Convert string times to HEC time integers
    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    currentAlt.addComputeMessage('add_or_subtract_flows - Looking from {0} to {1}'.format(starttime_str, endtime_str))
    
    # Open the main DSS file
    dssFm = HecDss.open(dss_file)

    # Initialize containers for combined flows and timestamps
    inflows = []
    times = []

    # Default multipliers to 1.0 if none provided
    if multiplier is None:
        multiplier = [1.0 for i in range(len(inflow_records))]

    # Notify that inflow reading has started
    print('Reading inflows')
    
    # Loop over each DSS path provided for inflow
    for j, inflow_record in enumerate(inflow_records):
        
        pathname = inflow_record
        currentAlt.addComputeMessage('reading' + str(pathname))
        print('\nreading' + str(pathname))
        
        # Attempt to read the inflow time series
        try:
            print(starttime_str, endtime_str)
            
            # If a different DSS file is referenced, open and read from it
            if '::' in inflow_record:
  
                dss_file_alt,inflow_rec_alt = inflow_record.split('::')
                dssFm_alt = HecDss.open(dss_file_alt)
                ts = dssFm_alt.read(inflow_rec_alt, starttime_str, endtime_str, False)
                dssFm_alt.close()
                print(dss_file_alt)
                
            # Otherwise read from the main DSS file    
            else:                
                ts = dssFm.read(pathname, starttime_str, endtime_str, False)                
                print(dss_file)
            
            # Extract values, times, units, and type from the time series object            
            ts_data = ts.getData()
            values = ts_data.values
            hectimes = ts_data.times
            units = ts_data.units
            tstype = ts_data.type

            # Trim leading data earlier than the window start
            #if startdate is before the timewindow..
            if hectimes[0] < starttime_hectime: 
                
                print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], starttime_hectime))
                st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
                values = values[st_offset:]
                hectimes = hectimes[st_offset:]
            
            # Trim trailing data beyond the window end
            if hectimes[-1] > endtime_hectime:
                
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
                values = values[:(len(hectimes) - st_offset)]
                hectimes = hectimes[:(len(hectimes) - st_offset)]
                
        # Exit if the record cannot be read
        except HecMathException:
            
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))
            sys.exit(-1)

        # Convert metric flows (cms) to US customary (cfs)
        if units.lower() == 'cms':
            
            currentAlt.addComputeMessage('Converting cms to cfs')
            convvals = []
            
            for flow in values:
                
                convvals.append(flow * 35.314666213)
            values = convvals

        # Initialize the inflow accumulator
        if len(inflows) == 0:
            
            if multiplier[j] != 1.0:
                inflows = []
                
                for vi, v in enumerate(values):
                    inflows.append(v * multiplier[j])
                    
            else:
                inflows = values
            times = hectimes #TODO: check how this handles missing values
        
        # For subsequent records, add or subtract according to operation flag
        else:
            
            if operation[j]:
                
                for vi, v in enumerate(values):
                    inflows[vi] += v * multiplier[j]
                    
            else:
                
                for vi, v in enumerate(values):
                    inflows[vi] -= v * multiplier[j]

    # Close the main DSS file
    dssFm.close()
                    
    # Prepare a time series container for output
    tsc = TimeSeriesContainer()
    tsc.times = times
    tsc.fullName = output_dss_record_name
    tsc.values = inflows
    
    #Assign metadata for output time series
    tsc.units = 'CFS'
    tsc.type = tstype
    tsc.numberValues = len(inflows)
    
    # Write the computed record to output DSS
    dssFm_out = HecDss.open(output_dss_file)
    dssFm_out.write(tsc)

    # Close output DSS file
    dssFm_out.close()


def hec_str_time_to_dt(hec_str_time):
    '''
    Convert HEC date time format to python datetime object
    '''
    # Handles special case where time ends in 2400 (rolls to next day)
    dt_format = '%d%b%Y %H%M'
    add_day = False
    
    # If the time ends in '2400', replace it with '0000' so strptime can parse it,
    # and flag that one day needs to be added after parsing
    if hec_str_time.endswith('2400'):
        my_hec_str_time = hec_str_time[:-4] + '0000'
        add_day = True
    
    # Otherwise, use the original string as-is
    else:
        my_hec_str_time = hec_str_time

    # Parse the (possibly modified) HEC string into a datetime object
    dt = datetime.datetime.strptime(my_hec_str_time,dt_format)
    
    # If the original time was 2400, advance the datetime by one full day
    if add_day:
        dt = dt + datetime.timedelta(days=1)
    
    # Return converted datetime object
    return dt


def create_constant_dss_rec(currentAlt, timewindow, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS', fpart='ZEROS'):
    '''Create and write a dss record with a constant in it for the given time windows.
       what={'flow','temp-water'}
       period={'1HOUR','1DAY'}
    '''

    # Creates a synthetic time series (e.g., constant flow or temperature)
    if what.lower()=='flow':
        units = 'cfs'
        parameter = 'flow'
        
    elif what.lower()=='temp-water':
        units = 'C'
        parameter = 'temp-water'
        
    elif what.lower()=='gate':
        units = 'n/a'
        parameter = 'gate'
        
    elif what.lower()=='evap':
        units = 'ft'
        parameter = 'evap'
        
    elif what.lower()=='elev':
        units = 'ft'
        parameter = 'elev'
    
    # Validate the requested type of constant
    else:
        currentAlt.addComputeMessage('create_zero_dss_rec: what not known: %s'%what)
        return False

    if period.lower()=='1hour':
        pass
        
    elif period.lower()=='1day':
        pass
        
    else:
        
        # Ensure the period is recognized
        currentAlt.addComputeMessage('create_zero_dss_rec: period not known: %s'%period)
        return False

    # Format used for DTS strings
    dt_format = '%d%b%Y %H%M'
    
    # Retrieve time window bounds
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()

    # Extend window for padding 
    starttime_dt = hec_str_time_to_dt(starttime_str) - datetime.timedelta(days=1)    
    endtime_dt = hec_str_time_to_dt(endtime_str) + datetime.timedelta(days=1)
    
    # Convert padded times back to strings
    starttime_str_pad = starttime_dt.strftime(dt_format)
    endtime_str_pad = endtime_dt.strftime(dt_format)    
 
    # Log the main time window
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))

    ########################
    # Zero-Flow Time Series
    ########################

    # Generate a constant-value time series
    tsmath_zero_flow_day = tsmath.generateRegularIntervalTimeSeries(
        starttime_str_pad,
        endtime_str_pad,
        period, "0M", constant)
        
    # Assign required DSS metadata       
    tsmath_zero_flow_day.setUnits(units)
    tsmath_zero_flow_day.setType(dss_type)
    tsmath_zero_flow_day.setTimeInterval(period)
    tsmath_zero_flow_day.setLocation(cpart)
    tsmath_zero_flow_day.setParameterPart(parameter)
    tsmath_zero_flow_day.setVersion(fpart)

    #  Write and close the resulting record
    dssFm = HecDss.open(output_dss_file)
    dssFm.write(tsmath_zero_flow_day)
    dssFm.close()

    return True


def calculate_relative_humidity(air_temp, dewpoint_temp):
    """
    Calculate Relative Humidity given the air temperature and dewpoint temperature - August-Roche-Magnus approximation

    :param air_temp: Air Temperature in degrees Celsius
    :param dewpoint_temp: Dew Point Temperature in degrees Celsius
    :return: Relative Humidity in percentage
    """
    
    # Compute the numerator and denominator for the approximatio
    numerator = (112.0 - 0.1 * dewpoint_temp + air_temp)
    denominator = (112.0 + 0.9 * air_temp)
    
    # Compute the exponential term comparing saturation pressures
    exponent = ((17.62 * dewpoint_temp) / (243.12 + dewpoint_temp)) - ((17.62 * air_temp) / (243.12 + air_temp))
    
    # Apply the full relative humidity equation and clamp to valid range
    relative_humidity = 100.0 * (numerator / denominator) * math.exp(exponent)
    
    return max(0.01, min(100.0, relative_humidity))

def calculate_dewpoint(air_temp, relative_humidity):
    """
    Calculate Dew Point Temperature given the air temperature and relative humidity,
    using the algebraic inversion of the simplified August-Roche-Magnus approximation.
    
    Parameters:
        air_temp (float): Air temperature in C.
        relative_humidity (float): Relative Humidity in percent (0-100).
    
    Returns:
        float: Dew Point Temperature in C.
    """
    
    # Compute the gamma term used in the inverted Magnus formula
    gamma = math.log(relative_humidity / 100.0) + (17.62 * air_temp) / (243.12 + air_temp)
    
    # Compute dewpoint using the rearranged formula
    dewpoint = 243.12 * gamma / (17.62 - gamma)
    
    return dewpoint


def relhum_from_at_dp(met_dss_file, at_path, dp_path):
    
    # Open DSS file and read air temperature series
    dss = HecDss.open(met_dss_file)
    tsc = dss.read(at_path).getData()
    
    # Retrieve dewpoint data using helper function
    dp_data = data_from_dss(met_dss_file, dp_path, None, None)
    
    # Loop through values and compute relative humidity for each timestep
    for i in range(tsc.numberValues):
        tsc.values[i] = calculate_relative_humidity(tsc.values[i], dp_data[i])
    
    # Adjust pathname parts to reflect derived record
    parts = tsc.fullName.split('/')
    parts[2] = parts[2][:5]
    parts[3] = 'RELHUM-FROM-AT-DP'
    parts[6] = parts[6] + '-DERIVED'
    
    new_pathname = '/'.join(parts)
    tsc.fullName = new_pathname
    tsc.units = '%'
    
    #  Write the new relative humidity record
    print('writing: ', new_pathname)
    dss.write(tsc)
    dss.close()


def dp_from_at_relhum(met_dss_file, at_path, rh_path):
    
    # Open the DSS file and read air temperature
    dss = HecDss.open(met_dss_file)
    tsc = dss.read(at_path).getData()
    
    # Retrieve relative humidity data
    rh_data = data_from_dss(met_dss_file, rh_path, None, None)
    
    # Compute dewpoint at each timestep
    for i in range(tsc.numberValues):
        tsc.values[i] = calculate_dewpoint(tsc.values[i], rh_data[i])
    
    # Update naming for dewpoint output record
    parts = tsc.fullName.split('/')
    parts[2] = parts[2][:5]
    parts[3] = 'temp-dewpoint'
    parts[6] = parts[6] + '-DERIVED'
    
    new_pathname = '/'.join(parts)
    tsc.fullName = new_pathname
    
    print('writing: ', new_pathname)
    dss.write(tsc)
    dss.close()

def check_start_and_end(values, times, startime, endtime):
    
    # Trim the beginning if earlier than the allowable time range
    if times[0] < startime:  
        # Notify the user that the DSS start date precedes the target window
        print('start date ({0}) from DSS before timewindow ({1})..'.format(times[0], startime))
        
        # Calculate how many time steps to skip to reach the desired start time
        st_offset = (startime - times[0]) / (times[1] - times[0])
        
        # Slice both arrays to drop everything before the target start
        values = values[st_offset:]
        times = times[st_offset:]
    
    # Trim the end if extending beyond the target time window
    if times[-1] > endtime:
        
        # Notify the user that the DSS end date exceeds the target window
        print('end date ({0}) from DSS after timewindow ({1})..'.format(times[-1], endtime))
        
        # Calculate how many time steps from the end fall outside the target window
        st_offset = (times[-1] - endtime) / (times[1] - times[0])
        
        # Slice both arrays to drop everything after the target end
        values = values[:(len(times) - st_offset)]
        times = times[:(len(times) - st_offset)]
    
    # Return the trimmed values and times arrays aligned to the target window
    return values, times


def replace_data(currentAlt, timewindow, pairs, dss_file, dss_outfile, months, standard_interval=None):
    
    # Compute and log start and end times for the window
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    
    # Convert start/end strings to HEC integer time values for numeric comparisons
    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    
    # Log the active time window to the compute message log
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    
    
    for pair in pairs:
        
        # Open DSS file and read base record set
        dssFm = HecDss.open(dss_file)
        currentAlt.addComputeMessage('Replacing data for {0} with {1} during {2}'.format(pair[0], pair[1], months))
        
        # Read the base (target) time series over the specified window
        base = dssFm.read(pair[0], starttime_str, endtime_str, False)
        
        # Optionally resample the base series to a standard time interval
        if standard_interval is not None:
            base = standardize_interval(base,standard_interval)
        
        # Extract raw arrays and metadata from the base series container
        base_data = base.getData()
        base_values = base_data.values
        base_hectimes = base_data.times
        
        # Save metadata for later use
        base_units = base_data.units
        base_interval = base_data.interval
        base_type = base_data.type
        
        # Ensure base series is trimmed to relevant window
        base_values, base_hectimes = check_start_and_end(base_values, base_hectimes, starttime_hectime, endtime_hectime)

        
        # Read alternate source series for replacement
        alt = dssFm.read(pair[1], starttime_str, endtime_str, False)
        
        # Optionally resample the alternate series to the same standard interval
        if standard_interval is not None:
            alt = standardize_interval(alt,standard_interval)
        
        # Extract raw arrays and metadata from the alternate series container
        alt_data = alt.getData()
        alt_values = alt_data.values
        alt_hectimes = alt_data.times
        
        # Confirm alignment of replacement data window
        alt_units = alt_data.units
        alt_interval = alt_data.interval
        
        # Trim the alternate series to the same window as the base series
        alt_values, alt_hectimes = check_start_and_end(alt_values, alt_hectimes, starttime_hectime, endtime_hectime)
        dssFm.close()

        # Safety check: units and time intervals must match
        if base_units != alt_units:
            
            # Mismatched units would corrupt the merged output abort
            currentAlt.addComputeMessage('Units do not match for {0} and {1}, skipping'.format(pair[0], pair[1]))
            dssFm.close()
            sys.exit(1)
        
        if base_interval != alt_interval:
            
            # Mismatched intervals cannot be safely merged abort
            currentAlt.addComputeMessage('Intervals do not match for {0} and {1}, changing interval...'.format(pair[0], pair[1]))
            dssFm.close()
            sys.exit(1)

        # Perform month-filtered replacement
        for i in range(len(base_values)):
            
            # Only overwrite base values for time steps falling in the target months
            if base_data.getHecTime(i).month() in months:
                base_values[i] = alt_values[i]
                
        # Construct modified pathname reflecting merged data
        new_pathname = base_data.fullName.split('/')
        alt_pathname = alt_data.fullName.split('/')
        
        # Update the version/part F field to indicate the data origin
        new_pathname[-2] = 'MergedFrom_{0}'.format(alt_pathname[1])
        new_pathname = '/'.join(new_pathname)
        
        print('writing: ',new_pathname)
        
        # Assemble output container
        tsc = TimeSeriesContainer()
        tsc.times = base_hectimes
        tsc.fullName = new_pathname
        tsc.values = base_values
        tsc.units = base_units
        tsc.type = base_type
        
         # Set the record length so the container knows how many values to write
        tsc.numberValues = len(base_values)

        # Write merged result
        dssFmOut = HecDss.open(dss_outfile)
        dssFmOut.write(tsc)
        dssFm.close()
        dssFmOut.close()

def airtemp_lapse(dss_file,dss_rec,lapse_in_C,dss_outfile,f_part):
    
    # Open DSS file and read the temperature record
    dss = HecDss.open(dss_file)
    tsm = dss.read(dss_rec)
    
    # Convert lapse adjustment if data is in Fahrenheit
    lapse = lapse_in_C
    if 'f' in tsm.getUnits().lower():
        lapse = lapse*9.0/5.0+32.0
        
    # Apply lapse adjustment to the entire time series    
    tsm = tsm.add(lapse)
    tsc = tsm.getData()
    dss.close()

    # Update output pathname with modified F-part
    pathparts = dss_rec.split('/')
    pathparts[-2] = f_part
    tsc.fullName = '/'.join(pathparts)
    
    # Write the modified temperature record to output DSS
    dss_out = HecDss.open(dss_outfile)
    dss_out.write(tsc)
    dss_out.close()

def preprend_first_value_on_ts(dss_file,dss_rec,prepend_n):
    '''Sometimes ResSim needs some lookback values, or whatever
    
    Be careful that the first record is where you want to start - sometimes these things change
    '''
    
    # Open the DSS file and fetch the requested time series container
    dss = HecDss.open(dss_file)
    tsc = dss.get(dss_rec,True)

    # Compute the time interval between the first two timesteps
    time_delta = tsc.times[1] - tsc.times[0]
    
    # Convert the HEC time array into a mutable Python list
    times = [tsc.times[i] for i in range(len(tsc.times))] # convert to list - annoying
    
    # Prepend new timestamps spaced backward using the computed timestep
    tsc.times = [times[0] - time_delta*i for i in range(prepend_n,0,-1)] + times
    tsc.startTime = tsc.times[0]
    
    # Convert values array to list to allow inserting repeated values
    values = [tsc.values[i] for i in range(len(tsc.values))] # convert to list - annoying
    
    # Prepend the first value 'prepend_n' times to extend the time series
    tsc.values = [values[0]]*prepend_n + values
    tsc.numberValues = len(tsc.values)

    # Write the updated time series back into the DSS file
    dss.put(tsc)
    dss.close()
   

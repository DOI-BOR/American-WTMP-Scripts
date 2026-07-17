#version 2.0
#modified 03-28-2023 by Scott Burdick-Yahya
#modifed Dec 2023 by Ben Saenz

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
    """
    Copies a DSS time-series record to a new DSS pathname, either by replacing
    only the F-part of the existing pathname or by assigning a fully specified
    new pathname.

    Accepts either an open DSS file handle or a file path string. If a file path
    is provided, the function opens and closes the file itself. If a handle is
    provided, it is used directly and left open after the call.

    Inputs:
      dss_rec         -- DSS pathname string of the source record to copy
      new_fpart       -- string; replacement F-part for the output record pathname
                         (used when new_dss_rec is not provided)
      new_dss_rec     -- string; complete DSS pathname for the output record
                         (takes precedence over new_fpart when both are provided)
      dss_file_path   -- full path to the DSS file (used when dss_file_handle is None)
      dss_file_handle -- open HecDss file handle (used instead of dss_file_path when provided)

    Output:
      No return value. Writes the copied record to the DSS file under the new pathname.
      Raises ValueError if neither a file path nor a handle is supplied, or if neither
      a new F-part nor a complete new pathname is supplied.
    """                

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
    """
    Converts the HEC timestamps in a TimeSeriesContainer to a list of decimal
    day-of-year (Julian day) values.

    Internally converts HEC integer times to Python datetime objects via
    hectime_to_datetime, then applies decimal_doy to each.

    Inputs:
      tsc -- HEC TimeSeriesContainer object containing a populated time axis

    Output:
      Returns a list of floats, one per time step, representing the decimal
      day-of-year (e.g., noon on Jan 1 = 1.5).
    """
    
    # Convert DSS HEC times into Python datetime objects.
    dtt = hectime_to_datetime(tsc)
    
    # Return decimal day-of-year values for each timestamp.
    return [decimal_doy(dt) for dt in dtt]        


def decimal_doy(dt):
    """
    Computes the decimal day-of-year for a given Python datetime object.

    The integer day-of-year (1 = Jan 1) is combined with the fractional portion
    of the day derived from the hour, minute, second, and microsecond components.

    Inputs:
      dt -- Python datetime object

    Output:
      Returns a float representing the decimal day-of-year
      (e.g., 1.0 = midnight Jan 1, 1.5 = noon Jan 1, 32.0 = midnight Feb 1).
    """
    
    # Extract integer day-of-year from the datetime object.
    doy = dt.timetuple().tm_yday
    
    # Compute fractional portion of the day from the time components.
    fractional_day = (dt.hour / 24.0) + (dt.minute / 1440.0) + (dt.second / 86400.0) + (dt.microsecond / 86400000000.0)
    
    # Return day-of-year including fractional day contribution.
    return doy + fractional_day


def organizeLocations(curAlt, location_objs, loc_names, return_dss_paths=False):
    """
    Returns a list of WAT data location objects (or their DSS path strings)
    ordered according to the sequence of names specified in loc_names.

    For each name in loc_names, the matching location is found in location_objs
    using findLocationOrder. Computation is aborted via sys.exit if any requested
    name is not found.

    When return_dss_paths is True, the function resolves each location to its
    DSS pathname string, handling both model-linked locations (via loadTimeSeries
    and fixInputLocationFpart) and externally linked locations (via getDssPath).

    Inputs:
      curAlt           -- WAT scripting alternative object for logging and resolving
                          linked time-series paths
      location_objs    -- list of WAT DataLocation objects available in this alternative
      loc_names        -- list of strings; the location names to extract, in the desired order
      return_dss_paths -- if True, return DSS pathname strings instead of location objects
                          (default False)

    Output:
      Returns a list of WAT DataLocation objects or DSS pathname strings,
      ordered to match loc_names. Length equals len(loc_names).
      Calls sys.exit(1) if any requested location name is not found.
    """
    
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

# I tihnk this version is old ... keeping around in case behavior is broken
#
#def organizeLocations(curAlt, location_objs, loc_names, return_dss_paths=False):
#    locations_list = []
#    print('num_locs:',len(location_objs))
#    for name in loc_names:
#        i_loc = findLocationOrder(curAlt,location_objs,name)
#        if return_dss_paths:
#            tspath = str(curAlt.loadTimeSeries(location_objs[i_loc]))
#            tspath = fixInputLocationFpart(curAlt, tspath)
#            locations_list.append(tspath)
#        else:
#            locations_list.append(location_objs[i_loc])
#    return locations_list


def organizeLocationsPaired(curAlt, location_objs, loc_names_paired, return_dss_paths=False):   
    """
    Applies organizeLocations to each group of location names in a paired (nested)
    list, returning a list of ordered location lists or DSS path lists.

    Useful when locations must be organized into multiple independent groups
    (e.g., paired inflow/outflow sets) that each need to be ordered separately.

    Inputs:
      curAlt            -- WAT scripting alternative object passed through to organizeLocations
      location_objs     -- list of WAT DataLocation objects available in this alternative
      loc_names_paired  -- list of lists of strings; each inner list is a group of location
                           names to organize independently
      return_dss_paths  -- if True, return DSS pathname strings instead of location objects
                           (default False); passed through to organizeLocations

    Output:
      Returns a list of lists, one per group in loc_names_paired. Each inner list contains
      WAT DataLocation objects or DSS pathname strings ordered to match the corresponding
      inner list in loc_names_paired.
      Calls sys.exit(1) via organizeLocations if any requested location name is not found.
    """
    
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
    """
    Searches a list of WAT DataLocation objects for one matching a given name
    and returns its index.

    Performs a linear scan through location_objs, comparing each location's name
    to the requested name. Logs an error and calls sys.exit(1) if the name is
    not found, because a missing location is a fatal configuration error.

    Inputs:
      curAlt        -- WAT scripting alternative object used for logging error messages
      location_objs -- list of WAT DataLocation objects to search
      name          -- string; the location name to find

    Output:
      Returns the integer index of the matching location in location_objs.
      Calls sys.exit(1) if no matching location is found.
    """
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
    """
    Reads a DSS time-series record and returns its first value.

    When no time window is specified, reads the entire record. When both start_str
    and end_str are provided, reads only the specified time window.

    NOTE: This function is defined twice in the original source file. Both definitions
    are preserved here as-is. The second definition supersedes the first at runtime.

    Inputs:
      dss_file  -- full path to the DSS file containing the record
      dss_rec   -- DSS pathname string of the record to read
      start_str -- start date/time string in HEC format (e.g., '01Jan2014 0000');
                   pass None to read the full record (default None)
      end_str   -- end date/time string in HEC format; pass None to read the full
                   record (default None)

    Output:
      Returns the first float value in the retrieved time-series record.
    """
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
    """
    Resamples a TimeSeriesMath object to a target time interval using time-averaging,
    optionally setting its data type to PER-AVER before resampling.

    If the time series is already at the target interval, it is returned unchanged.
    Supported target intervals are '1hour', '1day', and '1week'.

    Inputs:
      tsm         -- HEC TimeSeriesMath object to standardize
      interval    -- target interval string: '1hour', '1day', or '1week' (case-insensitive)
      makePerAver -- if True and a transform is needed, sets the record type to PER-AVER
                     before resampling (default True)

    Output:
      Returns a TimeSeriesMath object at the target interval.
      Calls sys.exit(-1) if the interval string is not recognized.
    """
    
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
    """
    Returns a de-duplicated list of DSS pathnames from a DSS file with the D-part
    (date token) cleared, working around a known DSS library issue where
    getPathnameList() can return date-specific paths that do not exist or cannot
    be read.

    Strips the D-part (index 4) from each pathname and discards duplicates so that
    each unique time-series record is represented only once in the output list.

    Inputs:
      dss_file_path -- full path to the DSS file to query

    Output:
      Returns a list of unique DSS pathname strings with empty D-parts, suitable
      for safe iteration and reading.
    """
    
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
    """
    Reads a DSS time-series record safely, returning either a TimeSeriesContainer
    or a TimeSeriesMath object depending on the returnTSM flag.

    When no time window is specified, uses dss.get() to read the full record
    (preferred over read() which can truncate data). When a time window is
    provided, uses dss.read() for the specified range.

    Inputs:
      dssFilePath -- full path to the DSS file to read from
      dssRec      -- DSS pathname string of the record to read
      start_date  -- start date/time string in HEC format (e.g., '01Jan2014 0000');
                     pass None to read the full record (default None)
      end_date    -- end date/time string in HEC format; pass None to read the full
                     record (default None)
      returnTSM   -- if True, return a TimeSeriesMath object; if False, return a
                     TimeSeriesContainer (default False)
      debug       -- if True, print file path and record name to stdout (default False)

    Output:
      Returns either a TimeSeriesContainer or a TimeSeriesMath object containing
      the requested record data, depending on returnTSM.
      Raises an exception if the DSS file does not exist.
    """
    
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
    """
    Reads a DSS time-series record and returns its values array.

    When starttime_str and endtime_str are both None, reads the full record using
    dss.get(). Otherwise reads only the specified time window using dss.read().

    Inputs:
      dss_file      -- full path to the DSS file to read from
      dss_rec       -- DSS pathname string of the record to read
      starttime_str -- start date/time string in HEC format (e.g., '01Jan2014 0000'),
                       or None to read the full record
      endtime_str   -- end date/time string in HEC format, or None to read the full record

    Output:
      Returns the values array from the TimeSeriesContainer for the requested record.
    """
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
    """
    Converts the HEC integer timestamps in a TimeSeriesContainer to a list of
    Python datetime objects.

    Constructs an equivalent Java Date object for each time step by extracting
    year, month, day, hour, and minute from the HEC time, then converts the
    Java Date to a Python datetime via a Unix timestamp.

    Inputs:
      tsc -- HEC TimeSeriesContainer object with a populated time axis and
             getHecTime() method returning HEC time objects per index

    Output:
      Returns a list of Python datetime objects, one per time step in tsc.
      Length equals tsc.numberValues.
    """
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
    """
    Converts a single HEC time object to a decimal day-of-year (Julian day) value.

    Constructs a Java Date from the HEC time's year, month, day, hour, and minute
    fields, converts it to a Python datetime, and then applies decimal_doy.

    Inputs:
      time_string -- HEC time object with year(), month(), day(), hour(), and
                     minute() accessor methods (e.g., from HecTime.value())

    Output:
      Returns a float representing the decimal day-of-year corresponding to the
      input HEC time (e.g., 1.5 = noon on Jan 1).
    """
    
    # Convert the HEC Time to a Java time
    java_date = Date(time_string.year() - 1900, time_string.month() - 1, time_string.day(), time_string.hour(), time_string.minute())
    
    # Convert from the Java date into the python date
    python_date = datetime.datetime.fromtimestamp(java_date.getTime() / 1000)
    
    # Calculate the day of the year
    doy = decimal_doy(python_date)  
    
    # Return to the calling function
    return doy


def fixInputLocationFpart(currentAlternative, tspath):
    """
    Updates the F-part of a DSS pathname to match the current alternative's input
    F-part prefix, preserving the last colon-separated token of the original F-part.

    Used to correct linked time-series pathnames when the WAT alternative context
    changes between model runs (e.g., when the run ID embedded in the F-part changes).

    Inputs:
      currentAlternative -- WAT scripting alternative object providing getInputFPart()
      tspath             -- DSS pathname string whose F-part is to be updated

    Output:
      Returns the updated DSS pathname string with the corrected F-part.
    """
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
    """
    Appends a suffix to the A-part of a DSS pathname, separated by an underscore.
    If the A-part is empty, the suffix becomes the entire A-part.

    Inputs:
      current_path  -- DSS pathname string whose A-part is to be modified
                       NOTE: contains a pre-existing bug - uses the undefined variable
                       'tspath' instead of 'current_path' on its first line, which will
                       raise a NameError at runtime
      ApartAppend   -- string to append to the existing A-part

    Output:
      Returns the updated DSS pathname string with the modified A-part.
    """
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
    """
    Resolves the DSS pathname and DSS file path for a given WAT data location,
    handling both model-linked and externally linked locations.

    For model-linked locations, the DSS pathname is resolved from the alternative's
    linked time series and the F-part is corrected via fixInputLocationFpart. The
    file path is taken from the compute options DSS filename. For externally linked
    locations, the DSS path and file path are read directly from the linked location.

    Inputs:
      location           -- WAT DataLocation object describing the data source
      currentAlternative -- WAT scripting alternative object for resolving linked time series
      computeOptions     -- WAT compute options object providing the compute DSS filename

    Output:
      Returns a tuple (tspath, dsspath) where:
        tspath  -- the resolved DSS pathname string for the record
        dsspath -- the full path string to the DSS file containing the record
    """
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
    """
    Strips the first 4 characters (template ID prefix) from the F-part of all DSS
    record pathnames in a file, renaming them in place. A backup copy of the original
    file is created before any changes are made.

    If any record's F-part does not contain a hyphen, the function returns early
    without making any changes, assuming no template ID is present.

    Inputs:
      dssFilePath -- full path to the DSS file whose records are to be renamed
      currentAlt  -- WAT scripting alternative object used for logging renamed paths

    Output:
      No return value. Renames records in-place in the DSS file.
      A backup of the original file is saved at dssFilePath + '.bak'.
    """
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
    """
    Reads multiple DSS time-series records from a single DSS file and writes
    their element-wise sum to a new record in the same file.

    Units and data type are taken from the last record read. All input records
    must share the same time axis for the result to be meaningful.

    Inputs:
      currentAlt   -- WAT scripting alternative object used for logging
      dssFile      -- full path to the DSS file containing both input and output records
      timewindow   -- WAT run time window object providing start/end time strings and HecTime values
      input_data   -- list of DSS pathname strings to read and sum
      output_path  -- DSS pathname string for the output summed record

    Output:
      Returns 0 on successful completion.
      Writes the summed time-series record to dssFile at output_path.
    """
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
    """
    Resamples a DSS time-series record to a new time interval using time-averaging.

    Handles both upsampling (e.g., 1DAY -> 1HOUR) and downsampling. Optionally
    extends the read window by one month on each end (pad_1mon) to prevent the
    HEC math library from clipping incomplete day cycles at the boundaries.
    Optionally pads the source series with a synthetic value at both ends (pad_value)
    before resampling.

    Inputs:
      inputDSSFile -- full path to the source DSS file containing the record to resample
      inputRec     -- DSS pathname string of the record to resample
      timewindow   -- WAT run time window object providing start/end time strings,
                      or None to read the full record
      outputDSSFile-- full path to the destination DSS file for the resampled record
      newPeriod    -- target interval string (e.g., '1HOUR', '1DAY')
      pad_1mon     -- if True, extends the read window by 31 days on each end to avoid
                      boundary clipping (default False)
      pad_value    -- if not None, a float value prepended and appended to the source
                      series before resampling to prevent edge discontinuities (default None)

    Output:
      No return value. Writes the resampled time-series record to outputDSSFile.
    """
    
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
    """
    Applies an air temperature lapse rate correction to a DSS temperature record
    and writes the result with an updated F-part to an output DSS file.

    If the source record's units are in Fahrenheit, the lapse value is converted
    from Celsius to Fahrenheit before being applied.

    NOTE: This function is defined twice in the original source file. Both definitions
    are identical and are preserved here as-is. The second definition supersedes the
    first at runtime.

    Inputs:
      dss_file   -- full path to the source DSS file containing the temperature record
      dss_rec    -- DSS pathname string of the temperature record to adjust
      lapse_in_C -- lapse rate correction value in degrees Celsius (added to each value)
      dss_outfile-- full path to the destination DSS file for the corrected record
      f_part     -- new F-part string to apply to the output DSS pathname

    Output:
      No return value. Writes the lapse-corrected temperature record to dss_outfile.
    """
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
    """
    Enforces a minimum value threshold on a DSS time-series record and writes the
    result with an updated F-part to an output DSS file.

    Unlike min_ts_flow_cfs in earlier modules, this function does not perform unit
    conversion; the min_value is applied directly in the record's native units.

    Inputs:
      dss_file  -- full path to the source DSS file containing the record
      dss_rec   -- DSS pathname string of the record to process
      min_value -- minimum allowable value (in the record's native units)
      dss_outfile -- full path to the destination DSS file for the corrected record
      f_part    -- new F-part string to apply to the output DSS pathname

    Output:
      No return value. Writes the minimum-constrained record to dss_outfile.
    """
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
    """
    Reads multiple DSS flow records and writes their element-wise sum to a new DSS record.

    Supports reading records from alternate DSS files using the '::' separator syntax
    (e.g., 'alt_file.dss::dss_pathname'). Trims values outside the run time window,
    converts CMS to CFS where needed, and writes the summed output in CFS.

    Inputs:
      currentAlt             -- WAT scripting alternative object for logging
      timewindow             -- WAT run time window object providing start/end time strings
      inflow_records         -- list of DSS pathname strings to read and sum; entries may
                                use the 'filepath::dssrecord' syntax for alternate DSS files
      dss_file               -- full path to the primary DSS file for records without '::'
      output_dss_record_name -- DSS pathname string for the output summed record
      output_dss_file        -- full path to the DSS file where the output record is written

    Output:
      No return value. Writes the summed flow record (CFS, preserving source data type)
      to output_dss_file. Calls sys.exit(-1) if any record cannot be read.
    """
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
    """
    Reads multiple DSS flow records and combines them element-wise using per-record
    add or subtract operations, with optional per-record scalar multipliers and an
    optional time delay applied to the start of the read window.

    Supports reading records from alternate DSS files using the '::' separator syntax.
    Trims values outside the run time window, converts CMS to CFS where needed, and
    writes the combined output in CFS.

    Inputs:
      currentAlt             -- WAT scripting alternative object for logging
      timewindow             -- WAT run time window object providing start/end time strings
      inflow_records         -- list of DSS pathname strings to combine; entries may use
                                the 'filepath::dssrecord' syntax for alternate DSS files
      dss_file               -- full path to the primary DSS file for records without '::'
      operation              -- list of booleans/None parallel to inflow_records:
                                  None or True  = add this record to the accumulator
                                  False         = subtract this record from the accumulator
                                  (first record always initializes the accumulator)
      output_dss_record_name -- DSS pathname string for the output record
      output_dss_file        -- full path to the DSS file where the output record is written
      multiplier             -- list of floats parallel to inflow_records; each value is
                                applied to its record before adding/subtracting
                                (default None, treated as all 1.0)
      delay_days             -- integer number of days to advance the start of the read window,
                                used to account for routing travel time (default 0)
      delay_hours            -- integer number of hours to advance the start of the read window
                                in addition to delay_days (default 0)

    Output:
      No return value. Writes the combined flow record (CFS, preserving source data type)
      to output_dss_file. Calls sys.exit(-1) if any record cannot be read.
    """
    
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
    """
    Converts a HEC-formatted date/time string to a Python datetime object.

    Handles the HEC convention of '2400' end-of-day timestamps by converting
    them to midnight of the following day.

    Inputs:
      hec_str_time -- HEC date/time string in the format 'ddMmmYYYY HHMM'
                      (e.g., '01Jan2014 2400')

    Output:
      Returns a Python datetime object representing the equivalent date and time.
    """
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
    """
    Creates and writes a DSS time-series record filled with a constant value
    for the run time window, padded by one day on each end.

    Supports flow, water temperature, gate, evaporation, and elevation parameter types.
    The time window is extended by one day before and after to support lookback and
    balance flow calculations that require data outside the strict compute window.

    Inputs:
      currentAlt      -- WAT scripting alternative object for logging error messages
      timewindow      -- WAT run time window object providing start/end time strings
      output_dss_file -- full path to the DSS file where the constant record is written
      constant        -- constant value to fill the record with (default 0.0)
      what            -- parameter type string: 'flow', 'temp-water', 'gate', 'evap',
                         or 'elev' (default 'flow')
      dss_type        -- DSS data type string (e.g., 'PER-AVER', 'INST-VAL') (default 'PER-AVER')
      period          -- time interval string: '1HOUR' or '1DAY' (default '1HOUR')
      cpart           -- C-part (location) string for the DSS pathname (default 'ZEROS')
      fpart           -- F-part (version) string for the DSS pathname (default 'ZEROS')

    Output:
      Returns True on successful write.
      Returns False if 'what' or 'period' is not recognized (also logs an error message).
    """

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
    """
    Computes relative humidity from air temperature and dewpoint temperature DSS records
    using the August-Roche-Magnus approximation, and writes the derived record back to
    the meteorology DSS file.

    The output pathname is derived from the air temperature record's pathname:
      - B-part (location) is truncated to 5 characters
      - C-part (parameter) is set to 'RELHUM-FROM-AT-DP'
      - F-part (version) has '-DERIVED' appended

    Inputs:
      met_dss_file -- full path to the meteorology DSS file containing both input records
      at_path      -- DSS pathname string for the air temperature record (deg C)
      dp_path      -- DSS pathname string for the dewpoint temperature record (deg C)

    Output:
      No return value. Writes the derived relative humidity record (%, same time axis as
      air temperature) back to met_dss_file.
    """
    
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
    """
    Computes dewpoint temperature from air temperature and relative humidity DSS records
    using the inverted August-Roche-Magnus approximation, and writes the derived record
    back to the meteorology DSS file.

    The output pathname is derived from the air temperature record's pathname:
      - B-part (location) is truncated to 5 characters
      - C-part (parameter) is set to 'temp-dewpoint'
      - F-part (version) has '-DERIVED' appended

    Inputs:
      met_dss_file -- full path to the meteorology DSS file containing both input records
      at_path      -- DSS pathname string for the air temperature record (deg C)
      rh_path      -- DSS pathname string for the relative humidity record (%)

    Output:
      No return value. Writes the derived dewpoint temperature record (deg C, same time
      axis as air temperature) back to met_dss_file.
    """
    
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
    """
    Trims a parallel pair of values and times arrays to fit within a specified
    HEC integer time window, removing leading values before startime and trailing
    values after endtime.

    Assumes uniform time spacing; uses the difference between the first two time
    steps to compute the number of steps to trim.

    Inputs:
      values  -- list or array of data values to trim
      times   -- list or array of HEC integer time values corresponding to values
      startime-- HEC integer time value for the desired start of the window (inclusive)
      endtime -- HEC integer time value for the desired end of the window (inclusive)

    Output:
      Returns a tuple (values, times) with both arrays trimmed to the specified window.
    """
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
    """
    Replaces values in one or more base DSS time-series records with values from
    corresponding alternate records, but only for time steps falling within specified
    calendar months. The merged result is written to an output DSS file under a new
    pathname indicating the data source.

    For each pair, both the base and alternate records are read, optionally resampled
    to a standard interval, and trimmed to the run time window. Values in the base
    record for time steps in the specified months are overwritten with the corresponding
    alternate values. The output pathname's F-part is set to 'MergedFrom_<A-part of alt>'.

    Inputs:
      currentAlt        -- WAT scripting alternative object for logging
      timewindow        -- WAT run time window object providing start/end time strings
      pairs             -- list of 2-element lists [base_path, alt_path]; each pair defines
                           one base DSS record and one replacement source record
      dss_file          -- full path to the DSS file containing both base and alt records
      dss_outfile       -- full path to the DSS file where merged results are written
      months            -- list of integer month numbers (1=January, 12=December) during
                           which base values will be replaced with alternate values
      standard_interval -- optional interval string (e.g., '1hour') to resample both base
                           and alt records before merging (default None = no resampling)

    Output:
      No return value. Writes the merged time-series record for each pair to dss_outfile.
      Calls sys.exit(1) if units or time intervals do not match between base and alt records.
    """
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
            
            # Mismatched units would corrupt the merged output to abort
            currentAlt.addComputeMessage('Units do not match for {0} and {1}, skipping'.format(pair[0], pair[1]))
            dssFm.close()
            sys.exit(1)
        
        if base_interval != alt_interval:
            
            # Mismatched intervals cannot be safely merged - abort
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
    """
    Applies an air temperature lapse rate correction to a DSS temperature record
    and writes the result with an updated F-part to an output DSS file.

    If the source record's units are in Fahrenheit, the lapse value is converted
    from Celsius to Fahrenheit before being applied.

    NOTE: This is the second definition of airtemp_lapse in the original source file.
    Both definitions are identical and are preserved as-is. This definition supersedes
    the first at runtime.

    Inputs:
      dss_file   -- full path to the source DSS file containing the temperature record
      dss_rec    -- DSS pathname string of the temperature record to adjust
      lapse_in_C -- lapse rate correction value in degrees Celsius (added to each value)
      dss_outfile-- full path to the destination DSS file for the corrected record
      f_part     -- new F-part string to apply to the output DSS pathname

    Output:
      No return value. Writes the lapse-corrected temperature record to dss_outfile.
    """
    
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
    """
    Prepends the first value of a DSS time-series record to itself n times,
    extending the series backward in time by n evenly spaced time steps.

    Useful when a model (e.g., ResSim) requires lookback values before the
    first available data point. The time step size is inferred from the difference
    between the first two existing timestamps.

    NOTE: The function name contains a pre-existing typo ('preprend' instead of 'prepend').
    This is preserved as-is from the original source.

    Inputs:
      dss_file  -- full path to the DSS file containing the record to extend
      dss_rec   -- DSS pathname string of the record to prepend values to
      prepend_n -- integer number of additional time steps to prepend at the start

    Output:
      No return value. Writes the extended time-series record back to dss_file
      in-place under the same DSS pathname.
    """
    
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
   

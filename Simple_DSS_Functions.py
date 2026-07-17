from hec.heclib.dss import HecDss
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE
from hec.hecmath import HecMathException
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
import hec.hecmath.TimeSeriesMath as tsmath

def add_DSS_Data(currentAlt, dssFile, timewindow, input_data, output_path):
    """
    Sum multiple DSS time-series records and write the result to a new path.
    Opens a DSS file, reads each time-series path listed in
    ``input_data`` over the given time window, accumulates them
    element-by-element into a single aggregated series, and writes
    the combined result back to the same DSS file under
    ``output_path``.
    Parameters
    ----------
    currentAlt : Alternative object
        The current scripting alternative.  Used for logging via
        ``.addComputeMessage()``.
    dssFile : str
        File path to the HEC-DSS file that contains the input
        records and will receive the output record.
    timewindow : TimeWindow object
        Run time window that defines the period of interest.
        Exposes ``.getStartTimeString()``, ``.getEndTimeString()``,
        ``.getStartTime()``, and ``.getEndTime()``.
    input_data : list of str
        List of fully qualified DSS pathname strings to read and
        sum (e.g. ``['//LOC1/FLOW//1Day/TAG/',
        '//LOC2/FLOW//1Day/TAG/']``).
    output_path : str
        Fully qualified DSS pathname where the aggregated time
        series will be written.
    Returns
    -------
    status : int
        Always 0, indicating completion.
    """
    
    # Get the time window strings for the DSS read operation.
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    
    # Log the time range being processed.
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    
    # Open the DSS file for reading input records.
    dssFm = HecDss.open(dssFile)
    
    # Initialize the aggregated output data container.
    output_data = []
    
    # Loop through each input DSS path and accumulate the time series values.
    for dsspath in input_data:
        print('reading', str(dsspath))
        ts = dssFm.read(dsspath, starttime_str, endtime_str, False)
        ts = ts.getData()
        values = ts.values
        times = ts.times
        units = ts.units
        
        # Use the first series as the base, then add subsequent series element-by-element.
        if len(output_data) == 0:
            output_data = values
        
        else:
            for vi, val in enumerate(values):
                output_data[vi] += val
    
    # Create the output time series container for writing the aggregated result.    
    tsc = TimeSeriesContainer()
    tsc.times = times
    tsc.fullName = output_path
    tsc.values = output_data
    tsc.startTime = times[0]
    tsc.units = units
    tsc.endTime = times[-1]
    tsc.numberValues = len(output_data)
    tsc.startHecTime = timewindow.getStartTime()
    tsc.endHecTime = timewindow.getEndTime()
    
    # Write the combined time series back to DSS and close the file.
    dssFm.write(tsc)
    dssFm.close()
    
    # Log how many values were written to the output record.
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(output_data)))
    
    # Return zero to indicate completion.
    return 0

def resample_dss_ts(inputDSSFile, inputRec, timewindow, outputDSSFile, newPeriod):
    """
    Resample a regular-interval DSS time series to a new period.
    Reads a single DSS record, pads the time window to full-day
    boundaries (midnight-to-midnight) to prevent HecMath from
    clipping incomplete days, transforms the series to the
    requested period using period-average interpolation, and writes
    the resampled result to an output DSS file.
    Supports both upsampling (e.g. 1Day to 1Hour) and downsampling
    (e.g. 1Hour to 1Day).  Because the time window is padded, the
    function assumes the source DSS record contains data that
    covers the padded range; if it does not, the read may return
    undefined values without raising an error.
    Parameters
    ----------
    inputDSSFile : str
        File path to the HEC-DSS file containing the source record.
    inputRec : str
        Fully qualified DSS pathname of the record to resample
        (e.g. ``'//LOC/FLOW//1Day/TAG/'``).
    timewindow : TimeWindow object
        Run time window that defines the period of interest.
        Exposes ``.getStartTimeString()`` and
        ``.getEndTimeString()``.  Both are padded internally to
        midnight boundaries before reading.
    outputDSSFile : str
        File path to the HEC-DSS file where the resampled record
        will be written.  May be the same as ``inputDSSFile``.
    newPeriod : str
        Target time-step string recognized by HecMath's
        ``transformTimeSeries`` (e.g. ``'1Hour'``, ``'1Day'``,
        ``'1Month'``).
    Returns
    -------
    None
        The resampled series is written to ``outputDSSFile`` as a
        side effect.  No explicit value is returned.
    """
    
    # Open the input DSS file and get the original time window.
    dssFm = HecDss.open(inputDSSFile)
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    
    # Pad the time window to full-day boundaries to avoid clipped-day issues during resampling.
    #if newPeriod.lower() == '1day':  # some computes don't end on 2400, causes problems when last day doesn't get produced in this func
    starttime_str = starttime_str[:-4] + '0000'
    endtime_str = endtime_str[:-4] + '2400' # clipped days don't work in computes ... hope the downloaded DMS data is long enough to do this.
    
    # Read the source record and transform it to the requested time period.
    print('Resampling',newPeriod, inputRec,starttime_str,endtime_str)
    tsm = dssFm.read(inputRec, starttime_str, endtime_str, False)
    tsm_new = tsm.transformTimeSeries(newPeriod,"","AVE")
    dssFm.close()

    # Open the output DSS file, write the resampled series, and close the file.
    dssFmout = HecDss.open(outputDSSFile)
    dssFmout.write(tsm_new)
    dssFmout.close()


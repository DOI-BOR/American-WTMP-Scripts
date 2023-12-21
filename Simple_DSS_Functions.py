from hec.heclib.dss import HecDss
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE
from hec.hecmath import HecMathException
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
import hec.hecmath.TimeSeriesMath as tsmath

def add_DSS_Data(currentAlt, dssFile, timewindow, input_data, output_path):
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)
    output_data = []
    for dsspath in input_data:
        print('reading', str(dsspath))
        ts = dssFm.read(dsspath, starttime_str, endtime_str, False)
        ts = ts.getData()
        values = ts.values
        times = ts.times
        units = ts.units
        if len(output_data) == 0:
            output_data = values
        else:
            for vi, val in enumerate(values):
                output_data[vi] += val
                
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
    dssFm.write(tsc)
    dssFm.close()
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(output_data)))
    return 0

def resample_dss_ts(inputDSSFile, inputRec, timewindow, outputDSSFile, newPeriod):
    '''Can upsample an even period DSS timeseries, e.g. go from 1DAY -> 1HOUR, or downsample.  However, hecmath likes to
    clip of days that don't have the complete 24 hour cycle.  So, we pad here, but there is a chance we ask for data not
    available. The read gives garbage data and doens't complain.  
    TODO: fiugre out how to check for bounds for non-midnight start and end times.
    '''
    dssFm = HecDss.open(inputDSSFile)
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    #if newPeriod.lower() == '1day':  # some computes don't end on 2400, causes problems when last day doesn't get produced in this func
    starttime_str = starttime_str[:-4] + '0000'
    endtime_str = endtime_str[:-4] + '2400' # clipped days don't work in computes ... hope the downloaded DMS data is long enough to do this.
    print('Resampling',newPeriod, inputRec,starttime_str,endtime_str)
    tsm = dssFm.read(inputRec, starttime_str, endtime_str, False)
    tsm_new = tsm.transformTimeSeries(newPeriod,"","AVE")
    dssFm.close()

    dssFmout = HecDss.open(outputDSSFile)
    dssFmout.write(tsm_new)
    dssFmout.close()


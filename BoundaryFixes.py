from hec.heclib.dss import HecDss
from hec.io import TimeSeriesContainer

def replaceValuesOverThresh(currentAlt, dssFile, timewindow, primary_data_dsspath, secondary_data_dsspath, tertiary_data_dsspath, threshold, ):
    '''
    When Primary file under threshold, use established data values from tertiary_data_dsspath record
    when over threshold, use values from secondary data dsspath
    '''
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)
    
    PrimaryTS = dssFm.read(primary_data_dsspath, starttime_str, endtime_str)
    PrimaryTS = PrimaryTS.getData()

    Primary_units = PrimaryTS.units
    PrimaryTS_values = PrimaryTS.values
    if Primary_units.lower() == 'cms':
        currentAlt.addComputeMessage('Converting cms to cfs')
        PrimaryTS_values = []
        for val in PrimaryTS.values:
            PrimaryTS_values.append(val * 35.314666213)

    SecondaryTS = dssFm.read(secondary_data_dsspath, starttime_str, endtime_str)
    SecondaryTS = SecondaryTS.getData()
    SecondaryTS_values = SecondaryTS.values

    ExistingTS = dssFm.read(tertiary_data_dsspath, starttime_str, endtime_str)
    ExistingTS = ExistingTS.getData()
    ExistingTS_values = ExistingTS.values
    ExistingTS_times = ExistingTS.times
    ExistingTS_units = ExistingTS.units

    for i, primary_val in enumerate(PrimaryTS_values):
        if primary_val > threshold:
            ExistingTS[i] = SecondaryTS_values[i]

    tsc = TimeSeriesContainer()
    tsc.times = ExistingTS_times
    tsc.fullName = tertiary_data_dsspath
    tsc.values = ExistingTS
    tsc.startTime = ExistingTS_times[0]
    tsc.units = tempunits
    tsc.endTime = ExistingTS_times[-1]
    tsc.numberValues = len(ExistingTS)
    tsc.startHecTime = timewindow.getStartTime()
    tsc.endHecTime = timewindow.getEndTime()
    dssFm.write(tsc)
    dssFm.close()
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(ExistingTS)))
    return 0    
         
         
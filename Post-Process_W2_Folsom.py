import os
from com.rma.io import DssFileManagerImpl
from hec.heclib.dss import HecDss
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from hec.heclib.util import HecTime

import Forecast_preprocess as fpp
reload(fpp)

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

    run_dir = computeOptions.getRunDirectory() # runs/model_alt_name/scripting
    model_name = 'W2 Folsom'  # base Folsom forecast model (iterative w/ bypass)
    if "FixedATSP" in run_dir:
        model_name += " FixedATSP"
    if "NoBypass" in run_dir:
        model_name += " NoBypass"
    model_dir = fpp.model_dir_from_run_dir(run_dir,'Folsom',model_name)
    run_base_dir,_ = os.path.split(run_dir)
    model_run_dir_Folsom = os.path.join(run_base_dir,'CeQual-W2','Folsom',model_name)
    targt_temp_npt_filepath = os.path.join(model_dir,'TargetSchedulesA.npt') # overwrite what's there    
    shared_dir = os.path.join(fpp.study_dir_from_run_dir(run_dir),'shared')
    forecast_dss = os.path.join(shared_dir,'WTMP_American_forecast.dss')
    schedule_csv = os.path.join(model_dir,'SchedulesA.csv')

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
    with open(os.path.join(model_run_dir_Folsom,'TEMP_LOG.OPT'),'r') as fp:
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

    
    dssOut = HecDss.open(dss_file)

    # write schedule number to output DSS
    sched_parts = tsc_sched.fullName.split('/')
    out_parts = outputpaths[0].split('/')
    out_parts[5] = sched_parts[5] # get period from original record
    tsc_sched.fullName = '/'.join(out_parts)
    tsc_sched.numberValues = len(tsc_sched.values)
    #tsc_sched.setStoreAsDoubles(True)
    print('writing W2 Schedule Targt Temps: '+tsc_sched.fullName)
    dssOut.put(tsc_sched)  # TODO:  Why does this never show up in output file???
    
    # use above tsc to write constant nSchedule to dss, I guess
    print('writing W2 Schedule Number: '+outputpaths[1])
    write_constant_1day_ts(dssOut,outputpaths[1],rtw,nSchedule)
    dssOut.close()
    
    return True


def rectify_tsc_dates_to_model_year(tsc,model_year):

    ystr = str(model_year)

    new_hec_times = []
    for j in range(tsc.numberValues):
        # Assuming hectime can be converted to Java Date or has method to get the equivalent
        date_str = tsc.getHecTime(j).dateAndTime()  # 01Jan2010 0000
        date_str = date_str[:5]+ystr+date_str[9:]
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
    

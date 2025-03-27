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

# The goal here is to copy records back to the W2 Folsom model alt F-part, so that plotting
# can be more easily used.  We link a folsom record as the first input, and get the F-part
# from that.  The remaining inputs are copied to the results DSS with the W2 fpart.

W2_Folsom_linked_rec = 'W2_Folsom_link' # only used for get f part, could be any W2 output

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

    # pull W2 loc, organize rest into list
    locations_obj = currentAlternative.getInputDataLocations()
    locations_paths = []
    print('num_locs:',len(locations_obj))
    for loc in locations_obj:
        if loc.getName() == W2_Folsom_linked_rec:
            W2_Folsom_path = str(currentAlternative.loadTimeSeries(loc))
        else:
            locations_paths.append(str(currentAlternative.loadTimeSeries(loc)))
    
    # copy other recs to W2 fpart - there is potential for overwirintg here, but seems
    # unlikely that two models would have the same rec name that you want to copy?
    W2_fpart = W2_Folsom_path.split('/')[6]
    for loc_path in locations_paths:
        DSS_Tools.copy_dss_ts(loc_path,new_fpart=W2_fpart,dss_file_path=dss_file)
    
    return True


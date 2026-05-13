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

### Determine the alternative name ###
    # This step is necessary to set the model up correctly based on American configurations
    # Get the current folder
    study_root = computeOptions.getRunDirectory()  # returns watershed directory
    
    # Split the path and strip off the scripts, alternative-group, and run folder to obtain the WAT directory
    study_wat = study_root.split(os.sep)
    study_wat = os.sep.join(study_wat[:-3])
    
    # Create the directory to the WAT simulation group
    study_wat = study_wat + os.sep + 'wat' + os.sep + 'simGroups'
    
    # Get the contents of the directory
    all_groups = os.listdir(study_wat)
    
    # Filter to only forecast groups
    forecast_groups = [x for x in all_groups if ".fsimgrp" in x]
    
    # Remove backup files not caught by the previous filter
    forecast_groups = [x for x in forecast_groups if '.bak' not in x]
    
    for x in forecast_groups:
        currentAlternative.addComputeMessage("group name: " +  x)
    
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
    currentAlternative.addComputeMessage("final group name: " +  x)
    alt_name = sim_name.split(group_name)[0]
    
    # Remove the final dash that was not removed
    alt_name = alt_name[:-1]
    currentAlternative.addComputeMessage("Alternative name: " +  alt_name)  # log message  

    ### Adjust the model name based on the alternative name ###
    # Get the current run directory
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
    
    ### Setup the paths ###
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


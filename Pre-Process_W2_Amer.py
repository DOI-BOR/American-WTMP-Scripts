import os, sys
from com.rma.model import Project
from hec.heclib.util import HecTime
from java.util import Date

# print current path
print("Current paths: ", sys.path)

# create list of unwanted folders in sys.path
search_list = ["SacTrn", "Sacramento", "American", "Stanislaus"]

# initialize and search for unwanted paths
matching_paths = []
for p in sys.path:
    if any(phrase in p for phrase in search_list):
        matching_paths.append(p)

# print paths containing unwanted phrases
print("Paths to be removed:")
for path in matching_paths:
    print(path)

# remove matching paths from sys.path
for path in matching_paths:
    if path in sys.path:
        sys.path.remove(path)

# append path
sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))

import DSS_Tools
reload(DSS_Tools)

import DMS_preprocess
reload(DMS_preprocess)

import Forecast_preprocess as fpp
reload(fpp)

import create_balance_flow_jython as cbfj
reload(cbfj)

shutter_links = ['init_elev_shutter_1','init_elev_shutter_2','init_elev_shutter_3']

def computeAlternative(currentAlternative, computeOptions):
    """
    Orchestrate the CE-QUAL-W2 Folsom forecast preprocessing pipeline.
    Prepares all input files, DSS records, boundary conditions, and
    restart configurations needed to run a CE-QUAL-W2 temperature
    forecast of Folsom Reservoir on the American River.  Steps include
    resolving the alternative name, fixing DSS metadata, converting
    storage to elevation, building target-temperature schedules,
    generating outlet flow files, writing constant placeholder records,
    deriving initial shutter positions, and updating W2 control files.
    Parameters
    ----------
    currentAlternative : Alternative object
        The scripting alternative being computed.  Exposes
        ``.getName()``, ``.addComputeMessage()``, and
        ``.getInputDataLocations()`` for identification, logging,
        and reading configured input data links.
    computeOptions : ComputeOptions object
        Run configuration carrying the runtime window, run directory,
        and simulation name.  Exposes ``.getRunTimeWindow()``,
        ``.getRunDirectory()``, and ``.getSimulationName()``.
    Returns
    -------
    result : bool or None
        True if both the target-temperature schedule write and the
        seven-outlet flow generation succeed.  None (implicit) if
        either step fails.
    """
    
    # Log the start of computation for this scripting alternative.
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')

    # Read the runtime window and extract the start/end times and forecast year.
    rtw = computeOptions.getRunTimeWindow()
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()
    startyear_str = starttime_str[5:9]
    year = int(startyear_str)
    
    ### Determine the alternative name ###
    # This step is necessary to set the model up correctly based on American configurations
    # Get the current watershed run directory.
    study_root = computeOptions.getRunDirectory()  # returns watershed directory
    
    # Strip off the scripts, alternative-group, and run folder to get the study root.
    study_main = study_root.split(os.sep)
    study_main = os.sep.join(study_main[:-3])
    
    # Build the path to the WAT simulation group directory.
    study_wat = study_main + os.sep + 'wat' + os.sep + 'simGroups'
    
    # List all simulation group files.
    all_groups = os.listdir(study_wat)
    
    # Keep only forecast simulation groups.
    forecast_groups = [x for x in all_groups if ".fsimgrp" in x]
    
    # Exclude backup files that still match the simulation-group extension.
    forecast_groups = [x for x in forecast_groups if '.bak' not in x]
    
    # Remove the file extension so only group names remain.
    forecast_groups = [x.split(".fsimgrp")[0] for x in forecast_groups]
    
    # Sort descending so later-named groups are checked first.
    forecast_groups.sort(reverse=True)
    
    # Get the simulation name, expected in alternative-group format.
    sim_name = computeOptions.getSimulationName()
    
    # Find which forecast group name appears inside the simulation name.
    for group_name in forecast_groups:
        
        # Once a matching group is found, stop searching.
        if group_name in sim_name:
            
            # Valid group name component has been found. Break to continue processing 
            break
    
    # Remove the group suffix from the simulation name to isolate the alternative name.
    alt_name = sim_name.split(group_name)[0]
    
    # Remove the trailing dash from the alternative name.
    alt_name = alt_name[:-1]
    
    # Construct a disk-friendly simulation name using underscores.
    alt_name_underscore = alt_name.replace(' ', '_')
    sim_name_underscore = alt_name_underscore + '-' + group_name

    ### Adjust the model name based on the alternative name ###
    # Get the current run directory.
    run_dir = computeOptions.getRunDirectory()

    # Start with the base model name. This is an assumed convention. 
    model_name = "W2 Folsom"  # base Folsom forecast model (iterative w/ bypass)
    
    # Append the FixedATSP suffix when that mode is encoded in the alternative name.
    if "FixedATSP" in alt_name:
        
        # FixedATSP is found in the alternative name. Append it to the run directory
        model_name += " FixedATSP"
    
    # Check the alternative if the model needs adjustment for the no bypass mode
    if "NoBypass" in alt_name:
        
        # NoBypass is found in the alternative name. Append it to the run directory
        model_name += " NoBypass"
    
    ### Setup the paths ###
    # Build model and shared-data file paths used by the forecast workflow.
    model_dir = os.path.join(study_main, 'cequal-w2', 'Folsom', model_name)
    targt_temp_npt_filepath = os.path.join(model_dir,'TargetSchedulesA.npt') # overwrite what's there    
    shared_dir = os.path.join(study_main,'shared')
    forecast_dss = os.path.join(shared_dir,'WTMP_American_forecast.dss')
    schedule_csv = os.path.join(model_dir,'SchedulesA.csv')

    # Normalize forecast meteorological DSS data in case units/types are inconsistent.
    DMS_preprocess.fix_DMS_types_units(forecast_dss)

    # Convert monthly Folsom storage to elevation using the storage-area lookup file.
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_Folsom.csv'), 'Folsom')
    fpp.storage_to_elev('Folsom',elev_stor_area,forecast_dss,'//FOLSOM/STORAGE//1Month/AMER_BC_SCRIPT/',conic=False)

    # invent Natoma elevation record from folsom storage rec
    fpp.invent_elevation('Natoma',forecast_dss,'//FOLSOM/STORAGE//1Month/AMER_BC_SCRIPT/',123.0)

    # Predict daily elevations so the model can start on an arbitrary day.
    fpp.write_forecast_elevations(currentAlternative, rtw, forecast_dss, shared_dir)

    # Adjust release accounting so flow through the dam excludes municipal/pump flow.
    fpp.subtract_muni_pump(forecast_dss)
    
    # Split evaporation terms into the format expected by the W2 model.
    fpp.split_folsom_evap(forecast_dss,'/AMERICAN RIVER/FOLSOM LAKE/FLOW-ACC-DEP//1Day/AMER_BC_SCRIPT/')

    # Split Nimbus dam flow into three gate-specific series for the Natoma W2 model.
    fpp.split_nimbus_outflow(forecast_dss,'/AMERICAN RIVER/LAKE NATOMA/FLOW-NIMBUS ACTUAL//1Day/AMER_BC_SCRIPT/')

    # Get the downstream target temperature location from DSS.
    TT_loc = fpp.get_downstream_loc(forecast_dss)

    # Build time values and load target-temperature driver data.
    # load Tair and Folsom FaveFlow
    start_date = HecTime(starttime_str) 
    start_doy = DSS_Tools.hectime_to_julian(start_date) 
    
    endtime_date = HecTime(endtime_str)
    endtime_doy = DSS_Tools.hectime_to_julian(endtime_date)
    
    doys,FaveFlow,Tair = fpp.load_tt_data(forecast_dss, starttime_str, endtime_str) # day-of-year,CMS,C
    
    # Write the target temperature schedule file used by the W2 model.
    target_temp_write = fpp.write_target_temp_npt(year,TT_loc,doys,Tair,FaveFlow,schedule_csv,targt_temp_npt_filepath,lagWatt=True)

    # Generate flow time series for the seven Folsom outlets.
    folsom_outlets = fpp.write_qot_7outlets_flows(forecast_dss, starttime_str, endtime_str)

    # Create small constant DSS records used for linking or default values.
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='TinyFlow',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='TinyFlow',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=10.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='TENS',fpart='TENS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    
    # Create constant gate elevation thresholds used by the outlet/gate logic.
    print('Making ELEV-INITIAL-GATE records...')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=401.0, what='elev', 
                        dss_type='PER-AVER', period='1DAY',cpart='ELEV-INITIAL-GATE',fpart='FULL-HEIGHT')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=362.0, what='elev', 
                        dss_type='PER-AVER', period='1DAY',cpart='ELEV-INITIAL-GATE',fpart='1-OUT')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=336.0, what='elev', 
                        dss_type='PER-AVER', period='1DAY',cpart='ELEV-INITIAL-GATE',fpart='2-OUT')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=307.0, what='elev', 
                        dss_type='PER-AVER', period='1DAY',cpart='ELEV-INITIAL-GATE',fpart='ALL-OUT')

    # Find configured shutter input locations and derive initial shutter elevations.
    # find initial shutter elevations
    locations_obj = currentAlternative.getInputDataLocations()
    shutter_objs = DSS_Tools.organizeLocations(currentAlternative, locations_obj, shutter_links, return_dss_paths=False)
    w2_elevs = fpp.get_initial_shutter_positions(rtw, currentAlternative, computeOptions, shutter_objs)
    print('Initial W2 shutter elevs:',w2_elevs)

    # Update W2 restart/control files with the restart date and initial shutter elevations.
    # update W2 interative control files for restart date - wait this functionality is already in the w2 plugin - did you check to see if folsom iterative compute option is checked??
    # also update the inital gate elevations in the w2_con file
    fpp.update_W2_Folsom_iterative_restart_date_and_shutters(rtw,model_dir,w2_elevs)

    # Update the schedule guess/ensemble selection used by the W2 run.
    # if this is a non-iterative simulation, put the correct ensemble number into the 'Set Guess' file 
    #if 'FixedATSP' in model_name:
    # -- edit, 2025-01-27, B. Saenz.  We want to set the  correct ensemble for all runs, as it also
    # tells the iterative version of Folsom W2 to start with a parituclar schedule, to save time if you have a good guess    w2run_base,_ = os.path.split(run_dir)
    w2run_base,_ = os.path.split(run_dir)
    fpp.update_W2_Folsom_iterative_schedule_number(w2run_base,model_dir)

    # Adjust the selective withdrawal file so auto shutter logic works from the chosen start date.
    # The W2 selective file must be modified depending on the start gate to get the auto shutter to work correctly
    fpp.modify_w2_selective_start_date(rtw,os.path.join(model_dir,'w2_selective.npt'))

    # Return success only if both target temperature writing and outlet generation succeeded.
    if target_temp_write and folsom_outlets:
        return True

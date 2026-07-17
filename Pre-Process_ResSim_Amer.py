import os, sys
from com.rma.model import Project
from hec.heclib.util import HecTime

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



def computeAlternative(currentAlternative, computeOptions):
    """
    Entry point for the WAT scripting alternative compute that pre-processes
    forecast DSS inputs for the American River HEC-ResSim Folsom model.

    This script prepares all input DSS records needed before ResSim executes,
    including elevation forecasts, flow disaggregation, equilibrium temperature,
    resampled hourly inputs, constant placeholder records, and gate configuration
    records. Optionally disables the lower river bypass outlet when the NoBypass
    model variant is selected.

    Workflow:
      1. Determines the forecast year and model variant (base ResSim or NoBypass)
         from the run directory name.
      2. Normalizes data types and units in the forecast DSS file via
         DMS_preprocess.fix_DMS_types_units.
      3. Computes and writes forecast reservoir elevations for Folsom and Natoma
         using fpp.write_forecast_elevations (mass-balance from starting elevation
         and inflow/outflow records).
      4. Subtracts municipal pumping from total Folsom release to produce the
         net downstream release via fpp.subtract_muni_pump.
      5. Computes Folsom outlet flow components (penstock, leakage, spill) for
         the 7-outlet model via fpp.write_qot_7outlets_flows.
      6. Computes the hourly equilibrium water surface temperature at Fair Oaks
         from air temperature, cloud cover, wind speed, solar irradiance, and
         dewpoint temperature via fpp.eq_temp.
      7. Resamples daily forecast inputs (pumping, NF/SF inflows, evaporation)
         to hourly resolution using DSS_Tools.resample_dss_ts.
      8. Writes constant DSS placeholder records for flow (2 CFS, 0 CFS),
         water temperature (10 deg C), evaporation (0), and gate elevation
         thresholds used by the selective withdrawal outlet logic.
      9. When the NoBypass variant is active, disables the lower river outlet
         usage record by overwriting it with -1.0 via
         fpp.remove_folsom_lower_river_use.

    Inputs:
      currentAlternative -- WAT scripting alternative object providing compute
                            messages and context
      computeOptions     -- WAT compute options object providing the run time
                            window, run directory, and simulation name

    Output:
      Returns True on successful completion.
      Writes the following categories of DSS records to the shared forecast DSS
      file (WTMP_American_forecast.dss):
        - Forecast reservoir elevations (Folsom monthly, daily, and predicted;
          Natoma constant 123.0 ft)
        - Daily Folsom release resampled to hourly
        - Net Folsom release to Natoma (total minus pumping)
        - Folsom outlet flow components (penstock, leakage, spill; CMS)
        - Hourly equilibrium temperature at Fair Oaks (deg C, INST-VAL);
          also resampled to daily and weekly
        - Hourly resampled pumping, NF/SF inflows, and evaporation records
        - Constant placeholder records: 2 CFS (tiny flow), 0 CFS (zeros),
          10 deg C (TENS), 0 evap (ZEROS)
        - Gate elevation constant records: 401 ft (FULL-HEIGHT), 362 ft (1-OUT),
          336 ft (2-OUT), 307 ft (ALL-OUT)
        - Lower river outlet usage set to -1.0 (NoBypass variant only)
    """
    
    # Log the start of the computation for this alternative.
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')

    # Read the runtime window and derive the forecast year from the start time.
    rtw = computeOptions.getRunTimeWindow()
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()
    startyear_str = starttime_str[5:9]
    year = int(startyear_str)

    # Determine the run directory and select the model variant.
    run_dir = computeOptions.getRunDirectory()
    model_name = 'ResSim'  # base Folsom forecast model
    if "nobypass" in run_dir.lower():
        model_name += " NoBypass"  
    
    # Build shared paths and locate the forecast DSS file.
    shared_dir = os.path.join(fpp.study_dir_from_run_dir(run_dir),'shared')
    forecast_dss = os.path.join(shared_dir,'WTMP_American_forecast.dss')

    # Normalize forecast met data types/units in case upstream data has inconsistencies.
    DMS_preprocess.fix_DMS_types_units(forecast_dss)

    # Generate forecast elevations so the model can start from an arbitrary day.
    fpp.write_forecast_elevations(currentAlternative, rtw, forecast_dss, shared_dir)

    # Adjust release accounting so dam flow excludes municipal/pump withdrawals.
    fpp.subtract_muni_pump(forecast_dss)

    # get target temp location from dummy dss file somehow?
    # Load basic time context for the forecast period.
    start_doy = HecTime(starttime_str).dayOfYear() #.value()
    endtime_doy = HecTime(endtime_str).dayOfYear() #.value()

    # Compute outlet flow series used by the downstream temperature workflow.
    # I don't think we need this split out; everything is dynamic including river outlets?
    folsom_outlets = fpp.write_qot_7outlets_flows(forecast_dss, starttime_str, endtime_str)

    # Create equilibrium temperature input needed for ResSim target temperature calculations.
    currentAlternative.addComputeMessage("Computing equilibrium temperature, this may take a while...")
    
    # eq_temp(rtw,at,cl,ws,sr,td,eq_temp_out)
    fpp.eq_temp(rtw,
            [forecast_dss,"/MR AM.-NATOMA LAKE/FAIR OAKS/TEMP-AIR//1HOUR/251.40.53.1.1/"],
            [forecast_dss,"/MR AM.-LOWER AMERICAN R./MATHER AFB/%-CLOUD COVER-FRAC//1Hour/252.50.129.1.1/"],
            [forecast_dss,"/MR AM.-NATOMA LAKE/FAIR OAKS/SPEED-WIND//1HOUR/251.40.133.1.1/"],
            [forecast_dss,"/MR AM.-FOLSOM LAKE/THEORETICAL FOLSOM SOLAR/IRRAD-SOLAR//1HOUR/250.15.135.1.1/"],
            [forecast_dss,"/MR AM.-NATOMA LAKE/FAIR OAKS/TEMP-DEWPOINT//1HOUR/251.40.51.1.1/"],
            [forecast_dss,"/MR AM.-NATOMA LAKE/FAIR OAKS/Temp-Equil//1Hour/amer_bc_script/"]
           )

    # Resample daily forecast inputs to hourly resolution for elevation and mass-balance calculations.
    # write an hourly forecast elevation based on starting elevation and flows
    DSS_Tools.resample_dss_ts(forecast_dss,'//FOLSOM/PUMPING (FP)//1Day/AMER_BC_SCRIPT/',rtw,forecast_dss,'1HOUR')
    DSS_Tools.resample_dss_ts(forecast_dss,'//Folsom-NF-in/FLOW-IN//1Day/AMER_BC_SCRIPT/',rtw,forecast_dss,'1HOUR')
    DSS_Tools.resample_dss_ts(forecast_dss,'//Folsom-SF-in/FLOW-IN//1Day/AMER_BC_SCRIPT/',rtw,forecast_dss,'1HOUR')
    DSS_Tools.resample_dss_ts(forecast_dss,'/AMERICAN RIVER/FOLSOM LAKE/FLOW-ACC-DEP//1Day/AMER_BC_SCRIPT/',rtw,forecast_dss,'1HOUR')
    
    # Define inflow and outflow records used in subsequent water balance logic.
    inflow_records = ['//Folsom-NF-in/FLOW-IN//1Hour/AMER_BC_SCRIPT/',
                      '//Folsom-SF-in/FLOW-IN//1Hour/AMER_BC_SCRIPT/',
                      '/AMERICAN RIVER/FOLSOM LAKE/FLOW-ACC-DEP//1Hour/AMER_BC_SCRIPT/']  # this actually evap, but negative already, so it goes as inflow
    outflow_records = ['//FOLSOM/FLOW-RELEASE//1Hour/AMER_BC_SCRIPT/']
    
    # Create constant DSS records used as placeholders, defaults, or linking aids.
    # make constant records for linking
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=2.0, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='2cfs',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=2.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='2cfs',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=10.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='TENS',fpart='TENS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='evap', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    
    # Create gate elevation threshold records used by outlet logic.
    print('Making ELEV-INITIAL-GATE records...')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=401.0, what='elev', 
                        dss_type='PER-AVER', period='1DAY',cpart='ELEV-INITIAL-GATE',fpart='FULL-HEIGHT')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=362.0, what='elev', 
                        dss_type='PER-AVER', period='1DAY',cpart='ELEV-INITIAL-GATE',fpart='1-OUT')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=336.0, what='elev', 
                        dss_type='PER-AVER', period='1DAY',cpart='ELEV-INITIAL-GATE',fpart='2-OUT')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=307.0, what='elev', 
                        dss_type='PER-AVER', period='1DAY',cpart='ELEV-INITIAL-GATE',fpart='ALL-OUT')
    
    # If the selected model is the no-bypass variant, remove lower river bypass usage.
    if 'nobypass' in model_name.lower():
        print('Setting: Do not use Folsom River Bypass in forecast')
        fpp.remove_folsom_lower_river_use(forecast_dss,"/FOLSOM/LOWER RIVER OUTLET USEAGE/COUNT//1Day/AMER_BC_SCRIPT/")

    # Signal successful completion.
    return True

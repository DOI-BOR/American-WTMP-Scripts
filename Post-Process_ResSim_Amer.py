import sys
print(sys.path)

from hec.heclib.dss import HecDss

import flowweightaverage
reload(flowweightaverage)

import DSS_Tools
reload(DSS_Tools)

# computeAlternative function is called when the ScriptingAlternative is computed.
#
# Arguments:
#   currentAlternative - the ScriptingAlternative.
#                        hec2.wat.plugin.java.impl.scripting.model.ScriptPluginAlt
#
#   computeOptions     - compute configuration supplied by HEC-WAT.
#                        hec.wat.model.ComputeOptions
#
# Return:
#   True  = successful execution
#   False = execution failure
#
# Note:
#   If no explicit return value is provided, HEC-WAT treats the
#   script as having completed successfully.
#
# Purpose
# -------
# This script post-processes ResSim Folsom release temperatures and
# generates flow-weighted outlet temperature products for downstream
# forecasting and visualization.
#
# Major Outputs
# -------------
# 1. Flow-weighted Penstock 1-3 release temperature.
# 2. Daily flow-weighted Penstock 1-3 release temperature.
# 3. Daily full-dam release temperature.
# 4. Target temperature records under ResSim F-part.
# 5. Penstock flow-minus-leakage records for buzz plots.
#
##

# Input location configuration.
#
# Each entry contains:
#   [Flow Location, Temperature Location]
#

# These paired records are used by flow-weighting routines.
ResSimFolsomInputs = [
  ['P1 Flow','P1 Temp'],
  ['P2 Flow','P2 Temp'],
  ['P3 Flow','P3 Temp'],
  ['Folsom Flow','Folsom Temp']
]


def computeAlternative(currentAlternative, computeOptions):
    """
    Entry point for the WAT scripting alternative compute for the HEC-ResSim Folsom
    release temperature post-processing workflow.

    Orchestrates five sequential operations using ResSim Folsom penstock flow and
    temperature input DSS locations:

      1. Flow-weighted Penstock 1-3 release temperature (hourly, FWA)
         Computes a combined hourly release temperature for penstocks 1-3 using
         flow-weighted averaging. Flows below 50 CFS are excluded.

      2. Daily flow-weighted Penstock 1-3 release temperature (FWA2_Daily)
         Creates a daily version of the penstock release temperature for use by
         downstream temperature forecasting workflows. A 1-day + 1-hour delay
         is applied to skip garbage values at the start of the record.

      3. Daily full-dam release temperature (FWA2_Daily)
         Computes a daily flow-weighted average temperature using the combined
         Folsom release record rather than individual penstocks. A 1-day + 1-hour
         delay is applied because the source record contains unusable values until
         approximately 01:00 on Day 2.

      4. Target temperature record under ResSim F-part
         Reads the standard American River forecast target temperature record
         (AMER_BC_SCRIPT F-part), converts Fahrenheit to Celsius if needed, and
         writes a duplicate copy under the ResSim F-part for compatibility with
         plotting tools that locate records by F-part.

      5. Penstock flow-minus-leakage diagnostic records (buzz plots)
         For each of penstocks 1-3, subtracts the leakage flow (Penstock 1 only;
         zero record used for penstocks 2 and 3) from the penstock flow to produce
         adjusted flow records used by operational buzz plots. A 1-day + 1-hour
         delay is applied for timing alignment.

    All output DSS records are written under the ResSim F-part (extracted from the
    first input location) so they appear alongside native ResSim outputs in DSS
    plotting and reporting tools.

    Inputs:
      currentAlternative -- WAT scripting alternative object providing input/output
                            data locations, compute messages, and DSS path creation
      computeOptions     -- WAT compute options object providing the DSS filename,
                            run time window, and run directory

    Output:
      Returns True on successful completion.
      Writes the following DSS records to the forecast DSS file (all under the
      ResSim F-part):
        - Folsom Penstock 1-3 flow-weighted temperature      (hourly, CFS-weighted)
        - Folsom Penstock 1-3 flow-weighted temperature      (1-day average)
        - Folsom full-dam flow-weighted temperature          (1-day average)
        - Folsom Downstream TEMP-WATER-TARGET                (1-hour, duplicate under ResSim F-part)
        - Folsom Lake-Penstock 1/2/3 minus leakage flow      (1-hour, per-penstock diagnostic)
    """
    
    
    #######################################################################
    # Script initialization
    #
    # Record script execution in the HEC-WAT compute log.
    #######################################################################
    
    
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())

    # Retrieve configured input DSS locations.
    locations_obj = currentAlternative.getInputDataLocations()

    # Organize locations into flow/temperature pairs and return DSS paths.
    #
    # Result structure:
    #
    # [
    #   [P1 Flow Path, P1 Temp Path],
    #   [P2 Flow Path, P2 Temp Path],
    #   [P3 Flow Path, P3 Temp Path],
    #   [Full Dam Flow Path, Full Dam Temp Path]
    # ]
    locations_paths = DSS_Tools.organizeLocationsPaired(currentAlternative, locations_obj, ResSimFolsomInputs, return_dss_paths=True)
    
    # Diagnostic output to console.
    print('locations_paths:')
    print(locations_paths)

    
    #######################################################################
    # Determine ResSim F-part
    #
    # Many output records are intentionally written under the same F-part
    # as the ResSim source records so they can be displayed alongside
    # native ResSim outputs.
    #######################################################################
    
    
    # Use the first input record as the source of the ResSim F-part.
    ressim_fpart = locations_paths[0][0].split('/')[6]
    
    # Echo all discovered DSS paths to the compute log.
    currentAlternative.addComputeMessage('Found DSS paths:')
    
    for location in locations_paths:
        
        # Each location contains a flow and temperature DSS record.
        for path in location:
            currentAlternative.addComputeMessage(str(path))
    
    # Blank line for readability in compute messages.
    currentAlternative.addComputeMessage('\n')
    
    #######################################################################
    # Runtime configuration
    #######################################################################
    
    # Forecast DSS output file.
    dss_file = computeOptions.getDssFilename()
    
    # Forecast runtime window.
    rtw = computeOptions.getRunTimeWindow()

    #######################################################################
    # Build output DSS path
    #
    # The first output location is assumed to represent the
    # flow-weighted Penstock 1-3 release temperature.
    #######################################################################
    
    outputlocations = currentAlternative.getOutputDataLocations()
    
    # Create DSS pathname from configured output location.
    outputpath = currentAlternative.createOutputTimeSeries(outputlocations[0])

    # Warn if more output locations exist than expected.
    if len(outputlocations) > 1:
        
        currentAlternative.addComputeMessage(
            "Found more than 1 output datapath locations. Using the first, {0}".format(outputlocations[0]))
    
    # Convert pathname into editable DSS components.
    tspath = str(outputpath)
    print('tspath: '+tspath)
    
    tspath = tspath.split('/')
    
    # Existing F-part.
    fpart = tspath[6]
   
    # Replace scripting-alternative F-part with ResSim F-part.
    #
    # This ensures resulting records appear with ResSim outputs.
    tspath[6] = ressim_fpart
    
    outputpath = '/'.join(tspath)
    print('outpath: '+outputpath)
    
    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpath))

    #######################################################################
    # Flow-weighted Penstock 1-3 temperature
    #
    # Generates an hourly temperature representing the combined release
    # from the three powerhouse penstocks.
    #######################################################################
    
    # Ignore very small flows when computing weighted averages.
    cfs_limit = 50.0  # float
    
    flowweightaverage.FWA(currentAlternative, dss_file, rtw, locations_paths[0:3], outputpath, cfs_limit)

    #######################################################################
    # Daily flow-weighted Penstock 1-3 temperature
    #
    # Creates a daily version of the release temperature used by
    # downstream temperature forecasting workflows.
    #######################################################################
    
    tspath = outputpath.split('/')
    
    # Force DSS interval to daily.
    tspath[5] = '1DAY'
    
    dailyoutputpath = '/'.join(tspath)
    
    flowweightaverage.FWA2_Daily(currentAlternative, dss_file, rtw, locations_paths[0:3], dailyoutputpath, 
                                cfs_limit, delay_days=1,delay_hours=1)

    #######################################################################
    # Full-dam outflow temperature
    #
    # Uses the combined Folsom release record instead of individual
    # penstocks.
    #######################################################################
    
    # flow-weight ave full dam outflow temp
    outputpath_full_dam = currentAlternative.createOutputTimeSeries(outputlocations[1])    
   
    tspath = str(outputpath_full_dam)
    print('tspath: '+tspath)
   
    tspath = tspath.split('/')
    fpart = tspath[6]
    
    # Preserve ResSim F-part.
    tspath[6] = ressim_fpart
    
    # Daily interval output.
    tspath[5] = '1DAY'
    
    outputpath_full_dam = '/'.join(tspath)
    
    print('outputpath_full_dam: '+outputpath_full_dam)
    
    # Process the full-dam release record.
    #
    # Historical note:
    # Source record contains unusable values until approximately
    # 01:00 on Day 2, hence the delay parameters.
    
    flowweightaverage.FWA2_Daily(currentAlternative, dss_file, rtw, [locations_paths[3]], outputpath_full_dam,
                                0.0, delay_days=1,delay_hours=1) # there is garbage in this record until 1:00 on second day


    #######################################################################
    # Save Target Temperature Record Under ResSim F-Part
    #
    # The downstream target temperature record is normally produced
    # elsewhere in the forecasting workflow. This section creates a
    # duplicate copy using the ResSim F-part so it appears alongside
    # ResSim-derived outputs in DSS plotting tools.
    #
    # This is primarily a visualization and reporting convenience.
    #######################################################################
    
    # Open the forecast DSS file for reading and writing.
    dssOut = HecDss.open(dss_file)
    
    # Standard target temperature record used throughout the American
    # River forecasting workflow.
    #
    # Assumption:
    # This pathname remains stable across forecast configurations.
    target_t_sched = "/WTMP_American/Folsom Downstream/TEMP-WATER-TARGET//1Hour/AMER_BC_SCRIPT/"
    
    # Retrieve the target temperature time series.
    tsc_sched = dssOut.get(target_t_sched,True)
    
    #######################################################################
    # Unit normalization
    #
    # Forecast workflows generally operate in Celsius.
    #
    # Some DSS records may still be stored in Fahrenheit depending on
    # source model configuration.
    #######################################################################
    
    if tsc_sched.units.lower() == 'f':
        # Temporary Celsius storage.
        tc = []
        
        # Convert each value individually.
        for tf in tsc_sched.values:
            tc.append((tf-32.)*5.0/9.0)
        
        # Replace original values with Celsius values.
        tsc_sched.values = tc
        
        # Update metadata.
        tsc_sched.units = 'C'    
    
    #######################################################################
    # Create ResSim-version pathname
    #
    # The record contents remain identical. Only the F-part changes.
    #######################################################################

    # Break DSS pathname into components.
    out_parts = target_t_sched.split('/')
    
    # Replace original F-part with ResSim F-part.
    out_parts[6] = ressim_fpart
    
    # Assign new pathname.
    tsc_sched.fullName = '/'.join(out_parts)
    
    # Store duplicate record.
    dssOut.put(tsc_sched) # write out copy under ressim_fpart
    
    #######################################################################
    # Finished target temperature duplication
    #######################################################################
    
    dssOut.close()    

    #######################################################################
    # Create Penstock Flow Minus Leakage Records
    #
    # Purpose:
    # Generate diagnostic flow records used by operational buzz plots.
    #
    # Historical Assumption:
    # All modeled leakage is assigned to Penstock 1.
    #
    # Therefore:
    #
    #   Penstock 1 adjusted flow =
    #       Penstock 1 flow - leakage
    #
    #   Penstock 2 adjusted flow =
    #       Penstock 2 flow
    #
    #   Penstock 3 adjusted flow =
    #       Penstock 3 flow
    #
    # Penstocks 2 and 3 use a zero-valued subtraction record so that
    # all three penstocks can be processed using identical logic.
    #######################################################################

    # Loop through Penstocks 1 through 3.
    for j in range(1,4):
        
        ###################################################################
        # Output record name
        ###################################################################
        out_rec = '//Folsom Lake-Penstock %i minus leakage/Flow//1Hour/%s/'%(j,ressim_fpart)
        
        ###################################################################
        # Determine subtraction record
        #
        # Penstock 1 receives leakage correction.
        # Penstocks 2 and 3 subtract a zero record.
        ###################################################################
        if j==1:
            subtraction_rec = '//Folsom_leakage/FLOW-LEAKAGE//1Hour/%s/'%(ressim_fpart)
        else:
            
            # Zero-valued placeholder record.
            #
            # Historical note:
            # This record is reportedly created automatically by
            # ResSim from a daily zero input.
            #
            # Future maintainers should verify that this pathname
            # remains stable across ResSim versions.
            
            subtraction_rec = '//ZEROS/flow//1Hour/ZEROS-%s/'%(ressim_fpart)  # made by ResSim by upscaling a daily input.  will it always be in results with this name?
        
        ###################################################################
        # Source records
        #
        # Record[0] = Penstock flow
        # Record[1] = Leakage or zero record
        ###################################################################
        
        flow_Records = ['//Folsom Lake-Penstock %i/Flow//1Hour/%s/'%(j,ressim_fpart), subtraction_rec]
        
        ###################################################################
        # Generate adjusted flow record
        #
        # add_or_subtract_flows() is being used here as a generic
        # arithmetic utility.
        #
        # Result:
        #
        #     output = flow - subtraction
        #
        # Delay parameters align timing with other forecast products.
        ###################################################################
        
        DSS_Tools.add_or_subtract_flows(currentAlternative, rtw, flow_Records, dss_file, [None,False], 
                                        out_rec, dss_file, multiplier=[1.0,1.0],delay_days=1.0,delay_hours=1.0)

    
    #######################################################################
    # Successful completion
    #
    # Outputs generated:
    #
    #   1. Hourly flow-weighted Penstock 1-3 temperature
    #   2. Daily flow-weighted Penstock 1-3 temperature
    #   3. Daily full-dam release temperature
    #   4. Target temperature record under ResSim F-part
    #   5. Penstock flow-minus-leakage diagnostic records
    #
    # Returning True indicates successful completion to the
    # HEC-WAT scripting framework.
    #######################################################################
    
    return True





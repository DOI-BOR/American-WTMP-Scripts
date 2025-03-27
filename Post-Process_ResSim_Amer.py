import sys
print(sys.path)

from hec.heclib.dss import HecDss

import flowweightaverage
reload(flowweightaverage)

import DSS_Tools
reload(DSS_Tools)

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

ResSimFolsomInputs = [
  ['P1 Flow','P1 Temp'],
  ['P2 Flow','P2 Temp'],
  ['P3 Flow','P3 Temp'],
  ['Folsom Flow','Folsom Temp']
]


def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    #locations = currentAlternative.getInputDataLocations()
    #locations_paths = flowweightaverage.organizeLocations(currentAlternative, locations)

    locations_obj = currentAlternative.getInputDataLocations()
    locations_paths = DSS_Tools.organizeLocationsPaired(currentAlternative, locations_obj, ResSimFolsomInputs, return_dss_paths=True)
    
    print('locations_paths:')
    print(locations_paths)

    # get ResSim F-part from input locations
    ressim_fpart = locations_paths[0][0].split('/')[6]
    
    currentAlternative.addComputeMessage('Found DSS paths:')
    for location in locations_paths:
        for path in location:
            currentAlternative.addComputeMessage(str(path))

    currentAlternative.addComputeMessage('\n')
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()

    outputlocations = currentAlternative.getOutputDataLocations()
    outputpath = currentAlternative.createOutputTimeSeries(outputlocations[0])

    if len(outputlocations) > 1:
        currentAlternative.addComputeMessage(
            "Found more than 1 output datapath locations. Using the first, {0}".format(outputlocations[0]))
    tspath = str(outputpath)
    print('tspath: '+tspath)
    tspath = tspath.split('/')
    fpart = tspath[6]
    #    new_fpart = fpart.lower().split(':scripting-')[0] #fpart will have scripting in it, but we want w2 version
    #    new_fpart += ':ResSim-' #replace scripting with ressim
    #    new_fpart += '-'.join(computeOptions.getSimulationName().split('-')[:-1]) #sim name will have the sim group, so snip that
    tspath[6] = ressim_fpart
    outputpath = '/'.join(tspath)
    print('outpath: '+outputpath)
    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpath))

    # flow-weight ave three penstock temps
    cfs_limit = 50.0  # float
    flowweightaverage.FWA(currentAlternative, dss_file, rtw, locations_paths[0:3], outputpath, cfs_limit)

    tspath = outputpath.split('/')
    tspath[5] = '1DAY'
    dailyoutputpath = '/'.join(tspath)
    flowweightaverage.FWA2_Daily(currentAlternative, dss_file, rtw, locations_paths[0:3], dailyoutputpath, 
                                cfs_limit, delay_days=1,delay_hours=1)

    # flow-weight ave full dam outflow temp
    outputpath_full_dam = currentAlternative.createOutputTimeSeries(outputlocations[1])    
    tspath = str(outputpath_full_dam)
    print('tspath: '+tspath)
    tspath = tspath.split('/')
    fpart = tspath[6]
    #    new_fpart = fpart.lower().split(':scripting-')[0] #fpart will have scripting in it, but we want w2 version
    #    new_fpart += ':ResSim-' #replace scripting with ressim
    #    new_fpart += '-'.join(computeOptions.getSimulationName().split('-')[:-1]) #sim name will have the sim group, so snip that
    tspath[6] = ressim_fpart
    tspath[5] = '1DAY'
    outputpath_full_dam = '/'.join(tspath)
    print('outputpath_full_dam: '+outputpath_full_dam)
    flowweightaverage.FWA2_Daily(currentAlternative, dss_file, rtw, [locations_paths[3]], outputpath_full_dam,
                                0.0, delay_days=1,delay_hours=1) # there is garbage in this record until 1:00 on second day


    # Save Target Temp record with ResSim F-part for plotting
    # ----------------------------------------------------------------
    dssOut = HecDss.open(dss_file)
    target_t_sched = "/WTMP_American/Folsom Downstream/TEMP-WATER-TARGET//1Hour/AMER_BC_SCRIPT/"
    tsc_sched = dssOut.get(target_t_sched,True)
    
    # convert to C just so we don't go crazy using DSS
    if tsc_sched.units.lower() == 'f':
        tc = []
        for tf in tsc_sched.values:
            tc.append((tf-32.)*5.0/9.0)
        tsc_sched.values = tc
        tsc_sched.units = 'C'    
    out_parts = target_t_sched.split('/')
    out_parts[6] = ressim_fpart
    tsc_sched.fullName = '/'.join(out_parts)
    dssOut.put(tsc_sched) # write out copy under ressim_fpart
    
    dssOut.close()    

    # Create Penstock Flow records minus leakage for buzz plots
    # Currently, all leakage is added to penstock one
    for j in range(1,4):
        out_rec = '//Folsom Lake-Penstock %i minus leakage/Flow//1Hour/%s/'%(j,ressim_fpart)
        if j==1:
            subtraction_rec = '//Folsom_leakage/FLOW-LEAKAGE//1Hour/%s/'%(ressim_fpart)
        else:
            subtraction_rec = '//ZEROS/flow//1Hour/ZEROS-%s/'%(ressim_fpart)  # made by ResSim by upscaling a daily input.  will it always be in results with this name?
        flow_Records = ['//Folsom Lake-Penstock %i/Flow//1Hour/%s/'%(j,ressim_fpart), subtraction_rec]
        DSS_Tools.add_or_subtract_flows(currentAlternative, rtw, flow_Records, dss_file, [None,False], 
                                        out_rec, dss_file, multiplier=[1.0,1.0],delay_days=1.0,delay_hours=1.0)

    return True





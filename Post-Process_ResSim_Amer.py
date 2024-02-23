import sys
print(sys.path)

import flowweightaverage
reload(flowweightaverage)

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
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    locations = currentAlternative.getInputDataLocations()
    locations_paths = flowweightaverage.organizeLocations(currentAlternative, locations)

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
    flowweightaverage.FWA_Daily(currentAlternative, dss_file, rtw, locations_paths[0:3], dailyoutputpath, 
                                cfs_limit, delay_days=1.0)

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
    flowweightaverage.FWA_Daily(currentAlternative, dss_file, rtw, [locations_paths[3]], outputpath_full_dam,
                                1.0, delay_days=1.0)
    
    return True





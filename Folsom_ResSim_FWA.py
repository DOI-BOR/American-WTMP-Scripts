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
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
    locations = currentAlternative.getInputDataLocations()
    locations_paths = flowweightaverage.organizeLocations(currentAlternative, locations)
    currentAlternative.addComputeMessage('Found DSS paths:')
    for location in locations_paths:
        for path in location:
            currentAlternative.addComputeMessage(str(path))
    
    currentAlternative.addComputeMessage('\n')
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    outputFpart = 'PostProcessed'
    outputlocations = currentAlternative.getOutputDataLocations()
    outputpath = currentAlternative.createOutputTimeSeries(outputlocations[0])


    if len(outputlocations) > 1:
        currentAlternative.addComputeMessage("Found more than 1 output datapath locations. Using the first, {0}".format(outputlocations[0]))
    tspath =str(outputpath)
    tspath = tspath.split('/')
    fpart = tspath[6]
#    new_fpart = fpart.lower().split(':scripting-')[0] #fpart will have scripting in it, but we want w2 version
#    new_fpart += ':ResSim-' #replace scripting with ressim
#    new_fpart += '-'.join(computeOptions.getSimulationName().split('-')[:-1]) #sim name will have the sim group, so snip that
    tspath[6] = outputFpart
    outputpath = '/'.join(tspath)
    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpath))

    cfs_limit = 50.0 #float
    flowweightaverage.FWA(currentAlternative, dss_file, rtw, locations_paths, outputpath, cfs_limit)
    
    tspath = outputpath.split('/')
    tspath[5] = '1DAY'
    dailyoutputpath = '/'.join(tspath)
    flowweightaverage.FWA_Daily(currentAlternative, dss_file, rtw, locations_paths, dailyoutputpath, cfs_limit)
    return True


            


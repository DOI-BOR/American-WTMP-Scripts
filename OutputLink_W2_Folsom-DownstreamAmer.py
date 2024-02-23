import sys
print(sys.path)
from hec.heclib.dss import HecDss
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE
import math
import datetime as dt
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

def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
 
    locations = currentAlternative.getInputDataLocations()
    locations = flowweightaverage.organizeLocations(currentAlternative, locations)
    currentAlternative.addComputeMessage('Found DSS paths:')
    for location in locations:
        for path in location:
            currentAlternative.addComputeMessage(str(path))
    
    currentAlternative.addComputeMessage('\n')
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()

    outputlocations = currentAlternative.getOutputDataLocations()

    # Add Folsom outflows into Natoma
    # ------------------------------------------------------------------------------------
    flow_locations = [loc[0] for loc in locations]

    outputpath = currentAlternative.createOutputTimeSeries(outputlocations[0])
    tspath = str(outputpath)
    
    currentAlternative.addComputeMessage("Original outpath 0 {0}".format(outputpath))
    tspath = tspath.split('/')
    fpart = tspath[6]
#    new_fpart = fpart.lower().split(':scripting-')[0] #fpart will have scripting in it, but we want w2 version
#    new_fpart += ':cequalw2-' #replace scripting with W2
#    new_fpart += '-'.join(computeOptions.getSimulationName().split('-')[:-1]) #sim name will have the sim group, so snip that
    if '|' in fpart:
        # remove everything before | including |
        tspath[6] = fpart[fpart.find('|')+1:]
    outputpath = '/'.join(tspath)    
    
    DSS_Tools.add_flows(currentAlternative, rtw, flow_locations, dss_file, outputpath, dss_file)
    
    # Flow-weight average Folsom outflow temps for linking
    # ------------------------------------------------------------------------------------

    # first outpath is flow-weighted dam temperature
    outputpath = currentAlternative.createOutputTimeSeries(outputlocations[1])
    #currentAlternative.addComputeMessage("Original output location 0 {0}".format(outputlocations[0]))
    #currentAlternative.addComputeMessage("Original tsc fullname 0 {0}".format(outputpath.fullName))
    tspath =str(outputpath)
    
    currentAlternative.addComputeMessage("Original outpath 0 {0}".format(outputpath))
    tspath = tspath.split('/')
    fpart = tspath[6]
#    new_fpart = fpart.lower().split(':scripting-')[0] #fpart will have scripting in it, but we want w2 version
#    new_fpart += ':cequalw2-' #replace scripting with W2
#    new_fpart += '-'.join(computeOptions.getSimulationName().split('-')[:-1]) #sim name will have the sim group, so snip that
    if '|' in fpart:
        # remove everything before | including |
        tspath[6] = fpart[fpart.find('|')+1:]
    outputpath = '/'.join(tspath)
    #outputpath = str(outputpath)
    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpath))

    cfs_limit = 1.0 #float
    flowweightaverage.FWA2(currentAlternative, dss_file, rtw, locations, outputpath, cfs_limit, 10.0)

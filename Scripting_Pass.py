import random
from com.rma.io import DssFileManagerImpl

def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
    currentAlternative.addComputeMessage("--simulation scripting bypass--" )
    return True


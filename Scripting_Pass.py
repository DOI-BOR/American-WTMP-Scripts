import random
from com.rma.io import DssFileManagerImpl

def computeAlternative(currentAlternative, computeOptions):
    # Log the start of computation for this scripting alternative.
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
    
    # Write a message indicating this path is using the simulation scripting bypass.
    currentAlternative.addComputeMessage("--simulation scripting bypass--" )
    
    # Return success immediately since no additional processing is performed.
    return True


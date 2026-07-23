import random
from com.rma.io import DssFileManagerImpl

def computeAlternative(currentAlternative, computeOptions):
    """
    Bypass the full simulation pipeline and return immediate success.
    A lightweight scripting alternative that performs no data
    preprocessing or model setup.  It logs a start message, flags
    the run as a simulation scripting bypass, and returns True
    without executing any computation steps.
    Parameters
    ----------
    currentAlternative : Alternative object
        The scripting alternative being computed.  Exposes
        ``.getName()`` for identification and
        ``.addComputeMessage()`` for logging.
    computeOptions : ComputeOptions object
        Run configuration (time windows, simulation parameters,
        rule sets).  Accepted for interface compatibility but not
        used in this bypass implementation.
    Returns
    -------
    result : bool
        Always True, signaling unconditional success to the
        calling framework.
    """
    
    # Log the start of computation for this scripting alternative.
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
    
    # Write a message indicating this path is using the simulation scripting bypass.
    currentAlternative.addComputeMessage("--simulation scripting bypass--" )
    
    # Return success immediately since no additional processing is performed.
    return True



from hec.heclib.dss import HecDss
from hec.hecmath import HecMathException
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
from hec.heclib.util import HecTime
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
import hec.hecmath.TimeSeriesMath as tsmath
from rma.util.RMAConst import MISSING_DOUBLE
import math
import sys
import datetime as dt
import os, sys

from com.rma.io import DssFileManagerImpl
from com.rma.model import Project
#import hec.hecmath.TimeSeriesMath as tsmath

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

from com.rma.io import DssFileManagerImpl
from java.util import TimeZone

import Acc_Dep_ResSim_American
reload(Acc_Dep_ResSim_American)

import DMS_preprocess
reload(DMS_preprocess)


def computeAlternative(currentAlternative, computeOptions):
    """
    Orchestrate the full computation pipeline for a scripting alternative.
    Executes two sequential steps to data preprocessing and accumulated
    depletion computation to for the American River ResSim workflow.
    Logs the start of the run, delegates to external modules, and
    returns a success flag only if both steps complete without error.
    Parameters
    ----------
    currentAlternative : Alternative object
        The scripting alternative being computed. Exposes `.getName()`
        for identification and `.addComputeMessage()` for logging.
    computeOptions : ComputeOptions object
        Run configuration (time windows, simulation parameters, rule
        sets) passed through to each sub-step.
    Returns
    -------
    result : bool or None
        True if both preprocessing and accumulated depletion succeed.
        None (implicit) if either step fails, since no explicit
        ``return False`` is defined.
    """
    
    # Log the start of computation for this scripting alternative.
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')

    # Run preprocessing for the American ResSim workflow.
    data_preprocess = DMS_preprocess.preprocess_ResSim_American(currentAlternative, computeOptions)

    # Run accumulated depletion computation for the same workflow.
    acc_dep = Acc_Dep_ResSim_American.computeAlternative(currentAlternative, computeOptions)

    # Return success only if both preprocessing and accumulated depletion complete successfully.
    if data_preprocess and acc_dep:
        return True


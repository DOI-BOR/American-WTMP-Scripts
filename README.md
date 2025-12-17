# American WTMP Scripts
The HEC-WAT scripts folder for the American River within the Bureau of Reclamation (USBR) Water Temperature Modeling Platform (WTMP).
This repository is a dependency of the American WTMP Study repository.
These scripts are called at different points in the WTMP workflow.
These scripts are all in Jython, the implementation of Python in Java.

## W2 files
### Forecast
* Post-Process_W2_Fol-Nat.py
  * The post-processing script for the forecast W2
* Post-Process_W2_Folsom.py
  * The post-processing script for the forecast W2
* Pre-Process_W2_Amer.py
  * The pre-processing script for the forecast W2

### Hindcast
### Planning
### Shared or not specified
* OutputLink_W2_Folsom-DownstreamAmer.py

## ResSim files
### Forecast
* Post-Process_ResSim_Amer-RivOnly.py
  * The post-processing script for the forecast ResSim
* Pre-Process_ResSim_Amer.py
  * The pre-processing script for the forecast ResSim
### Hindcast
### Planning
### Shared or not specified
* Acc_Dep_ResSim_American.py
* Folsom_ResSim_FWA.py
* Post-Process_ResSim_Amer.py
  * The post-processing script for ResSim
* Pre-Process_ResSim_AmerPO.py
  * The pre-processing script for ResSim


## Shared files between W2 and ResSim
These files contain various function that are shared by different models.
* BoundaryFixes.py
* create_balance_flow_jython.py
* DMS_preprocess.py
* DSS_Tools.py
* equilibrium_temp.py
* flowweightaverage.py
* Scripting_Pass.py
* Simple_DSS_Functions.py
* tz_offset.py

### Forecast
* Forecast_preprocess.py

## Dependencies
There are dependencies that come from the WAT and from external locations that must be brought in. 

- hec
  - hec.heclib
  - hec.io
  - hec.hecmath
- rma
  - com.rma
  - rma.util.RMAConst

## Usage
### Usage withing the build process
To be added later.

### Post build implementation
After making any desired changes to the code, files must be replaced in the scripts folder in a WTMP American Study folder. 
Scripts will be triggered as various points of the modeling workflow.

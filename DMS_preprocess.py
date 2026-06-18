
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

import DSS_Tools
reload(DSS_Tools)

import Simple_DSS_Functions as sdf
reload(sdf)

units_need_fixing = ['tenths','deg','kph','fract'] #'radians',]

def fix_DMS_types_units(dss_file):
    '''This method was implemented to change data types to PER-AVER that are not coming from the DMS that way'''
    
    # Get sanitized list of DSS records (filtered/cleaned path list)
    recs = DSS_Tools.get_sanitized_record_list(dss_file)

    # Open DSS file for reading/writing
    dss = HecDss.open(dss_file)
    
    # Loop through all DSS records
    for r in recs:
        rlow = r.lower()
        
        # Skip records that are not relevant to flow/unit fixes
        if not '/location info' in rlow and not '/temp-equil' in rlow and \
          not '/depth-temp' in rlow and not 'icpathsmap' in rlow and \
          not '/downstream_control_loc' in rlow and not 'temp-water-target' in rlow:
            
            # Read time series container (keep metadata + data)
            tsc = dss.get(r,True)

            # Convert flow-type records or daily records to PER-AVER if needed
            if "/flow" in rlow or "/1day/" in rlow:
                if not "/storage" and not "/stor" in rlow:
                    tsc.type = 'PER-AVER'
                 
                    # Debug output showing conversion
                    print('FixDMS Write original name: ',rlow)
                    print('FixDMS Write output name: ',tsc.fullName)
                    
                    # Write updated record back to DSS
                    dss.write(tsc)
            
            # Normalize units string for comparison
            units = str(tsc.units).lower()  # just to make sure
            
            # Fix known problematic unit types
            if units in units_need_fixing:
                
                # Convert "tenths" to fractional representation
                if units == 'tenths':
                    # save off a copy of cloud record in 0-1 for ResSim
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-FRAC'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'FRAC'
                    
                    # Scale values from tenths to fraction
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 10.0
                    
                    # Debug output showing conversion
                    print('FixDMS Write original name: ',rlow)
                    print('FixDMS Write output name: ',tsc.fullName)
                    
                    # Write updated record back to DSS
                    dss.write(tsc)
                    
                # Convert radians to degrees copy
                if units == 'radians':
                    # save off a copy in deg
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-DEG'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'deg'
                    
                    # Convert rad → deg
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / (2*3.141592653589793) * 360.0
                    
                    # Debug output showing conversion                    
                    print('FixDMS Write original name: ',rlow)
                    print('FixDMS Write output name: ',tsc.fullName)
                    
                    # Write updated record back to DSS
                    dss.write(tsc)
                    
                # Convert degrees to radians copy
                if units == 'deg':
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-RADIANS'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'radians'
                    
                    # Convert deg → rad
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 360.0 * (2*3.141592653589793)
                    
                    # Debug output showing conversion                    
                    print('FixDMS Write original name: ',rlow)
                    print('FixDMS Write output name: ',tsc.fullName)
                    
                    # Write updated record back to DSS
                    dss.write(tsc)
                    
                # Convert kilometer per hour to meter per second copy
                if units == 'kph':
                    # convert to m/s 
                    tsc.units = 'm/s'
                    
                    # Convert kph → m/s
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 3.6
                    
                    # Debug output showing conversion
                    print('FixDMS Write original name: ',rlow)
                    print('FixDMS Write output name: ',tsc.fullName)
                    
                    # Write updated record back to DSS
                    dss.write(tsc)
                    
                # Convert fractional cloud values (fract → FRAC + tenths copy)
                if units == 'fract':
                    # save off a copy of cloud record in 0-1 for ResSim, with proper naming, reset orignial to tenths
                    original_fullName = tsc.fullName
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-FRAC'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'FRAC'
                    
                    # Debug output showing conversion
                    print('FixDMS Write original name: ',rlow)
                    print('FixDMS Write output name: ',tsc.fullName)
                    
                    # Write updated record back to DSS
                    dss.write(tsc)
                
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] * 10.0
                    
                    tsc.units = 'tenths'
                    tsc.fullName = original_fullName
                    
                    # Debug output showing conversion
                    print('FixDMS Write original name: ',rlow)
                    print('FixDMS Write output name: ',tsc.fullName)
                    
                    # Write updated record back to DSS
                    dss.write(tsc)
    
    # Close DSS file
    dss.close()


def fix_DMS_types_units_old(dss_file):
    '''This method was implemented to change data types to PER-AVER that are not coming from the DMS that way'''
    
    # Get sanitized record list
    recs = DSS_Tools.get_sanitized_record_list(dss_file)
    
    # Open DSS file
    dss = HecDss.open(dss_file)
    
    # Loop through records
    for r in recs:
        rlow = r.lower()
        
        # Skip irrelevant metadata/control records
        if not '/location info' in rlow and not '/temp-equil' in rlow and \
          not '/depth-temp' in rlow and not 'icpathsmap' in rlow and \
          not '/downstream_control_loc' in rlow:
        
            # Read full DSS time series
            tsm = dss.read(r)

            # Force flow and daily records to PER-AVER
            if "/flow" in rlow or "/1day/" in rlow:
                tsm.setType('PER-AVER')
                dss.write(tsm)
            
            # Fix units if they are known problematic types
            if tsm.getUnits().lower() in units_need_fixing:
                
                # Convert tenths → FRAC copy
                if tsm.getUnits() == 'tenths':
                    # save off a copy of cloud record in 0-1 for ResSim
                    tsc = tsm.getData()
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-FRAC'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'FRAC'
                    
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 10.0  
                        
                    dss.write(tsc)
                    
                # Convert radians → degrees copy
                if tsm.getUnits() == 'radians':
                    # save off a copy in deg
                    tsc = tsm.getData()
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-DEG'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'deg'
                    
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / (2*3.141592653589793) * 360.0                
                    
                    dss.write(tsc)
                
                # Convert fract handling (duplicate + scale logic)
                if tsm.getUnits().lower() == 'fract':
                    
                    tsc = tsm.getData()
                    original_fullName = tsc.fullName
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-FRAC'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'FRAC'
                    dss.write(tsc)
                
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] * 10.0
                    tsc.units = 'tenths'
                    tsc.fullName = original_fullName
                    
                    dss.write(tsc)

                # Convert degrees → radians copy
                if tsm.getUnits() == 'deg':
                    # save off a copy in redians
                    tsc = tsm.getData()
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-RADIANS'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'radians'
                    
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 360.0 * (2*3.141592653589793) 
                        
                    dss.write(tsc)
                    
                # Convert kph → m/s    
                if tsm.getUnits() == 'kph':
                    # convert to m/s 
                    tsc = tsm.getData()
                    tsc.units = 'm/s'
                    
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 3.6
                        
                    dss.write(tsc)

                # Create W2-linked m/s variant
                if tsm.getUnits() == 'm/s':
                    # make a copy divied by kph conversion as a hack to get W2 linking the wind speed correctly 
                    tsc = tsm.getData()
                    rec_parts = tsc.fullName.split('/')
                    
                    if not "w2link" in rec_parts[3].lower():
                        rec_parts[3] += '-W2link'
                        tsc.fullName = '/'.join(rec_parts)
                        
                        for i in range(len(tsc.values)) :
                            tsc.values[i] = tsc.values[i] / 3.6
                        
                        dss.write(tsc)
    
    # Close DSS file
    dss.close()


def standardize_bc_temp_water_to_C(dss_file,output_dss_file):
    '''Make copies of temp-water records in C (standardizing on C) for ResSim linking'''
   
    # Get DSS record list
    recs = DSS_Tools.get_sanitized_record_list(dss_file)
    
    # Open input DSS
    dss = HecDss.open(dss_file)

    # Decide whether to overwrite or write to new DSS
    if dss_file == output_dss_file:
        dss_out = dss
    else:
        dss_out = HecDss.open(output_dss_file)
    
    # Loop through all records
    for r in recs:
        rlow = r.lower()
        
        # Only process temperature-water records
        if '/temp-water' in rlow:

            # Read record
            tsc = dss.get(r,True)

            incoming_units = tsc.units.lower()
            
            # Create new Celsius record copy
            tsc = dss.get(r,True)
            rec_parts = tsc.fullName.split('/')
            rec_parts[3] += '-C'
            tsc.fullName = '/'.join(rec_parts)
            tsc.units = 'C'
            
            # Convert Fahrenheit → Celsius if needed            
            if incoming_units == 'f' or incoming_units == 'degf':                
                for i in range(len(tsc.values)) :
                    tsc.values[i] = (tsc.values[i] - 32.0)*5.0/9.0             
            
            # Write standardized record
            dss_out.put(tsc)
    
    # Close DSS handles
    dss.close()
    if dss_file != output_dss_file:
        dss_out.close()


def DMS_fix_units_types(hydro_dss,met_dss_file):
    
    # Run unit/type fixes on both hydro and meteorology DSS files
    fix_DMS_types_units(hydro_dss)
    fix_DMS_types_units(met_dss_file)



def interp(x, xp, fp, left=None, right=None):
    """
    One-dimensional linear interpolation.

    Returns the one-dimensional piecewise linear interpolant to a function
    with given values at discrete data-points.

    Parameters
    ----------
    x : array_like
        The x-coordinates at which to evaluate the interpolated values.
    xp : 1-D sequence of floats
        The x-coordinates of the data points, must be increasing.
    fp : 1-D sequence of floats
        The y-coordinates of the data points, same length as `xp`.
    left : float, optional
        Value to return for `x < xp[0]`, default is `fp[0]`.
    right : float, optional
        Value to return for `x > xp[-1]`, default is `fp[-1]`.

    Returns
    -------
    y : float or ndarray
        The interpolated values, same shape as `x`.
    """

    # Handle list input by recursive mapping
    if isinstance(x, list):
        return [interp(point, xp, fp, left, right) for point in x]
    
    # Set boundary defaults
    else:
        if left is None:
            left = fp[0]
        if right is None:
            right = fp[-1]
        
        # Left extrapolation
        if x < xp[0]:
            return left
            
        # Right extrapolation    
        elif x > xp[-1]:
            return right
            
        # Linear interpolation within range    
        else:
            for i in range(len(xp) - 1):
                if x >= xp[i] and x <= xp[i+1]:
                    
                    # Perform the linear interpolation
                    t = (x - xp[i]) / (xp[i+1] - xp[i])
                    return fp[i] + t * (fp[i+1] - fp[i])


def interp_monthly_coeff(coeffs):
    '''
    Convert a set of 12 monthly regression coefficients into an
    hourly coefficient series for an entire leap year (8784 hours).

    The regression equations were developed on a monthly basis, but
    hourly temperature calculations require smoothly varying coefficients.
    This function linearly interpolates between monthly midpoint values
    to generate an hourly coefficient record.

    Note:
    - Assumes leap year length (366 days = 8784 hours).
    - Padding is added at beginning/end so interpolation wraps smoothly
      between December and January.
    '''
    # Approximate midpoint day-of-year for each month.
    month_midpoints = [16, 45, 75, 105, 136, 166, 197, 228, 259, 289, 320, 350]    
    
    # Convert midpoint days into hours from start of year.
    month_midpoints_hours = [(month_midpoints[i]-1)*24 for i in range(12)]
    
    # Months with 31 days get shifted 12 hours to better represent
    # the true center of the month.
    for i in  [1,3,5,7,8,10,12]:
        month_midpoints_hours[i-1] += 12
        
    # Pad endpoints so interpolation wraps continuously
    # from December back to January.    
    month_midpoints_hours_padded = [-15.5*24.0,] + month_midpoints_hours + [8784.0 + 15.5*24.0,]
    
    # Duplicate December and January coefficients at the ends
    # to support the wraparound interpolation.
    coeffs_padded = [coeffs[-1],] + coeffs + [coeffs[0],]
    
    # Store interpolated hourly coefficient values.
    hourly_coeff = []
    
    # Generate hourly coefficient values for every hour in a leap year.
    for i in range(8784):
        hourly_coeff.append(interp(i,month_midpoints_hours_padded,coeffs_padded))
    return hourly_coeff

def american_NF_temp_array_hourly(hour, NF_cms, MF_cms, T_air):
    '''
    CARDNO/Stantec regression equation for estimating
    North Fork American River inflow temperature to Folsom Lake.

    Inputs:
        hour    - hour index within leap year
        NF_cms  - North Fork flow (cms)
        MF_cms  - Middle Fork flow (cms)
        T_air   - air temperature (deg C)

    Returns:
        Estimated water temperature (deg C)

    Regression form:

        Temp =
            C0 +
            C1 * log10(NF_flow) +
            C2 * log10(MF_flow) +
            C3 * AirTemp

    Monthly regression coefficients are interpolated
    to hourly values before application.
    '''
    
    
    NF_coeff = ([[ 3.774,  5.013,  7.568, 13.929, 19.23 , 22.008, 27.481, 26.076, 19.876, 11.463,  7.827,  3.52],
       [ 1.266,  2.088,  3.042,  1.493, -4.149, -2.19 ,  0.461, -0.056, -2.334,  0.665,  0.685,  -0.27],
       [-0.123, -2.308, -4.644, -5.956, -2.651, -4.32 , -8.106, -7.756, -4.285, -2.909, -1.342, 1.59],
       [ 0.209,  0.289,  0.336,  0.278,  0.279,  0.182,  0.071,  0.064, 0.107,  0.355,  0.367,  0.30]])
    
    # Interpolate monthly coefficients to hourly resolution.
    coeff0 = interp_monthly_coeff(NF_coeff[0])
    coeff1 = interp_monthly_coeff(NF_coeff[1])
    coeff2 = interp_monthly_coeff(NF_coeff[2])
    coeff3 = interp_monthly_coeff(NF_coeff[3])
   
    # Output water temperatures.
    nf_temp = []
    
    for i in range(len(hour)):
        h_idx = hour[i]
        
        # Diagnostic print useful during calibration/debugging.
        print(hour[i],NF_cms[i],MF_cms[i],T_air[i])
        
        nf_temp.append( coeff0[h_idx] + coeff1[h_idx] * math.log10(NF_cms[i]) + coeff2[h_idx] * math.log10(MF_cms[i]) + coeff3[h_idx] * T_air[i] )        
        
    return nf_temp

def american_SF_temp_array_hourly(hour, SF_cms, T_air):
    '''
    CARDNO/Stantec regression equation for estimating
    South Fork American River inflow temperature to Folsom Lake.

    Inputs:
        hour    - hour index within leap year
        SF_cms  - South Fork flow (cms)
        T_air   - air temperature (deg C)

    Returns:
        Estimated water temperature (deg C)

    Regression form:

        Temp =
            C0 +
            C1 * log10(SF_flow) +
            C2 * AirTemp

    Monthly coefficients are interpolated to hourly values.
    '''
    
    # These will later be interpolated from monthly to hourly resolution
    SF_coeff = ([[ 1.956,  3.894,  8.456, 12.605, 19.374, 22.03 , 23.604, 21.761, 17.663, 11.832,  6.521,  3.430],
       [ 1.374,  0.221, -1.422, -3.05 , -5.815, -6.605, -5.623, -5.196, -4.067, -2.665, -0.366,  0.755],
       [ 0.29 ,  0.282,  0.224,  0.223,  0.204,  0.216,  0.114,  0.105, 0.155,  0.299,  0.374,  0.358 ]])
    
    # Convert monthly regression coefficients into hourly coefficients.
    coeff0 = interp_monthly_coeff(SF_coeff[0])
    coeff1 = interp_monthly_coeff(SF_coeff[1])
    coeff2 = interp_monthly_coeff(SF_coeff[2])
    
    # Prepare output list for computed temperatures    
    sf_temp = []
    
    # For each hour index, compute temperature via regression:
    for i in range(len(hour)):
        h_idx = hour[i]
        
        sf_temp.append( coeff0[h_idx] + coeff1[h_idx] * math.log10(SF_cms[i]) + coeff2[h_idx] * T_air[i] )
    
    # Return sequence of hourly temperatures (deg C)
    return sf_temp

            
def american_NF_temp_array(month, NF_cms, MF_cms, T_air):
    '''
    Monthly version of the CARDNO/Stantec North Fork American
    River temperature regression.

    Unlike the hourly version, this function uses discrete monthly
    coefficient sets and does not interpolate between months.

    Inputs:
        month   - month number (1-12)
        NF_cms  - North Fork flow (cms)
        MF_cms  - Middle Fork flow (cms)
        T_air   - air temperature (deg C)

    Returns:
        Estimated water temperature (deg C)
    '''
    
    
    # Coefficient tuples per month:
    # [C0, C1*(log10(NF)), C2*(log10(MF)), C3*(AirTemp)]

    NF_coeff = [
        [3.774, 1.266, -0.123, 0.209],
        [5.013, 2.088, -2.308, 0.289],
        [7.568, 3.042, -4.644, 0.336],
        [13.929, 1.493, -5.956, 0.278],
        [19.23, -4.149, -2.651, 0.279],
        [22.008, -2.190, -4.320, 0.182],
        [27.481, 0.461, -8.106, 0.071],
        [26.076, -0.056, -7.756, 0.064],
        [19.876, -2.334, -4.285, 0.107],
        [11.463, 0.665, -2.909, 0.355],
        [7.827, 0.685, -1.342, 0.367],
        [3.52, -0.27, 1.59, 0.30]
    ]
    
    # Initialize output list
    nf_temp = []
    
    # For each entry, select monthly coefficients and apply regression
    for i in range(len(month)):
        coeff = NF_coeff[month[i]-1]
        
        nf_temp.append( coeff[0] + coeff[1] * math.log10(NF_cms[i]) + coeff[2] * math.log10(MF_cms[i]) + coeff[3] * T_air[i] )
    
    # Return monthly temperature estimates (deg C)    
    return nf_temp

def american_SF_temp_array(month, SF_cms, T_air):
    '''
    Monthly version of the CARDNO/Stantec South Fork American
    River temperature regression.

    Uses month-specific coefficients directly rather than
    interpolating coefficients to hourly resolution.

    Inputs:
        month   - month number (1-12)
        SF_cms  - South Fork flow (cms)
        T_air   - air temperature (deg C)

    Returns:
        Estimated water temperature (deg C)
    '''
    
    # Monthly coefficients [C0, C1*(log10(SF)), C2*(AirTemp)]
    SF_coeff = [
        [1.956, 1.374, 0.290],
        [3.894, 0.221, 0.282],
        [8.456, -1.422, 0.224],
        [12.605, -3.050, 0.223],
        [19.374, -5.815, 0.204],
        [22.03, -6.605, 0.216],
        [23.604, -5.623, 0.114],
        [21.761, -5.196, 0.105],
        [17.663, -4.067, 0.155],
        [11.832, -2.665, 0.299],
        [6.521, -0.366, 0.374],
        [3.430, 0.755, 0.358]
    ]
    
    # Initialize output list
    sf_temp = []
    
    # Use month-specific coefficients without interpolation
    for i in range(len(month)):
        coeff = SF_coeff[month[i]-1]
        
        sf_temp.append( coeff[0] + coeff[1] * math.log10(SF_cms[i]) + coeff[2] * T_air[i] )
    
    # Return monthly temperature estimates (deg C)    
    return sf_temp

def american_SC_temp_array(month):
    '''
    South Canal inflow temperature estimate.

    Unlike NF and SF inflows, no regression equation is used.
    Instead, monthly average temperatures developed by
    CARDNO/Stantec are applied.

    Source temperatures are stored in degrees Fahrenheit and
    converted to Celsius before being returned.

    Inputs:
        month - month number (1-12)

    Returns:
        Monthly average South Canal temperature (deg C)
    '''
    
    # Monthly average temperatures (deg F) to be converted to deg C
    SC_ave_temp = [
        46.02,
        46.48,
        48.94,
        49.83,
        52.32,
        55.61,
        59.43,
        63.05,
        64.82,
        60.24,
        53.48,
        48.53
    ]
    
    # Initialize output list
    sc_temp = []
    
    # Convert each month's Fahrenheit average to Celsius
    for i in range(len(month)):
        sc_temp.append( (SC_ave_temp[month[i]-1] - 32.0)*5.0/9.0 )
    
    # Return monthly average SC temperatures (deg C)    
    return sc_temp

def calc_folsom_inflow_temps(currentAlt, rtw, hydro_dss, met_dss_file, output_dss_file,hourly=True):
    '''
    Compute estimated Folsom inflow temperatures and write them
    to DSS records.

    Historical note:
        This routine is currently unused because measured DMS
        temperature records are available. It is retained because
        implementing the regression workflow in Jython required
        substantial effort and may be useful in the future.

    Workflow:
        1. Read NF, MF, SF flow records.
        2. Read Fair Oaks air temperature record.
        3. Convert units as needed.
        4. Build month/hour indices.
        5. Calculate NF, SF, and SC inflow temperatures.
        6. Write temperature time series to DSS.

    Output records:
        /Folsom/NF Inflow/Temp-Water/
        /Folsom/SF Inflow/Temp-Water/
        /Folsom/SC Inflow/Temp-Water/
    '''
    
    # Gather time window (start/end) in both string and HEC integer time formats
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()    
    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    currentAlt.addComputeMessage('Calculating Folsom Inflow Temperatures...')

    # hardcoded paths ... yuck
    # Define default record paths depending on hourly vs daily resolution
    if hourly:
        NF = '::'.join([output_dss_file,'/MR Am.-Folsom Lake/NF American River-Flow/Flow//1Hour/250.400.125.1.1/']) # should change to lake clementime for future use
        MF = '::'.join([output_dss_file,'/MR Am.-Folsom Lake/MF calc/Flow//1Hour/ResSim_PreProcess/']) # should change to foresthill for use
        SF = '::'.join([output_dss_file,'/MR Am.-Folsom Lake/SF American River-Flow/Flow//1Hour/250.402.125.1.1/']) # should change to placerville for use
        AT = '::'.join([met_dss_file,'/MR Am.-Natoma Lake/Fair Oaks-Air Temp/Temp-Air//1hour/251.40.53.1.1/'])
        output_period = '1Hour'
        
    else:
        NF = '::'.join([hydro_dss,'/MR Am.-Folsom Lake/NF American River-Flow/Flow//1Day/250.400.125.1.1/'])
        MF = '::'.join([output_dss_file,'/MR Am.-Folsom Lake/MF calc/Flow//1Day/ResSim_PreProcess/'])
        SF = '::'.join([hydro_dss,'/MR Am.-Folsom Lake/SF American River-Flow/Flow//1Day/250.402.125.1.1/'])
        AT = '::'.join([output_dss_file,'/MR Am.-Natoma Lake/Fair Oaks-Air Temp/Temp-Air//1DAY/251.40.53.1.1/'])
        output_period = '1Day'

    # hardcoded paths ... yuck
    # Override paths with shared American_inflows_6 DSS—ensures consistent sources
    shared_path,_ = os.path.split(hydro_dss)
    American_inflows_6 = os.path.join(shared_path,'American_inflows_6.dss')
    NF = '::'.join([American_inflows_6,'/11427000/NF AMERICAN/FLOW//1HOUR/USGS-CARDNO-FROM-DAILY/'])
    MF = '::'.join([American_inflows_6,'/11433300/MF AMERICAN/FLOW//1HOUR/USGS-CARDNO-FROM-DAILY/'])
    SF = '::'.join([American_inflows_6,'/11444500/SF AMERICAN/FLOW//1HOUR/USGS-CARDNO-FROM-DAILY/'])
    AT = '::'.join([met_dss_file,'/MR Am.-Natoma Lake/Fair Oaks-Air Temp/Temp-Air//1hour/251.40.53.1.1/'])
    output_period = '1Hour'

    # Read inputs from DSS; convert units as needed (air temp F→C, flow cfs→cms)
    print('Reading inflows')
    inputs = []
    
    # Loop over NF, MF, SF flow records (i=0,1,2) and air temperature record (i=3)
    for i,input_rec in enumerate([NF,MF,SF,AT]):
        
        # All records are expected in 'filepath::dssrecord' format; parse and read accordingly
        if '::' in input_rec:
            dss_file_alt,inflow_rec_alt = input_rec.split('::')
            dssFm_alt = HecDss.open(dss_file_alt)
            ts_data = dssFm_alt.read(inflow_rec_alt, starttime_str, endtime_str, False).getData()
            dssFm_alt.close()
            print(dss_file_alt)
        else:
            currentAlt.addComputeMessage('ERROR reading' + str(input_rec) + '[does it have form <filepath>::<dssrec>?]')
            sys.exit(-1)
            
        # Final item is air temperature; handle unit conversion
        if i==3:
            
            # Convert Fahrenheit to Celsius if needed; error out on unrecognised units
            if ts_data.units.lower() == 'f':
                currentAlt.addComputeMessage('Converting airtemp F to C')
                convvals = []
                for t in ts_data.values:
                    convvals.append((t - 32.0)*5.0/9.0)
                inputs.append(convvals)
            elif ts_data.units.lower() == 'c':
                inputs.append(ts_data.values)
            else:
                currentAlt.addComputeMessage('ERROR unknown units for air temperature ' + str(input_rec) + ' units ' + ts_data.units)
                sys.exit(-1)                   
        else: 

            # Flow series: normalize to cms
            # Pass through cms values directly; convert cfs to cms if needed; error on unknown units
            if ts_data.units.lower() == 'cms':
                inputs.append(ts_data.values)
            elif ts_data.units.lower() == 'cfs':
                currentAlt.addComputeMessage('Converting cfs to cms')
                convvals = []
                for flow in ts_data.values:
                    convvals.append(flow / 35.314666213)
                inputs.append(convvals)
            else:
                currentAlt.addComputeMessage('ERROR unknown units for flow ' + str(input_rec) + ' units ' + ts_data.units)
                sys.exit(-1)               
        
        # Retain the time array from the last record read (all series share the same time axis)
        hectimes = ts_data.times

    # Build month/day/hour indices from HEC times for regression inputs
    time_post = HecTime(HecTime.MINUTE_INCREMENT)
    month = []
    day = []
    hour = []
    
    # Decompose each HEC time stamp into calendar components used by the regression functions
    for time in hectimes:
        time_post.set(time)
        month.append(time_post.month())
        day.append(time_post.dayOfYear() + time_post.hour()/24.0)
        hour.append((time_post.dayOfYear()-1)*24 + time_post.hour())

    # Compute inflow temperatures for NF, SF, and SC; write per output period
    dssFm_out = HecDss.open(output_dss_file)
    
    for i in range(3):
        
        # Select output pathname and call the appropriate regression function for each tributary
        if i==0: # NF
            outrec = '/Folsom/NF Inflow/Temp-Water//%s/ResSim_PreProcess/'%output_period
            outdata = american_NF_temp_array_hourly(hour, inputs[0], inputs[1], inputs[3])
        elif i==1: # SF
            outrec = '/Folsom/SF Inflow/Temp-Water//%s/ResSim_PreProcess/'%output_period
            outdata = american_SF_temp_array_hourly(hour, inputs[2], inputs[3])
        else:  # SC temperatures (monthly averages converted to C)
            outrec = '/Folsom/SC Inflow/Temp-Water//%s/ResSim_PreProcess/'%output_period
            outdata = american_SC_temp_array(month)    
        
        # Package results into a TimeSeriesContainer and write to DSS
        tsc = TimeSeriesContainer()
        tsc.times = hectimes
        tsc.fullName = outrec
        tsc.values = outdata
        tsc.units = 'C'
        tsc.type = 'PER-AVER'
        tsc.numberValues = len(outdata)   

        # Finalize output file        
        dssFm_out.write(tsc)

    dssFm_out.close()


def compute_folsom_flows(currentAlternative, rtw, hydro_dss, output_dss_file):
    # Add Mormon Ravine and Newcastle PP flows (there is one inflow designed for these in ResSim)
    '''
    Generate derived flow records needed by downstream
    ResSim and CE-QUAL-W2 workflows.

    Major calculations:

    1. Combine Newcastle Powerplant and Mormon Ravine flows.

    2. Calculate North Arm inflow:
           Lake Clementine +
           Foresthill -
           American River Pump Station (ARPS)

    3. Sum upper-elevation Folsom outlet releases
       (G1-G4).

    4. Sum lower-elevation Folsom outlet releases
       (G5-G8).

    5. Sum Natoma powerhouse releases.

    Several records are written to the preprocessing DSS
    database for later model use.
    '''
    
    # Combine daily flows for the two upstream contributors
    inflow_records = ['/MR Am.-Folsom Lake/11425416 Newcastle PP-Daily Flow/Flow//1Day/250.114.125.1.1/',
                      '/MR Am.-Folsom Lake/11433930 Mormon Ravine-Daily Flow/Flow//1Day/250.115.125.1.1/']
    
    # Write combined record to preprocessing DSS
    DSS_Tools.add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Am.-Folsom Lake/MormonR_NewcastlePP_Sum/Flow//1Day/ResSim_PreProcess/', output_dss_file)
    
    # Resample monthly ARPS flows to daily; padding helps cover edges and gaps
    DSS_Tools.resample_dss_ts(hydro_dss,'/MR Am.-Folsom Lake/American River Pump Station (ARPS)-Flow/Flow//1Mon/250.401.125.1.1/',
                              rtw,output_dss_file,'1Day',pad_1mon=True)
                              
    out_rec = '/MR Am.-Folsom Lake/North Arm/Flow//1Day/ResSim_PreProcess/'
    
    shared_path,_ = os.path.split(hydro_dss)
    
    American_inflows_6 = os.path.join(shared_path,'American_inflows_6.dss')
    
    # Assemble input paths (two adds and one subtract) for North Arm computation
    flow_Records = ['/MR Am.-Folsom Lake/11427000 Lake Clementine Dam-Daily Flow/Flow//1Day/250.112.125.1.1/',
                    '/MR Am.-Folsom Lake/11433300 Foresthill-Daily Flow/Flow//1Day/250.113.125.1.1/',
                    #'/MR Am.-Folsom Lake-11433300 Foresthill-Daily Flow/Flow//1Day/250.113.125.1.1/',  # this one is not downloading yet
                    #'::'.join([American_inflows_6,'/11433300/MF AMERICAN/FLOW//1DAY/USGS-CARDNO/']),
                    '::'.join([output_dss_file,'/MR Am.-Folsom Lake/American River Pump Station (ARPS)-Flow/Flow//1Day/250.401.125.1.1/']),]
                    #'::'.join([American_inflows_6,'/FOLSOM/AMERICAN RIVER PUMP STATION/FLOW//1DAY/CDEC-CARDNO/']),] # this one is not downloading yet
    
    # Perform add/subtract combining with the specified operation flags    # [N.A.,ADD,SUBTRACT]
    DSS_Tools.add_or_subtract_flows(currentAlternative, rtw, flow_Records, hydro_dss, [None,True,False], out_rec, output_dss_file)
        
    # EID outflow - do we need to normalize to daily?
    # /MR Am.-Folsom Lake-EID Folsom Diversion-Diversion Flow/Flow/ --?

    # Sum upper-elevation river outlet releases (G1–G4), then enforce minimum
    out_rec = '/MR Am.-Folsom Lake/Upper_River_Outlets_Sum_min4/Flow//1Hour/ResSim_PreProcess/'
    outflow_records = ['/MR Am.-Folsom Lake/FOL-Outlet Release G1/Flow//1Hour/250.3.125.23.1/',
                       '/MR Am.-Folsom Lake/FOL-Outlet Release G2/Flow//1Hour/250.3.125.24.1/',
                       '/MR Am.-Folsom Lake/FOL-Outlet Release G3/Flow//1Hour/250.3.125.25.1/',
                       '/MR Am.-Folsom Lake/FOL-Outlet Release G4/Flow//1Hour/250.3.125.26.1/',]
    DSS_Tools.add_flows(currentAlternative, rtw, outflow_records, hydro_dss, out_rec, output_dss_file)
    DSS_Tools.min_ts(output_dss_file, out_rec, 4.0, output_dss_file, 'ResSim_PreProcess')

    # Sum lower-elevation river outlet releases (G5–G8), then enforce minimum
    out_rec = '/MR Am.-Folsom Lake/Lower_River_Outlets_Sum_min4/Flow//1Hour/ResSim_PreProcess/'
    outflow_records = ['/MR Am.-Folsom Lake/FOL-Outlet Release G5/Flow//1Hour/250.3.125.27.1/',
                       '/MR Am.-Folsom Lake/FOL-Outlet Release G6/Flow//1Hour/250.3.125.28.1/',
                       '/MR Am.-Folsom Lake/FOL-Outlet Release G7/Flow//1Hour/250.3.125.29.1/',
                       '/MR Am.-Folsom Lake/FOL-Outlet Release G8/Flow//1Hour/250.3.125.30.1/',]
    DSS_Tools.add_flows(currentAlternative, rtw, outflow_records, hydro_dss, out_rec, output_dss_file)
    DSS_Tools.min_ts(output_dss_file, out_rec, 4.0, output_dss_file, 'ResSim_PreProcess')

    # Sum Natoma powerhouse releases (units U1 + U2) to a combined series
    out_rec = '/MR Am.-Natoma Lake/NAT-Gen Release Sum/Flow//1Hour/ResSim_PreProcess/'
    outflow_records = ['/MR Am.-Natoma Lake/NAT-Gen Release U1/Flow//1Hour/251.4.125.3.1/',
                       '/MR Am.-Natoma Lake/NAT-Gen Release U2/Flow//1Hour/251.4.125.4.1/']
    DSS_Tools.add_flows(currentAlternative, rtw, outflow_records, hydro_dss, out_rec, output_dss_file)


def compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file):
    '''
    Placeholder for creation of plotting and diagnostic records.

    Currently no implementation exists.

    Intended future purpose:
        - Derived QA/QC records
        - Plot-friendly summary time series
        - Model diagnostic outputs
    '''
    
    pass


def preprocess_W2_American(currentAlternative, computeOptions):
    '''
    Preprocessing workflow used by the CE-QUAL-W2
    American River modeling framework.

    Responsibilities:
        - Establish project/run directories.
        - Locate hydro and meteorological DSS databases.
        - Standardize DSS units/types.
        - Generate derived Folsom flow records.
        - Generate plotting records.

    Output:
        DMS_American_Pre-Process.dss
    '''
    
    # Retrieve the run time window from compute options
    rtw = computeOptions.getRunTimeWindow()
    
    # Identify project directory and run directory for logging and file locations
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    
    # Get model timestep and shared directory path
    balance_period = currentAlternative.getTimeStep()
    shared_dir = os.path.join(project_dir, 'shared')

    # Define key DSS file paths used for preprocessing output
    output_dss_file = os.path.join(shared_dir,'DMS_American_Pre-Process.dss') 

    # Identify hydro and meteorological input DSS file
    hydro_dss = os.path.join(shared_dir, 'DMS_AmericanHydroTS.dss')
    fix_DMS_types_units(hydro_dss)
    met_dss_file = os.path.join(shared_dir,'DMS_AmericanMet.dss')
    fix_DMS_types_units(met_dss_file)

    # Create derived flow records used in CE-QUAL-W2
    compute_folsom_flows(currentAlternative, rtw, hydro_dss, output_dss_file)

    # Generate plotting-helper time series used by visualization scripts
    compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file)


def preprocess_ResSim_American(currentAlternative, computeOptions):
    '''
    Main preprocessing routine for the American River
    ResSim workflow.

    Responsibilities:

        1. Identify project and DSS file locations.
        2. Standardize DMS units and record types.
        3. Convert water temperatures to Celsius.
        4. Create standard zero/one DSS template records
           used by model logic.
        5. Generate derived Folsom flow records.
        6. Generate plotting records.

    Creates supporting records such as:
        - Zero flow series
        - Zero temperature series
        - Zero gate series
        - One gate series
        - Zero evaporation series

    Returns:
        True on successful completion.
    '''
    
    # Retrieve input DSS filename (main input source for ResSim)
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    # Determine run and project locations for assembling shared DSS paths
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    
    # Model timestep and location of shared storage directory
    balance_period = currentAlternative.getTimeStep()
    shared_dir = os.path.join(project_dir, 'shared')

    # Output DSS file used as the preprocessing data store
    output_dss_file = os.path.join(shared_dir,'DMS_American_Pre-Process.dss')

    # Identify hydro and meteorological DMS files and standardize them
    hydro_dss = os.path.join(shared_dir, 'DMS_AmericanHydroTS.dss')
    fix_DMS_types_units(hydro_dss)
    met_dss_file = os.path.join(shared_dir,'DMS_AmericanMet.dss')
    fix_DMS_types_units(met_dss_file)
    
    # Convert all temperature records to Celsius for consistency
    standardize_bc_temp_water_to_C(hydro_dss,output_dss_file)

    # Create reusable constant zero/one time series (gates, evaporation, flows)
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='temp-water', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0, what='gate', 
                        dss_type='INST-VAL', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=1, what='gate', 
                        dss_type='INST-VAL', period='1HOUR',cpart='ONES',fpart='ONES')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='evap', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')

    # Generate Folsom inflow/outflow derived flow records
    compute_folsom_flows(currentAlternative, rtw, hydro_dss, output_dss_file)

    # Build visualization/diagnostic DSS series for QA and plotting
    compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file)

    # Indicate successful preprocessing
    return True


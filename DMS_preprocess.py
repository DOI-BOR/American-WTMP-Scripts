
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
    recs = DSS_Tools.get_sanitized_record_list(dss_file)

    dss = HecDss.open(dss_file)
    
    for r in recs:
        rlow = r.lower()
        # things not to read: paired data, integer/scalar/text vars and some
        # other things that are causing trouble.
        if not '/location info' in rlow and not '/temp-equil' in rlow and \
          not '/depth-temp' in rlow and not 'icpathsmap' in rlow and \
          not '/downstream_control_loc' in rlow and not 'temp-water-target' in rlow:
        
            tsc = dss.get(r,True)

            if "/flow" in rlow or "/1day/" in rlow:
                if not "/storage" and not "/stor" in rlow:
                    tsc.type = 'PER-AVER'
                    #tsc.setStoreAsDoubles(True)
                    print('FixDMS Write original name: ',rlow)
                    print('FixDMS Write output name: ',tsc.fullName)
                    dss.write(tsc)

            units = str(tsc.units).lower()  # just to make sure
            
            if units in units_need_fixing:
                if units == 'tenths':
                    # save off a copy of cloud record in 0-1 for ResSim
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-FRAC'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'FRAC'
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 10.0
                    #tsc.setStoreAsDoubles(True)         
                    print('FixDMS Write original name: ',rlow)
                    print('FixDMS Write output name: ',tsc.fullName)
                    dss.write(tsc)
                if units == 'radians':
                    # save off a copy in deg
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-DEG'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'deg'
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / (2*3.141592653589793) * 360.0
                    #tsc.setStoreAsDoubles(True)         
                    print('FixDMS Write original name: ',rlow)
                    print('FixDMS Write output name: ',tsc.fullName)
                    dss.write(tsc)
                if units == 'deg':
                    # save off a copy in redians
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-RADIANS'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'radians'
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 360.0 * (2*3.141592653589793)
                    #tsc.setStoreAsDoubles(True)            
                    print('FixDMS Write original name: ',rlow)
                    print('FixDMS Write output name: ',tsc.fullName)
                    dss.write(tsc)
                if units == 'kph':
                    # convert to m/s 
                    tsc.units = 'm/s'
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 3.6
                    #tsc.setStoreAsDoubles(True)
                    print('FixDMS Write original name: ',rlow)
                    print('FixDMS Write output name: ',tsc.fullName)
                    dss.write(tsc)
                if units == 'fract':
                    # save off a copy of cloud record in 0-1 for ResSim, with proper naming, reset orignial to tenths
                    original_fullName = tsc.fullName
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-FRAC'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'FRAC'
                    print('FixDMS Write original name: ',rlow)
                    print('FixDMS Write output name: ',tsc.fullName)
                    dss.write(tsc)
                
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] * 10.0
                    tsc.units = 'tenths'
                    tsc.fullName = original_fullName
                    print('FixDMS Write original name: ',rlow)
                    print('FixDMS Write output name: ',tsc.fullName)
                    dss.write(tsc)

    dss.close()


def fix_DMS_types_units_old(dss_file):
    '''This method was implemented to change data types to PER-AVER that are not coming from the DMS that way'''
    recs = DSS_Tools.get_sanitized_record_list(dss_file)
    dss = HecDss.open(dss_file)
    
    for r in recs:
        rlow = r.lower()
        if not '/location info' in rlow and not '/temp-equil' in rlow and \
          not '/depth-temp' in rlow and not 'icpathsmap' in rlow and \
          not '/downstream_control_loc' in rlow:
        
            tsm = dss.read(r)

            if "/flow" in rlow or "/1day/" in rlow:
                tsm.setType('PER-AVER')
                dss.write(tsm)
            
            if tsm.getUnits().lower() in units_need_fixing:
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
                if tsm.getUnits().lower() == 'fract':
                    # save off a copy of cloud record in 0-1 for ResSim, with proper naming, reset orignial to tenths
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
                if tsm.getUnits() == 'kph':
                    # convert to m/s 
                    tsc = tsm.getData()
                    tsc.units = 'm/s'
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 3.6
                    dss.write(tsc)

                    # also, add w2link
                    #rec_parts = tsc.fullName.split('/')
                    #if not "w2link" in rec_parts[3].lower():
                    #    rec_parts[3] += '-W2link'
                    #    tsc.fullName = '/'.join(rec_parts)
                    #    for i in range(len(tsc.values)) :
                    #        tsc.values[i] = tsc.values[i] / 3.6
                    #    dss.write(tsc)

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
    dss.close()


def standardize_bc_temp_water_to_C(dss_file,output_dss_file):
    '''Make copies of temp-water records in C (standardizing on C) for ResSim linking'''
    recs = DSS_Tools.get_sanitized_record_list(dss_file)
    dss = HecDss.open(dss_file)

    if dss_file == output_dss_file:
        dss_out = dss
    else:
        dss_out = HecDss.open(output_dss_file)
    
    for r in recs:
        rlow = r.lower()
        if '/temp-water' in rlow:

            tsc = dss.get(r,True)

            incoming_units = tsc.units.lower()
        
            tsc = dss.get(r,True)
            rec_parts = tsc.fullName.split('/')
            rec_parts[3] += '-C'
            tsc.fullName = '/'.join(rec_parts)
            tsc.units = 'C'
                        
            if incoming_units == 'f' or incoming_units == 'degf':                
                for i in range(len(tsc.values)) :
                    tsc.values[i] = (tsc.values[i] - 32.0)*5.0/9.0             

            dss_out.put(tsc)

    dss.close()
    if dss_file != output_dss_file:
        dss_out.close()


def DMS_fix_units_types(hydro_dss,met_dss_file):
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

    if isinstance(x, list):
        return [interp(point, xp, fp, left, right) for point in x]
    else:
        if left is None:
            left = fp[0]
        if right is None:
            right = fp[-1]

        if x < xp[0]:
            return left
        elif x > xp[-1]:
            return right
        else:
            for i in range(len(xp) - 1):
                if x >= xp[i] and x <= xp[i+1]:
                    # Perform the linear interpolation
                    t = (x - xp[i]) / (xp[i+1] - xp[i])
                    return fp[i] + t * (fp[i+1] - fp[i])


def interp_monthly_coeff(coeffs):
    '''jython sucks'''
    month_midpoints = [16, 45, 75, 105, 136, 166, 197, 228, 259, 289, 320, 350]    
    month_midpoints_hours = [(month_midpoints[i]-1)*24 for i in range(12)]
    for i in  [1,3,5,7,8,10,12]:
        month_midpoints_hours[i-1] += 12
    month_midpoints_hours_padded = [-15.5*24.0,] + month_midpoints_hours + [8784.0 + 15.5*24.0,]
    coeffs_padded = [coeffs[-1],] + coeffs + [coeffs[0],]
    hourly_coeff = []
    for i in range(8784):
        hourly_coeff.append(interp(i,month_midpoints_hours_padded,coeffs_padded))
    return hourly_coeff

def american_NF_temp_array_hourly(hour, NF_cms, MF_cms, T_air):
    '''CARDNO/Stantec North Fork American water temperature regression into Folsom
    returns degrees C'''
    NF_coeff = ([[ 3.774,  5.013,  7.568, 13.929, 19.23 , 22.008, 27.481, 26.076, 19.876, 11.463,  7.827,  3.52],
       [ 1.266,  2.088,  3.042,  1.493, -4.149, -2.19 ,  0.461, -0.056, -2.334,  0.665,  0.685,  -0.27],
       [-0.123, -2.308, -4.644, -5.956, -2.651, -4.32 , -8.106, -7.756, -4.285, -2.909, -1.342, 1.59],
       [ 0.209,  0.289,  0.336,  0.278,  0.279,  0.182,  0.071,  0.064, 0.107,  0.355,  0.367,  0.30]])
    coeff0 = interp_monthly_coeff(NF_coeff[0])
    coeff1 = interp_monthly_coeff(NF_coeff[1])
    coeff2 = interp_monthly_coeff(NF_coeff[2])
    coeff3 = interp_monthly_coeff(NF_coeff[3])
   
    nf_temp = []
    for i in range(len(hour)):
        h_idx = hour[i]
        print(hour[i],NF_cms[i],MF_cms[i],T_air[i])
        nf_temp.append( coeff0[h_idx] + coeff1[h_idx] * math.log10(NF_cms[i]) + coeff2[h_idx] * math.log10(MF_cms[i]) + coeff3[h_idx] * T_air[i] )        
        #print(hour[i],NF_cms[i],MF_cms[i],T_air[i],nf_temp[-1])
    return nf_temp

def american_SF_temp_array_hourly(hour, SF_cms, T_air):
    '''CARDNO/Stantec South Fork American water temperature regression into Folsom
    returns degrees C'''
    SF_coeff = ([[ 1.956,  3.894,  8.456, 12.605, 19.374, 22.03 , 23.604, 21.761, 17.663, 11.832,  6.521,  3.430],
       [ 1.374,  0.221, -1.422, -3.05 , -5.815, -6.605, -5.623, -5.196, -4.067, -2.665, -0.366,  0.755],
       [ 0.29 ,  0.282,  0.224,  0.223,  0.204,  0.216,  0.114,  0.105, 0.155,  0.299,  0.374,  0.358 ]])
    coeff0 = interp_monthly_coeff(SF_coeff[0])
    coeff1 = interp_monthly_coeff(SF_coeff[1])
    coeff2 = interp_monthly_coeff(SF_coeff[2])
       
    sf_temp = []
    for i in range(len(hour)):
        h_idx = hour[i]
        sf_temp.append( coeff0[h_idx] + coeff1[h_idx] * math.log10(SF_cms[i]) + coeff2[h_idx] * T_air[i] )
        #print(month[i],SF_cms[i],T_air[i],sf_temp[-1])
    return sf_temp

            
def american_NF_temp_array(month, NF_cms, MF_cms, T_air):
    '''CARDNO/Stantec North Fork American water temperature regression into Folsom
    returns degrees C'''
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
    nf_temp = []
    for i in range(len(month)):
        coeff = NF_coeff[month[i]-1]
        nf_temp.append( coeff[0] + coeff[1] * math.log10(NF_cms[i]) + coeff[2] * math.log10(MF_cms[i]) + coeff[3] * T_air[i] )
        #print(month[i],NF_cms[i],MF_cms[i],T_air[i],nf_temp[-1])
    return nf_temp

def american_SF_temp_array(month, SF_cms, T_air):
    '''CARDNO/Stantec South Fork American water temperature regression into Folsom
    returns degrees C'''
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
    sf_temp = []
    for i in range(len(month)):
        coeff = SF_coeff[month[i]-1]
        sf_temp.append( coeff[0] + coeff[1] * math.log10(SF_cms[i]) + coeff[2] * T_air[i] )
        #print(month[i],SF_cms[i],T_air[i],sf_temp[-1])
    return sf_temp

def american_SC_temp_array(month):
    '''CARDNO/Stantec South Canal monthly average inflow temperature into Folsom
    returns degrees C
    '''
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
    sc_temp = []
    for i in range(len(month)):
        sc_temp.append( (SC_ave_temp[month[i]-1] - 32.0)*5.0/9.0 )
    return sc_temp

def calc_folsom_inflow_temps(currentAlt, rtw, hydro_dss, met_dss_file, output_dss_file,hourly=True):
    '''This is currently not needed - using measured temps fro DMS.
    This is was an implementation of the forecast inflow temp regressions - keeping this code aroud because it took a long time to write
    in jython, and might be useful.  The functions called below can do all the regression coefficient interpolations, etc.
    '''
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()    
    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    currentAlt.addComputeMessage('Calculating Folsom Inflow Temperatures...')

    # hardcoded paths ... yuck
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
    shared_path,_ = os.path.split(hydro_dss)
    American_inflows_6 = os.path.join(shared_path,'American_inflows_6.dss')
    NF = '::'.join([American_inflows_6,'/11427000/NF AMERICAN/FLOW//1HOUR/USGS-CARDNO-FROM-DAILY/'])
    MF = '::'.join([American_inflows_6,'/11433300/MF AMERICAN/FLOW//1HOUR/USGS-CARDNO-FROM-DAILY/'])
    SF = '::'.join([American_inflows_6,'/11444500/SF AMERICAN/FLOW//1HOUR/USGS-CARDNO-FROM-DAILY/'])
    AT = '::'.join([met_dss_file,'/MR Am.-Natoma Lake/Fair Oaks-Air Temp/Temp-Air//1hour/251.40.53.1.1/'])
    output_period = '1Hour'

    # Read inputs
    print('Reading inflows')
    inputs = []
    for i,input_rec in enumerate([NF,MF,SF,AT]):
        if '::' in input_rec:
            dss_file_alt,inflow_rec_alt = input_rec.split('::')
            dssFm_alt = HecDss.open(dss_file_alt)
            ts_data = dssFm_alt.read(inflow_rec_alt, starttime_str, endtime_str, False).getData()
            dssFm_alt.close()
            print(dss_file_alt)
        else:
            currentAlt.addComputeMessage('ERROR reading' + str(input_rec) + '[does it have form <filepath>::<dssrec>?]')
            sys.exit(-1)
        if i==3:
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
        hectimes = ts_data.times

    # generate month record needed for calcs
    time_post = HecTime(HecTime.MINUTE_INCREMENT)
    month = []
    day = []
    hour = []
    for time in hectimes:
        time_post.set(time)
        month.append(time_post.month())
        day.append(time_post.dayOfYear() + time_post.hour()/24.0)
        hour.append((time_post.dayOfYear()-1)*24 + time_post.hour())

    # find/write inflow temps
    dssFm_out = HecDss.open(output_dss_file)
    for i in range(3):
        if i==0: # NF
            outrec = '/Folsom/NF Inflow/Temp-Water//%s/ResSim_PreProcess/'%output_period
            #outdata = american_NF_temp_array(month, inputs[0], inputs[1], inputs[3])
            outdata = american_NF_temp_array_hourly(hour, inputs[0], inputs[1], inputs[3])
        elif i==1: # SF
            outrec = '/Folsom/SF Inflow/Temp-Water//%s/ResSim_PreProcess/'%output_period
            #outdata = american_SF_temp_array(month, inputs[2], inputs[3])
            outdata = american_SF_temp_array_hourly(hour, inputs[2], inputs[3])
        else: # SC
            outrec = '/Folsom/SC Inflow/Temp-Water//%s/ResSim_PreProcess/'%output_period
            outdata = american_SC_temp_array(month)    
        # Output record
        tsc = TimeSeriesContainer()
        tsc.times = hectimes
        tsc.fullName = outrec
        tsc.values = outdata
        #tsc.startTime = times[1]
        tsc.units = 'C'
        tsc.type = 'PER-AVER'
        tsc.numberValues = len(outdata)        
        dssFm_out.write(tsc)

    dssFm_out.close()


def compute_folsom_flows(currentAlternative, rtw, hydro_dss, output_dss_file):
    # Add Mormon Ravine and Newcastle PP flows (there is one inflow designed for these in ResSim)
    inflow_records = ['/MR Am.-Folsom Lake/11425416 Newcastle PP-Daily Flow/Flow//1Day/250.114.125.1.1/',
                      '/MR Am.-Folsom Lake/11433930 Mormon Ravine-Daily Flow/Flow//1Day/250.115.125.1.1/']
    DSS_Tools.add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Am.-Folsom Lake/MormonR_NewcastlePP_Sum/Flow//1Day/ResSim_PreProcess/', output_dss_file)

    # total North Arm Folsom inflow = NF (Lake Clementine) + MF (Foresthill) - ARPS (upsample monthly A. R. Pump stations) 

    # this mothnly record is a pain in the butt. the resampling always leaves off the start and end because the read typically excludes them,
    # and the DMS seems to be missing Dec 2021 at this point. so we expand dates but also end Nov 30 in 2021
    DSS_Tools.resample_dss_ts(hydro_dss,'/MR Am.-Folsom Lake/American River Pump Station (ARPS)-Flow/Flow//1Mon/250.401.125.1.1/',
                              rtw,output_dss_file,'1Day',pad_1mon=True)
    out_rec = '/MR Am.-Folsom Lake/North Arm/Flow//1Day/ResSim_PreProcess/'
    shared_path,_ = os.path.split(hydro_dss)
    American_inflows_6 = os.path.join(shared_path,'American_inflows_6.dss')
    flow_Records = ['/MR Am.-Folsom Lake/11427000 Lake Clementine Dam-Daily Flow/Flow//1Day/250.112.125.1.1/',
                    '/MR Am.-Folsom Lake/11433300 Foresthill-Daily Flow/Flow//1Day/250.113.125.1.1/',
                    #'/MR Am.-Folsom Lake-11433300 Foresthill-Daily Flow/Flow//1Day/250.113.125.1.1/',  # this one is not downloading yet
                    #'::'.join([American_inflows_6,'/11433300/MF AMERICAN/FLOW//1DAY/USGS-CARDNO/']),
                    '::'.join([output_dss_file,'/MR Am.-Folsom Lake/American River Pump Station (ARPS)-Flow/Flow//1Day/250.401.125.1.1/']),]
                    #'::'.join([American_inflows_6,'/FOLSOM/AMERICAN RIVER PUMP STATION/FLOW//1DAY/CDEC-CARDNO/']),] # this one is not downloading yet
                                                                                     # [N.A.,ADD,SUBTRACT]
    DSS_Tools.add_or_subtract_flows(currentAlternative, rtw, flow_Records, hydro_dss, [None,True,False], out_rec, output_dss_file)
        
    # EID outflow - do we need to normalize to daily?
    # /MR Am.-Folsom Lake-EID Folsom Diversion-Diversion Flow/Flow/ --?

    # Sum Folsom Dam river outlets - Upper elevation
    out_rec = '/MR Am.-Folsom Lake/Upper_River_Outlets_Sum_min4/Flow//1Hour/ResSim_PreProcess/'
    outflow_records = ['/MR Am.-Folsom Lake/FOL-Outlet Release G1/Flow//1Hour/250.3.125.23.1/',
                       '/MR Am.-Folsom Lake/FOL-Outlet Release G2/Flow//1Hour/250.3.125.24.1/',
                       '/MR Am.-Folsom Lake/FOL-Outlet Release G3/Flow//1Hour/250.3.125.25.1/',
                       '/MR Am.-Folsom Lake/FOL-Outlet Release G4/Flow//1Hour/250.3.125.26.1/',]
    DSS_Tools.add_flows(currentAlternative, rtw, outflow_records, hydro_dss, out_rec, output_dss_file)
    DSS_Tools.min_ts(output_dss_file, out_rec, 4.0, output_dss_file, 'ResSim_PreProcess')

    # Sum Folsom Dam river outlets - Lower elevation
    out_rec = '/MR Am.-Folsom Lake/Lower_River_Outlets_Sum_min4/Flow//1Hour/ResSim_PreProcess/'
    outflow_records = ['/MR Am.-Folsom Lake/FOL-Outlet Release G5/Flow//1Hour/250.3.125.27.1/',
                       '/MR Am.-Folsom Lake/FOL-Outlet Release G6/Flow//1Hour/250.3.125.28.1/',
                       '/MR Am.-Folsom Lake/FOL-Outlet Release G7/Flow//1Hour/250.3.125.29.1/',
                       '/MR Am.-Folsom Lake/FOL-Outlet Release G8/Flow//1Hour/250.3.125.30.1/',]
    DSS_Tools.add_flows(currentAlternative, rtw, outflow_records, hydro_dss, out_rec, output_dss_file)
    DSS_Tools.min_ts(output_dss_file, out_rec, 4.0, output_dss_file, 'ResSim_PreProcess')

    # Sum Natoma Gen releases
    out_rec = '/MR Am.-Natoma Lake/NAT-Gen Release Sum/Flow//1Hour/ResSim_PreProcess/'
    outflow_records = ['/MR Am.-Natoma Lake/NAT-Gen Release U1/Flow//1Hour/251.4.125.3.1/',
                       '/MR Am.-Natoma Lake/NAT-Gen Release U2/Flow//1Hour/251.4.125.4.1/']
    DSS_Tools.add_flows(currentAlternative, rtw, outflow_records, hydro_dss, out_rec, output_dss_file)


def compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file):
    pass


def preprocess_W2_American(currentAlternative, computeOptions):
    rtw = computeOptions.getRunTimeWindow()
    
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    balance_period = currentAlternative.getTimeStep()
    shared_dir = os.path.join(project_dir, 'shared')

    output_dss_file = os.path.join(shared_dir,'DMS_American_Pre-Process.dss') 

    hydro_dss = os.path.join(shared_dir, 'DMS_AmericanHydroTS.dss')
    fix_DMS_types_units(hydro_dss)
    met_dss_file = os.path.join(shared_dir,'DMS_AmericanMet.dss')
    fix_DMS_types_units(met_dss_file)

    compute_folsom_flows(currentAlternative, rtw, hydro_dss, output_dss_file)

    compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file)


def preprocess_ResSim_American(currentAlternative, computeOptions):
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    balance_period = currentAlternative.getTimeStep()
    shared_dir = os.path.join(project_dir, 'shared')

    output_dss_file = os.path.join(shared_dir,'DMS_American_Pre-Process.dss')

    hydro_dss = os.path.join(shared_dir, 'DMS_AmericanHydroTS.dss')
    fix_DMS_types_units(hydro_dss)
    met_dss_file = os.path.join(shared_dir,'DMS_AmericanMet.dss')
    fix_DMS_types_units(met_dss_file)
    standardize_bc_temp_water_to_C(hydro_dss,output_dss_file)

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

    # if template IDs exist still, remove them
    #DSS_Tools.strip_templateID_and_rename_records(hydro_dss,currentAlternative)
    #DSS_Tools.strip_templateID_and_rename_records(met_dss_file,currentAlternative)

    compute_folsom_flows(currentAlternative, rtw, hydro_dss, output_dss_file)

    # resampling for inflow temp calc
    #sdf.resample_dss_ts(met_dss_file,'/MR Am.-Natoma Lake/Fair Oaks-Air Temp/Temp-Air//1Hour/251.40.53.1.1/',rtw,output_dss_file,'1DAY')

    #calc_folsom_inflow_temps(currentAlternative, rtw, hydro_dss, met_dss_file, output_dss_file)   

    compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file)


    return True


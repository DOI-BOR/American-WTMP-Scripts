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
import os, sys, csv, fileinput

from com.rma.io import DssFileManagerImpl
from com.rma.model import Project
#import hec.hecmath.TimeSeriesMath as tsmath

import DSS_Tools
reload(DSS_Tools)

import DMS_preprocess
reload(DMS_preprocess)

import equilibrium_temp
reload(equilibrium_temp)

import create_balance_flow_jython as cbfj
reload(cbfj)

from com.rma.io import DssFileManagerImpl
from java.util import TimeZone

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


def interp_monthly_coeff_daily(monthly_coeffs):
    '''jython sucks'''
    month_midpoints = [16.5, 46, 75.5, 106, 136.5, 167, 197.5, 228.5, 259, 289.5, 320, 350.5]    
    month_midpoints_padded = [-14.5] + month_midpoints + [366.0+15.5]
    coeffs_padded = [monthly_coeffs[-1],] + monthly_coeffs + [monthly_coeffs[0],]
    daily_coeff = []
    for i in range(1,367):
        daily_coeff.append(interp(i,month_midpoints_padded,coeffs_padded))
    return daily_coeff

def eq_temp(rtw,at,cl,ws,sr,td,eq_temp_out):
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()

    # get at data and times in the formats needed
    dssFm = HecDss.open(at[0])        
    tsc = dssFm.read(at[1], starttime_str, endtime_str, False).getData()
    tsc_int_times = tsc.times
    dtt = DSS_Tools.hectime_to_datetime(tsc)
    at_data = tsc.values
    dssFm.close()
    # get the rest of the data over the same period
    cl_data = DSS_Tools.data_from_dss(cl[0],cl[1],starttime_str,endtime_str)
    ws_data = DSS_Tools.data_from_dss(ws[0],ws[1],starttime_str,endtime_str)
    sr_data = DSS_Tools.data_from_dss(sr[0],sr[1],starttime_str,endtime_str)
    td_data = DSS_Tools.data_from_dss(td[0],td[1],starttime_str,endtime_str)

    for i in range(48):
        print(dtt[i].strftime('%Y-%m-%d %H:%M'),at_data[i],cl_data[i],sr_data[i],td_data[i],ws_data[i])
    
    Te = equilibrium_temp.calc_equilibrium_temp(dtt,at_data,cl_data,sr_data,td_data,ws_data)
    (dtt, at, cl, sr, td, ws)
    print('writing: ',eq_temp_out[1])
    tsc = TimeSeriesContainer()
    tsc.times = tsc_int_times
    tsc.fullName = eq_temp_out[1]
    tsc.values = Te
    tsc.units = 'C'
    tsc.type = 'INST-VAL'
    tsc.numberValues = len(tsc.values)

    tsm = tsmath(tsc)
    tsm_day = DSS_Tools.standardize_interval(tsm,'1day')
    tsm_wk = DSS_Tools.standardize_interval(tsm,'1week')
        
    dssFmOut = HecDss.open(eq_temp_out[0])
    dssFmOut.write(tsc)
    dssFmOut.write(tsm_day)
    dssFmOut.write(tsm_wk)
    dssFmOut.close()

def read_temp_schedule_csv(csv_file_path):
    ''' Expecting file of this format:
    Sch#,May,June,July,Aug,Sept,Oct,Nov
    1,63,63,63,63,63,56,56
    2,63,63,63,63,63,57,56
    3,63,63,63,63,63,58,56
    4,63,63,63,63,63,59,56
    5,63,63,63,63,63,60,56
    ....
    
    returns a list of lists of temps by row
    '''
    schedReader = csv.reader(open(csv_file_path), delimiter=',', quotechar='|')
    monthlyTempsSched = []
    for i,row in enumerate(schedReader):
        if i>0: # header row?
            monthlyTempsSched.append(row)
    return monthlyTempsSched


def get_interpolated_coeffs(location):
    # regression parameters for Watt Ave temperature target location -> Folsom release temp - location == 1
    watt_coeffs = [[ 1.81876287,  3.53936881,  4.98729431,  8.85508504, 11.86438162,
        12.62006676, 13.08970612, 16.90056672, 15.26051768,  2.75799392,
        -1.81979657, -1.3513815 ],
       [ 0.11264133,  0.20514104,  0.1686952 ,  0.15834058,  0.13853404,
         0.0864804 ,  0.04903162,  0.07265319,  0.18274807,  0.24001519,
         0.23128533,  0.16759159],
       [ 0.73259158,  0.77538058,  0.9613015 ,  0.7189698 ,  0.74203971,
         0.80089022,  0.682267  ,  0.62275869,  0.454786  ,  0.6680434 ,
         0.82870536,  0.68563562],
       [-0.14771242, -1.49933353, -2.51250945, -3.03473477, -4.37180174,
        -4.43688051, -3.29314944, -5.22242857, -4.50893573, -0.34458292,
         0.70364452,  1.78619426],
       [ 0.57283027,  0.4516801 ,  0.9059361 ,  0.90501824,  0.9544556 ,
         0.94857614,  0.93717274,  0.91059113,  0.83642545,  0.81923008,
         0.95854366,  0.87062062]]

    # regression parameters for Havel Ave temperature target location -> Folsom release temp - location == 2
    hazel_coeffs = [
        [ 1.94637876,  3.06584481,  4.84020194,  8.04959568,  7.50079909, 9.25305813,  7.07595543,  9.12456928,  9.76617634,  9.75723367, 8.16619157,  1.53895622],
        [ 0.01337584,  0.06344339,  0.03040309,  0.0490183 ,  0.03563944, 0.02569569, -0.0156343 ,  0.03155582,  0.08446625,  0.05954057, 0.19014651,  0.02926542],
        [ 0.93981272,  0.80547868,  0.84083274,  0.60822924,  0.83867899, 0.79682248,  0.88547439,  0.76824479,  0.69468425,  0.59837608, 0.5110841 ,  0.89079026],
        [-0.61820524, -0.82588954, -1.53055365, -1.78212831, -2.49025734, -2.9129365 , -1.9228547 , -2.65828672, -3.27316002, -2.16631531, -1.814512  , -0.22356768]
    ]

    # TODO: regression parameters for RiverMile, location == 3
    rivermile_coeffs = [
        [ 2.23549682e+00,  3.56028367e+00,  6.52712530e+00, 1.02213747e+01,  1.05801274e+01,  1.42266170e+01, 1.36368749e+01,  1.66864436e+01,  1.49065332e+01, 8.96273105e+00,  7.70645096e+00,  9.70660073e-01],
        [ 5.54179380e-02,  1.36863661e-01,  9.33124620e-02, 1.00763965e-01,  1.09121502e-01,  4.82506310e-02, 1.90429540e-02,  6.24275800e-02,  1.33762465e-01, 1.33697106e-01,  2.69051024e-01,  9.57283810e-02],
        [ 8.25737722e-01,  8.11085616e-01,  8.58480132e-01, 6.27146407e-01,  8.45379413e-01,  7.53645265e-01, 7.60110529e-01,  6.14723949e-01,  5.48927681e-01, 5.69276475e-01,  4.53617803e-01,  8.36113881e-01],
        [-5.26446885e-01, -1.09530356e+00, -1.95242061e+00, -2.31341337e+00, -3.40187678e+00, -3.97142722e+00, -3.24758750e+00, -4.15399187e+00, -4.02447280e+00, -1.72353780e+00, -1.73660858e+00, -1.60383653e-01],
        [ 8.75064800e-03, -3.64423800e-02, -8.41306680e-02, -9.34771460e-02, -1.22967346e-01, -1.24829172e-01, -1.27904193e-01, -1.39371934e-01, -1.10916731e-01, -3.13301590e-02,  8.45384500e-03,  2.25493060e-02]
    ]

    if location==1: # Watt Ave
        coeffs = watt_coeffs
    elif location==3: # Hazel Ave
        coeffs = hazel_coeffs
    else:
        raise ValueError("Forecast Preprocess: W2 downstream target not supported:"+str(location))
        
    # TODO: needs work as the form of the regression is different below
    #elif location==3: # RiverMile
    #    coeffs = rivermile_coeffs
    #    cf4 = interp_monthly_coeff_daily(coeffs[4]) # x4

    cf0 = interp_monthly_coeff_daily(coeffs[0]) # int(tercpet)
    cf1 = interp_monthly_coeff_daily(coeffs[1]) # x1
    cf2 = interp_monthly_coeff_daily(coeffs[2]) # x2
    cf3 = interp_monthly_coeff_daily(coeffs[3]) # x3

    return cf0,cf1,cf2,cf3


def write_target_temp_npt(year,location,doys,Tair,FaveFlow,schedule_csv,targt_temp_npt_filepath,lagWatt=True):
    '''    
    Tair [C] 
    FaveFlow [CMS]
  
    '''
    cf0,cf1,cf2,cf3 = get_interpolated_coeffs(location)

    monthlyTempsSched = read_temp_schedule_csv(schedule_csv)

    ReleaseTemp = []
    n99_line = [-99.0] * len(monthlyTempsSched)

    ReleaseTemp_debug = []

    dt_base = dt.datetime(year,1,1) + dt.timedelta(days=doys[0]-1)
    one_day = dt.timedelta(days=1)
    n_days = len(doys)
    dtt = [dt_base + one_day*i for i in range(n_days)]

    for i,date in enumerate(dtt):
        mon = date.month
        dayofyear = date.timetuple().tm_yday 
        if mon > 4 and mon < 12:          
            lag = 0
            if lagWatt and location==1:
                lag = int(3966.8*(35.314*FaveFlow[i])**(-0.944))
                lag = max(0,lag)
            dateLag = date + dt.timedelta(days=lag)
            # schedule runs from may -> nov, index for monthy schedule temp, but first value is jday...
            # example downstream target schedule:
            # jday,may,jun,jul,aug,sep,oct,nov            
            #   20, 65, 65, 65, 65, 65, 63, 58  <-- deg F, except jday
            mlag = dateLag.month - 4  # index to schedule list above, accounting for lag date and jday as first index
            mlag = min(mlag,7)
            ilag = min(i + lag, n_days-1)  # lag data index, but don't run past end of data
            d = dayofyear - 1 # day-of-year, minus 1, as index to coeffs

            rt = [dayofyear]
            for k,sched in enumerate(monthlyTempsSched):
                # Below is original code, which has some negatives that don't make a lot of sense? Copied to python exactly anyway
                # FORTRAN: ReleaseTemp(k,i)=(-((TTarg(month(k),i)-32)/1.8)+int(m)+x1(m)*aveTair(k)+x3(m)*log10(FaveFlow(k)))/(-x2(m))                
                # JYTHON:               (-1.0*((sched[mlag]-32.0)/1.8)+cf0[d] + cf1[d]*Tair[ilag] + cf3[d]*math.log10(FaveFlow[ilag]))/(-1.0*cf2[d])

                target_T = (float(sched[mlag]) - 32.0) / 1.8 # regression works in Celsius

                release_T =  (-1.0*target_T + cf0[d] + cf1[d]*Tair[ilag] + cf3[d] * math.log10(FaveFlow[ilag])) / (-1.0*cf2[d])

                # TODO: insert RiverMile calc here, if location==3

                if release_T > target_T or FaveFlow[ilag] > 300.0:
                    release_T = target_T

                rt.append(release_T)
                
                # Debug
                #if k==0:
                #   print(rt[-1])
                #   print(float(sched[mlag]),cf0[d],cf1[d],Tair[ilag],cf3[d],FaveFlow[ilag],cf2[d])
                #   print(mlag,d,ilag)
                
            ReleaseTemp.append(rt)            
            ReleaseTemp_debug.append(rt+[dayofyear,Tair[ilag],FaveFlow[ilag],cf0[d],cf1[d],cf2[d],cf3[d]])
        else:
            d = dayofyear - 1
            ReleaseTemp.append([dayofyear]+n99_line)
            ReleaseTemp_debug.append([dayofyear]+n99_line+[dayofyear,-99.0,-99.0,cf0[d],cf1[d],cf2[d],cf3[d]])

    # got an EOF error in FORTRAN ... adding extra day to see if it solves this
    ReleaseTemp.append([dayofyear+1]+n99_line)
    
    # Write to npt format
    with open(targt_temp_npt_filepath,'w') as fp:
        fp.write("   Temperature Target based release temps (May 1st - Nov 30) -- WTMP generated\n")
        for rt_values in ReleaseTemp:
            fp.write("%8i,"%rt_values[0])
            for rt in rt_values[1:]:
                fp.write("  %6.2f,"%rt)
            fp.write("\n")

    with open(targt_temp_npt_filepath+'.debug.csv','w') as fp:
        fp.write("   Upstream Target Temp debug file -- WTMP generated\n")
        for rt_values in ReleaseTemp_debug:
            fp.write("%8i,"%rt_values[0])
            for rt in rt_values[1:-7]:
                fp.write("  %6.2f,"%rt)
            fp.write("%8i,"%rt_values[-7])
            for rt in rt_values[-6:]:
                fp.write("  %.12e,"%rt)            
            fp.write("\n")

    return True


def calc_downstream_temp_W2(year,location,doys,Tair,Tmodel,FaveFlow,lagWatt=False):
    '''
    if Tmodel is set (to the modeled Folson release temperature record), that indicates we are doing a reverse regression 
    back downstream, using model results.  Otherwise, a Folsom release temperature is calculated from the schedule temperatures.    
    
    Tair [C] 
    FaveFlow [CMS]

    REVERSE REGRESSIONS--------------------------------------------------
        WattTemp(k)=intW(m)+x1W(m)*aveTair(k)+x2W(m)*avetemp1(k-Lag)+x3W(m)*log10(FaveFlow(k))  
        HazelTemp(k)=intH(m)+x1H(m)*aveTair(k)+x2H(m)*avetemp1(k)+x3H(m)*log10(FaveFlow(k))
        RiverMileTemp(k)=intRM(m)+x1RM(m)*aveTair(k)+x2RM(m)*avetemp1(k)+x3RM(m)*log10(FaveFlow(k))+x4RM(m)*RM
    
    '''
    cf0,cf1,cf2,cf3 = get_interpolated_coeffs(location)

    DownstreamTemp = []

    dt_base = dt.datetime(year,1,1) + dt.timedelta(days=doys[0]-1)
    one_day = dt.timedelta(days=1)
    n_days = len(doys)
    dtt = [dt_base + one_day*i for i in range(n_days)]

    for i,date in enumerate(dtt):
        mon = date.month
        dayofyear = date.timetuple().tm_yday 
        
        lag = 0
        if lagWatt and location==1:
            lag = int(3966.8*(35.314*FaveFlow[i])**(-0.944))
            lag = max(0,lag)
        dateLag = date + dt.timedelta(days=lag)
        ilag = min(i + lag, n_days-1)  # lag data index, but don't run past end of data
        d = dayofyear - 1 # day-of-year, minus 1, as index to coeffs        
        iBackLag = max(0,i-lag) # needed for reverse regression

        # below is original reverse regression code of Watt Ave:
        # avetemp1 is the daily model result temperature at Watt (Tmodel)
        # FORTRAN: WattTemp(k)=intW(m)+x1W(m)*aveTair(k)+x2W(m)*avetemp1(k-Lag)+x3W(m)*log10(FaveFlow(k))
        # JYTHON:              cf0[d] + cf1[d]*Tair[i] + cf2[d]*Tmodel[iBackLag] + cf3[d] * math.log10(FaveFlow[i])
        backTemp = cf0[d] + cf1[d]*Tair[i] + cf2[d]*Tmodel[iBackLag] + cf3[d] * math.log10(FaveFlow[i])

        # debug
        #print(i,d,cf0[d],cf1[d],cf2[d],cf3[d],Tair[i],Tmodel[iBackLag],FaveFlow[i],backTemp)

        # TODO: insert RiverMile calc here, if location==3

        # no regression at high water
        if FaveFlow[i] > 300.0:
            backTemp = Tmodel[iBackLag]
            
        DownstreamTemp.append(backTemp)

    return dtt,DownstreamTemp



def update_W2_Folsom_iterative_restart_date(rtw,model_dir):    
    starttime_str = rtw.getStartTimeString()
    doy = HecTime(starttime_str).dayOfYear()
    restart_doy = 120 if doy < 120 else doy + 1 # write restart file either May 1 or start day + 1

    # fileinput should backup and write over file, with inplace=True

    # replace doy in line 403 (index 402) in w2_con.csv
    w2_con = os.path.join(model_dir,"w2_con.csv")
    for line in fileinput.input(w2_con, inplace=True):
        if fileinput.filelineno()==403:
            print("%i,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,"%restart_doy)
        else:
            print("%s" % line.rstrip())

    # replace doy in line 4 in folsom_in.npt
    folsom_in = os.path.join(model_dir,"folsom_in.npt")
    for line in fileinput.input(folsom_in, inplace=True):
        if fileinput.filelineno()==4:
            print("%i,,,,,,,,,,,,,,,,,,"%restart_doy)
        else:
            print("%s" % line.rstrip())


def update_W2_Folsom_iterative_schedule_number(run_dir,model_dir):
    with open(os.path.join(run_dir,'current_ensemble.txt'),'r') as fp:
        ensemble_number = int(fp.readline().strip())

    # replace doy in line 2 of "BestGuessTTS.csv"
    best_guess_file = os.path.join(model_dir,"BestGuessTTS.csv")
    with open(best_guess_file,'w') as fp:
        fp.write("Year,BG,CB\n1,%i,10\n"%ensemble_number)


def storage_to_elev(res_name,elev_stor_area,forecast_dss,storage_rec,conic=False):
    dssFmRec = HecDss.open(forecast_dss)
    tsc = dssFmRec.get(storage_rec,True) # read ALL data in record

    elev = []
    if conic:
        print('Conic interpolation of elevations from storage not supported yet.')
        sys.exit(-1)
    else:
        for j in range(tsc.numberValues):
            elev.append(cbfj.linear_interpolation(elev_stor_area['stor'], elev_stor_area['elev'], tsc.values[j]))
            print('stor2elev: ',j,tsc.times[j],tsc.values[j],elev[-1])

    #if tsc.type.lower != 'inst-val':
    #    # going from PER-CUM -> INST-VAL, so need to move the times up one to account for the end-of-period reporting of PER-CUM
    #    time_delta = tsc.times[1] - tsc.times[0]
    #    times = [tsc.times[0]-time_delta]
    #    for j in range(tsc.numberValues-1):
    #        times.append(tsc.times[j])
    #    tsc.times = times

    recparts = tsc.fullName.split('/')
    recparts[2] = res_name
    recparts[3] = 'ELEV'
    tsc.fullName = '/'.join(recparts)
    tsc.units = 'ft'
    tsc.type = 'INST-VAL'
    tsc.values = elev
    dssFmRec.write(tsc)
    dssFmRec.close()

def invent_elevation(res_name,forecast_dss,storage_rec,elev_constant_ft):
    dssFmRec = HecDss.open(forecast_dss)
    tsc = dssFmRec.get(storage_rec,True)
    recparts = tsc.fullName.split('/')
    recparts[2] = res_name
    recparts[3] = 'ELEV'
    tsc.fullName = '/'.join(recparts)
    tsc.units = 'ft'
    tsc.type = 'INST-VAL'
    tsc.values = [elev_constant_ft for j in range(tsc.numberValues)]

    #if tsc.type.lower != 'inst-val':
    #    # going from PER-CUM -> INST-VAL, so need to move the times up one to account for the end-of-period reporting of PER-CUM
    #    time_delta = tsc.times[1] - tsc.times[0]
    #    times = [tsc.times[0]-time_delta]
    #    for j in range(tsc.numberValues-1):
    #        times.append(tsc.times[j])
    #    tsc.times = times
    
    dssFmRec.write(tsc)
    dssFmRec.close()

def write_forecast_elevations(currentAlternative, rtw, forecast_dss, shared_dir):

    # get date for starting elevation - look for end-of-month before start time
    ht = HecTime(rtw.getStartTimeString())
    if ht.month() == 1:
        start_dt = dt.datetime(ht.year()-1,12,31)
        start_str = start_dt.strftime('%d%b%Y')+ ' 2400'
    else:
        start_dt = dt.datetime(ht.year(),ht.month(),1)
        start_dt = start_dt - dt.timedelta(days=1)
        start_str = start_dt.strftime('%d%b%Y')+ ' 2400'
    end_str = rtw.getEndTimeString()   

    currentAlternative.addComputeMessage('Forecast Elevations: '+start_str+' '+end_str)

    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_Folsom.csv'), 'Folsom')

    # convert Folsom storage to elevation
    storage_to_elev('Folsom',elev_stor_area,forecast_dss,'//FOLSOM/STORAGE//1Month/AMER_BC_SCRIPT/',conic=False)
    storage_to_elev('Folsom',elev_stor_area,forecast_dss,'/AMERICAN RIVER/FOLSOM LAKE/STORAGE-CVP//1Day/AMER_BC_SCRIPT/',conic=False)
    
    # write an hourly forecast elevation based on starting elevation and flows
    DSS_Tools.resample_dss_ts(forecast_dss,'//FOLSOM/FLOW-RELEASE//1Hour/AMER_BC_SCRIPT/',None,forecast_dss,'1DAY')
    inflow_records = ['//Folsom-NF-in/FLOW-IN//1Day/AMER_BC_SCRIPT/',
                      '//Folsom-SF-in/FLOW-IN//1Day/AMER_BC_SCRIPT/',
                      '/AMERICAN RIVER/FOLSOM LAKE/FLOW-ACC-DEP//1Day/AMER_BC_SCRIPT/']  # this actually evap, but negative already, so it goes as inflow
    outflow_records = ['//FOLSOM/FLOW-RELEASE//1Day/AMER_BC_SCRIPT/']
    starting_elevation = DSS_Tools.first_value(forecast_dss,'//Folsom/ELEV//1Month/AMER_BC_SCRIPT/',start_str,end_str)
    print('starting_elevation ',starting_elevation)

    elev_calc_start_dt = start_dt + dt.timedelta(days=1)
    elev_calc_start_str = elev_calc_start_dt.strftime('%d%b%Y')+ ' 2400'
    #cbfj.predict_elevation(currentAlternative, rtw, 'Folsom Lake', inflow_records, outflow_records, starting_elevation,
    cbfj.predict_elevation(currentAlternative, start_str, end_str, 'Folsom Lake', inflow_records, outflow_records, starting_elevation,
                         elev_stor_area, forecast_dss, '//Folsom Lake/ELEV-FORECAST//1DAY/AMER_BC_SCRIPT/', forecast_dss, shared_dir,
                         use_conic=False, alt_period=None, alt_period_string=None, balance_period_str='1Day')

    # invent flow-reg reservoir daily elevation record from folsom elev (used for timing only)
    invent_elevation('Natoma',forecast_dss,'//Folsom Lake/ELEV-FORECAST//1DAY/AMER_BC_SCRIPT/',123.0)
    invent_elevation('Natoma',forecast_dss,'//FOLSOM/STORAGE//1Month/AMER_BC_SCRIPT/',123.0)


def split_folsom_evap(forecast_dss,evap_rec):
    dssFmRec = HecDss.open(forecast_dss)
    tsc = dssFmRec.get(evap_rec,True)
    recparts = tsc.fullName.split('/')
    total_evap = tsc.values

    # 2/3 to NF branch
    recparts[2] = 'Folsom Evap NF Split'
    tsc.fullName = '/'.join(recparts)
    tsc.values = [0.667*total_evap[j] for j in range(tsc.numberValues)]
    dssFmRec.write(tsc)

    # 1/3 to SF branch, but using subtraction for numerical precision
    recparts[2] = 'Folsom Evap SF Split'
    tsc.fullName = '/'.join(recparts)
    sf_values =  [total_evap[j] - tsc.values[j] for j in range(tsc.numberValues)]
    tsc.values = sf_values
    dssFmRec.write(tsc)
    
    dssFmRec.close()

def split_nimbus_outflow(forecast_dss,nimbus_outflow_rec):
    #tsc_outflow = DSS_Tools.dss_read_ts_safe(forecast_dss,nimbus_outflow_rec)
    dssFm = HecDss.open(forecast_dss)
    tsc_outflow = dssFm.get(nimbus_outflow_rec,True)
    cfs2cms = 1.0/35.314666213 if tsc_outflow.units.lower() == 'cfs' else 1.0   # need CMS

    hatchery_constant = 1.7 # cms
    gate_max = 144.4 # cmm
    
    hatchery_flow = []
    spill_flow = []
    gated_spill_flow = []
      
    for vi, v in enumerate(tsc_outflow.values):
        
        gate_flow = tsc_outflow.values[vi]*cfs2cms - hatchery_constant
        hatchery_flow.append(hatchery_constant)
        if gate_flow > gate_max:
            spill_flow.append(gate_flow - gate_max)
            gated_spill_flow.append(gate_max)
        else:
            spill_flow.append(0.0)
            gated_spill_flow.append(gate_flow)

    tsc_outflow.units = 'cms'
    
    tsc_outflow.fullName = '/AMERICAN RIVER/LAKE NATOMA NIMBUS/FLOW-GATEDSPILLWAY//1Day/AMER_BC_SCRIPT/'
    tsc_outflow.values = gated_spill_flow
    dssFm.write(tsc_outflow)

    tsc_outflow.fullName = '/AMERICAN RIVER/LAKE NATOMA NIMBUS/FLOW-SPILLWAY//1Day/AMER_BC_SCRIPT/'
    tsc_outflow.values = spill_flow
    dssFm.write(tsc_outflow)

    tsc_outflow.fullName = '/AMERICAN RIVER/LAKE NATOMA NIMBUS/FLOW-FISHHATCHERY//1Day/AMER_BC_SCRIPT/'
    tsc_outflow.values = hatchery_flow
    dssFm.write(tsc_outflow)

    dssFm.close()


def write_qot_7outlets_flows(forecast_dss, starttime_str, endtime_str):
    '''Split Folsom forecast outflow: 
    1) subtract pumping/Sac Muni flow from total, since it seems to be included in total release in spreadsheet
    2) split remainder equally into 3 penstocks, but with max 32.5 cms/penstock + associated leakage, 
    3) any remaining flow over 3X32.5 + max leakage goes to spillway
    '''    
    #p3_max = 32.5 # cms - orignial code from Vanessa had this value
    p3_max = 49.0825 # back-calculated based on a total power outflow of 8000 cfs, per Vanessa 2025-03-06
    leakFraction = 0.35
    p_max_w_leakage = p3_max*3.0/(1.0-leakFraction)
    leakage_max = p_max_w_leakage*leakFraction
    
    dssFmRec = HecDss.open(forecast_dss)
    tsm_flow = dssFmRec.read('//FOLSOM/FLOW-RELEASE//1HOUR/AMER_BC_SCRIPT/')#, starttime_str, endtime_str)
    tsc_flow = DSS_Tools.standardize_interval(tsm_flow,'1day').getData()
    FaveFlow = tsc_flow.values
    tsc_muni_pump = dssFmRec.read('//FOLSOM/PUMPING (FP)//1Day/AMER_BC_SCRIPT/').getData()#, starttime_str, endtime_str).getData() # should end up w/ same timing as hourly above!
    muniPump = tsc_muni_pump.values
    
    cfs2cms = 1.0
    if tsc_flow.units.lower() == 'cfs':   # need CMS
        cfs2cms = 1.0/35.314666213

    # divide up total folson flow
    # 3 penstocks are evenly split, so we will link single DSS record to three outlets
    penstock = []
    leakage = []
    spill = []
    for j in range(tsc_flow.numberValues):
        dam_flow = (FaveFlow[j] - muniPump[j])*cfs2cms
        if dam_flow > p_max_w_leakage:
            leakage.append(leakage_max)
            penstock.append(p3_max)
            spill.append(dam_flow - p_max_w_leakage)
        else:
            leakage.append(dam_flow*leakFraction)
            penstock.append(dam_flow*(1.0-leakFraction)/3.0)
            spill.append(0.0)
        
    recparts = tsc_flow.fullName.split('/')
    recparts[5] = '1day'   # just to make sure path correct, maybe preventing DLL crash?
    tsc_flow.units = 'cms'
    
    recparts[3] = 'FLOW-RELEASE-SINGLEPENSTOCK'
    tsc_flow.fullName = '/'.join(recparts)
    tsc_flow.values = penstock
    dssFmRec.write(tsc_flow)

    recparts[3] = 'FLOW-RELEASE-SPILL'
    tsc_flow.fullName = '/'.join(recparts)
    tsc_flow.values = spill
    dssFmRec.write(tsc_flow)

    recparts[3] = 'FLOW-RELEASE-LEAKAGE'
    tsc_flow.fullName = '/'.join(recparts)
    tsc_flow.values = leakage
    dssFmRec.write(tsc_flow)

    dssFmRec.close()

    return True

def study_dir_from_run_dir(run_dir):
    w2sim,_ = os.path.split(run_dir)
    runs_dir,_ = os.path.split(w2sim)
    study_dir,_ = os.path.split(runs_dir)
    return study_dir

def model_dir_from_run_dir(run_dir,model_place,model_name):
    study_dir = study_dir_from_run_dir(run_dir)
    model_dir = os.path.join(study_dir,'cequal-w2',model_place,model_name)
    return model_dir

def remove_evap_from_inflows(forecast_dss, starttime_str, endtime_str):
    # these better all be inn the same units - should all be CFS from bc script...

    dssFmRec = HecDss.open(forecast_dss)
    tsc_nf = dssFmRec.read('//Folsom-NF-in/FLOW-IN//1Day/AMER_BC_SCRIPT/', starttime_str, endtime_str).getData()
    tsc_sf = dssFmRec.read('//Folsom-SF-in/FLOW-IN//1Day/AMER_BC_SCRIPT/', starttime_str, endtime_str).getData()
    tsc_evap = dssFmRec.read('/AMERICAN RIVER/FOLSOM LAKE/FLOW-ACC-DEP//1Day/AMER_BC_SCRIPT/', starttime_str, endtime_str).getData() # should be negative of evap

    nf_less_evap = []
    sf_less_evap = []

    for j in range(tsc_nf.numberValues):
        in_sum = tsc_nf.values[j]+tsc_sf.values[j]
        nf_less_evap.append( tsc_nf.values[j] + tsc_evap.values[j]*tsc_nf.values[j]/in_sum )
        sf_less_evap.append( tsc_sf.values[j] + tsc_evap.values[j]*tsc_nf.values[j]/in_sum )

    tsc_nf.fullName = '//Folsom-NF-in-minus-evap/FLOW-IN//1Day/AMER_BC_SCRIPT/'
    tsc_nf.values = nf_less_evap
    dssFmRec.write(tsc_nf)
    tsc_sf.fullName = '//Folsom-SF-in-minus-evap/FLOW-IN//1Day/AMER_BC_SCRIPT/'
    tsc_sf.values = sf_less_evap
    dssFmRec.write(tsc_sf)

    dssFmRec.close()

def subtract_muni_pump(forecast_dss):
    # not using DSS_Tools.add_or_subtract_flows b/c we are using full range of data, not run time window...
    DSS_Tools.resample_dss_ts(forecast_dss,'//FOLSOM/FLOW-RELEASE//1Hour/AMER_BC_SCRIPT/',None,forecast_dss,'1DAY')
    dssFmRec = HecDss.open(forecast_dss)
    tsc_outflow = dssFmRec.get('//FOLSOM/FLOW-RELEASE//1Day/AMER_BC_SCRIPT/',True)
    tsc_pump = dssFmRec.get('//FOLSOM/PUMPING (FP)//1Day/AMER_BC_SCRIPT/',True)
    out_to_natoma = []
    for j in range(tsc_outflow.numberValues):
        out_to_natoma.append( tsc_outflow.values[j] - tsc_pump.values[j])
    tsc_outflow.fullName = '//FOLSOM/FLOW-RELEASE-TO-NATOMA//1Day/AMER_BC_SCRIPT/'
    tsc_outflow.values = out_to_natoma
    dssFmRec.write(tsc_outflow)
    dssFmRec.close()


def load_tt_data(forecast_dss, starttime_str, endtime_str):

    dssFmRec = HecDss.open(forecast_dss)
    tsm_flow = dssFmRec.read('/AMERICAN RIVER/LAKE NATOMA/FLOW-NIMBUS ACTUAL//1Day/AMER_BC_SCRIPT/', starttime_str, endtime_str)
    tsm_at = dssFmRec.read('/MR Am.-Natoma Lake/Fair Oaks/Temp-Air//1Hour/251.40.53.1.1/', starttime_str, endtime_str)
    dssFmRec.close()

    #tsc_flow = DSS_Tools.standardize_interval(tsm_flow,'1day').getData()
    tsc_flow = tsm_flow.getData()
    tsc_at = DSS_Tools.standardize_interval(tsm_at,'1day').getData()

    FaveFlow = tsc_flow.values
    if tsc_flow.units.lower() == 'cfs':   # need CMS
        for j in range(tsc_flow.numberValues):
            FaveFlow[j] = FaveFlow[j] / 35.314666213

    Tair = tsc_at.values
    if tsc_at.units.lower() == 'f':   # need C
        for j in range(tsc_at.numberValues):
            Tair[j] = (Tair[j] - 32.0)* 5.0 / 9.0
    
    # get doy of year
    doys = []
    for j in range(tsc_flow.numberValues):
        doys.append(tsc_flow.getHecTime(j).dayOfYear())
    
    return doys,FaveFlow,Tair


def remove_folsom_lower_river_use(forecast_dss,lro_use_rec):
    dssFm = HecDss.open(forecast_dss)
    tsc = dssFm.get(lro_use_rec,True)
    values = []
    for i in range(tsc.numberValues):
        values.append(-1.0)
    tsc.values = values
    dssFm.write(tsc)
    dssFm.close()
    

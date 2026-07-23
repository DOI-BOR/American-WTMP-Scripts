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

# create list of unwanted folders in sys.path
search_list = ["SacTrn", "Sacramento", "American", "Stanislaus"]

# initialize list of paths to remove
matching_paths = []

# search sys.path for unwanted substrings
for p in sys.path:
    # check if any unwanted phrase appears in this path
    if any(phrase in p for phrase in search_list):
        matching_paths.append(p)

# print paths containing unwanted phrases
print("Paths to be removed:")
for path in matching_paths:
    print(path)

# remove unwanted paths from sys.path
for path in matching_paths:
    if path in sys.path:
        sys.path.remove(path)

# append project scripts directory to path
sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))

# Import and reload project modules
import DSS_Tools
reload(DSS_Tools)

import DMS_preprocess
reload(DMS_preprocess)

import equilibrium_temp
reload(equilibrium_temp)

import create_balance_flow_jython as cbfj
reload(cbfj)

# import DSS and timezone utilities
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

    #If x is a list, recursively interpolate each entry
    if isinstance(x, list):
        return [interp(point, xp, fp, left, right) for point in x]
    
    # default boundary behavior
    else:
        
        # Default left/right boundary values if not supplied
        if left is None:
            left = fp[0]
        
        if right is None:
            right = fp[-1]

        # left extrapolation
        if x < xp[0]:
            return left
        
        # right extrapolation        
        elif x > xp[-1]:
            return right
         
        # linear interpolation between points         
        else:
            
            # Search segment where xp[i] <= x <= xp[i+1]
            for i in range(len(xp) - 1):
                if x >= xp[i] and x <= xp[i+1]:
                    # Perform the linear interpolation
                    t = (x - xp[i]) / (xp[i+1] - xp[i])
                    # return interpolated value
                    return fp[i] + t * (fp[i+1] - fp[i])


def interp_monthly_coeff_daily(monthly_coeffs):
    """
    Converts a list of 12 monthly regression coefficients into a list of 366
    daily coefficient values covering an entire leap year.

    Monthly coefficients are anchored at approximate month midpoints (day-of-year).
    The sequence is padded at both ends so that interpolation wraps smoothly across
    the year boundary between December and January.

    Inputs:
      monthly_coeffs -- list of 12 floats; one coefficient per calendar month
                        (January = index 0, December = index 11)

    Output:
      Returns a list of 366 floats containing the linearly interpolated daily
      coefficient values for days 1 through 366 of a leap year.
    """
    
    # midpoints of each month, approximating 12 evenly spaced locations
    month_midpoints = [16.5, 46, 75.5, 106, 136.5, 167, 197.5, 228.5, 259, 289.5, 320, 350.5]    
    
    # pad endpoints to support interpolation wrap-around
    month_midpoints_padded = [-14.5] + month_midpoints + [366.0+15.5]
    
    # wrap coefficients to allow smooth interpolation across year boundary
    coeffs_padded = [monthly_coeffs[-1],] + monthly_coeffs + [monthly_coeffs[0],]
    
    # output daily coefficient array
    daily_coeff = []
    
    # interpolate for each day of year
    for i in range(1,367):
        daily_coeff.append(interp(i,month_midpoints_padded,coeffs_padded))
    return daily_coeff
    

def eq_temp(rtw,at,cl,ws,sr,td,eq_temp_out):
    """
    Computes the equilibrium water surface temperature time series from meteorological
    DSS inputs and writes the result (at hourly, daily, and weekly resolution) to a
    DSS output file.

    Reads air temperature, cloud cover, wind speed, shortwave radiation, and dewpoint
    temperature from DSS records over the run time window, then calls
    equilibrium_temp.calc_equilibrium_temp() to compute the equilibrium temperature
    at each time step. The hourly result is resampled to daily and weekly averages
    before writing.

    Inputs:
      rtw         -- WAT run time window object providing start/end time strings
      at          -- 2-element list [dss_file_path, dss_record]; air temperature (deg C or F)
      cl          -- 2-element list [dss_file_path, dss_record]; cloud cover fraction (0-1)
      ws          -- 2-element list [dss_file_path, dss_record]; wind speed (m/s)
      sr          -- 2-element list [dss_file_path, dss_record]; shortwave irradiance (W/m2)
      td          -- 2-element list [dss_file_path, dss_record]; dewpoint temperature (deg C)
      eq_temp_out -- 2-element list [output_dss_file_path, output_dss_record]; destination
                     for the equilibrium temperature output (deg C, INST-VAL)

    Output:
      No return value. Writes three DSS records to eq_temp_out[0]:
        - Hourly equilibrium temperature  (deg C, INST-VAL)
        - Daily average equilibrium temperature  (deg C, PER-AVER)
        - Weekly average equilibrium temperature (deg C, PER-AVER)
    """
    
    # read start and end time strings
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()

    # get at data and times in the formats needed
    dssFm = HecDss.open(at[0])        
    tsc = dssFm.read(at[1], starttime_str, endtime_str, False).getData()
    
    # store integer HEC times and datetime conversion
    tsc_int_times = tsc.times
    dtt = DSS_Tools.hectime_to_datetime(tsc)
    at_data = tsc.values
    dssFm.close()
    # get remaining meteorological variables using helper utility
    cl_data = DSS_Tools.data_from_dss(cl[0],cl[1],starttime_str,endtime_str)
    ws_data = DSS_Tools.data_from_dss(ws[0],ws[1],starttime_str,endtime_str)
    sr_data = DSS_Tools.data_from_dss(sr[0],sr[1],starttime_str,endtime_str)
    td_data = DSS_Tools.data_from_dss(td[0],td[1],starttime_str,endtime_str)

    # print first 48 rows for debugging (2 days of hourly data)
    for i in range(48):
        print(dtt[i].strftime('%Y-%m-%d %H:%M'),at_data[i],cl_data[i],sr_data[i],td_data[i],ws_data[i])
    
    # compute equilibrium temperature array using dedicated module
    Te = equilibrium_temp.calc_equilibrium_temp(dtt,at_data,cl_data,sr_data,td_data,ws_data)
    
    # unpack (unused tuple, likely legacy debug line)
    (dtt, at, cl, sr, td, ws)
    print('writing: ',eq_temp_out[1])
    
    # prepare DSS output container
    tsc = TimeSeriesContainer()
    tsc.times = tsc_int_times
    tsc.fullName = eq_temp_out[1]
    tsc.values = Te
    tsc.units = 'C'
    tsc.type = 'INST-VAL'
    tsc.numberValues = len(tsc.values)

    # convert to daily and weekly series
    tsm = tsmath(tsc)
    tsm_day = DSS_Tools.standardize_interval(tsm,'1day')
    tsm_wk = DSS_Tools.standardize_interval(tsm,'1week')
    
    # write all three series to DSS output file    
    dssFmOut = HecDss.open(eq_temp_out[0])
    dssFmOut.write(tsc)
    dssFmOut.write(tsm_day)
    dssFmOut.write(tsm_wk)
    dssFmOut.close()

def read_temp_schedule_csv(csv_file_path):
    """
    Reads a temperature schedule CSV file and returns its data rows as a list of lists.

    The expected file format is:
        Sch#,May,June,July,Aug,Sept,Oct,Nov
        1,63,63,63,63,63,56,56
        2,63,63,63,63,63,57,56
        ...

    The header row (row 0) is skipped. Each subsequent row represents one temperature
    schedule, with values corresponding to monthly target temperatures (deg F).

    Inputs:
      csv_file_path -- full path to the CSV temperature schedule file

    Output:
      Returns a list of lists of strings, one inner list per non-header CSV row.
      Each inner list contains the raw string fields from the row (including the
      schedule number in index 0 and monthly temperatures in indices 1 onward).
    """
    
     # initialize CSV reader
    schedReader = csv.reader(open(csv_file_path), delimiter=',', quotechar='|')
    
    monthlyTempsSched = []
    
     # loop over all rows in file
    for i,row in enumerate(schedReader):
        
        # skip header row
        if i>0: 
            monthlyTempsSched.append(row)
    
    return monthlyTempsSched


def get_interpolated_coeffs(location):
    """
    Returns four sets of daily regression coefficients (cf0, cf1, cf2, cf3)
    for the specified downstream temperature target location, interpolated
    from monthly values to daily resolution using interp_monthly_coeff_daily.

    Two locations are currently supported:
      - Location 1: Watt Avenue (American River below Nimbus Dam)
      - Location 3: Hazel Avenue (further downstream)

    The regression form for both locations is:
        T_downstream = cf0[d] + cf1[d]*T_air + cf2[d]*T_release + cf3[d]*log10(flow)

    where d is the zero-based day-of-year index.

    A third location (RiverMile) has placeholder coefficients defined but is
    marked with a TODO and is not returned by this function.

    Inputs:
      location -- integer; downstream target location identifier:
                    1 = Watt Avenue
                    3 = Hazel Avenue
                  (other values raise a ValueError)

    Output:
      Returns a tuple (cf0, cf1, cf2, cf3), each a list of 366 floats representing
      the daily interpolated regression coefficient for the selected location.
      Raises ValueError if the location identifier is not supported.
    """
    
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

    # select coefficient set based on location
    # Watt Ave
    if location==1: 
        coeffs = watt_coeffs
    
    # Hazel Ave
    elif location==3: 
        coeffs = hazel_coeffs
    
    else:
        raise ValueError("Forecast Preprocess: W2 downstream target not supported:"+str(location))
   
    # interpolate each monthly coefficient series to daily resolution
    cf0 = interp_monthly_coeff_daily(coeffs[0]) # int(tercpet)
    cf1 = interp_monthly_coeff_daily(coeffs[1]) # x1
    cf2 = interp_monthly_coeff_daily(coeffs[2]) # x2
    cf3 = interp_monthly_coeff_daily(coeffs[3]) # x3

    return cf0,cf1,cf2,cf3


def write_target_temp_npt(year,location,doys,Tair,FaveFlow,schedule_csv,targt_temp_npt_filepath,lagWatt=True):
    """
    Computes Folsom Dam release temperatures required to meet downstream temperature
    targets at the specified location, and writes the results to a .npt format file
    (plus an accompanying debug CSV file).

    For each day in May through November, the function applies a forward regression
    that inverts the downstream temperature target to a required Folsom release temperature:

        release_T = (-target_T + cf0[d] + cf1[d]*Tair + cf3[d]*log10(flow)) / (-cf2[d])

    An optional flow-dependent travel time lag (lagWatt) shifts the air temperature
    and flow inputs to account for routing time from Folsom to the target location.
    Release temperatures are capped at the target when the regression overshoots or
    when flow exceeds 300 CMS. Days outside May to November are filled with -99.0.

    One extra day is appended at the end to match the behavior of the original FORTRAN code.

    Inputs:
      year                    -- integer; the simulation year (used to build datetime objects)
      location                -- integer; downstream target location (1=Watt Ave, 3=Hazel Ave)
      doys                    -- list of floats; day-of-year values for each time step
      Tair                    -- list of floats; daily average air temperature (deg C)
      FaveFlow                -- list of floats; daily average flow (CMS)
      schedule_csv            -- full path to the temperature schedule CSV file
                                 (see read_temp_schedule_csv for format)
      targt_temp_npt_filepath -- full path to the output .npt file to be written
      lagWatt                 -- bool; if True, applies a flow-dependent lag when
                                 location == 1 (Watt Ave) (default True)

    Output:
      Returns True on successful completion.
      Writes two files:
        - targt_temp_npt_filepath         -- .npt format release temperature schedule
        - targt_temp_npt_filepath+'.debug.csv' -- debug CSV with inputs and coefficients
    """
    
    # load regression coefficients for selected location
    cf0,cf1,cf2,cf3 = get_interpolated_coeffs(location)

    # Read monthly temperature schedule CSV
    monthlyTempsSched = read_temp_schedule_csv(schedule_csv)

    # initialize output containers
    ReleaseTemp = []
    n99_line = [-99.0] * len(monthlyTempsSched)

    # Debug list with additional regression detai
    ReleaseTemp_debug = []

    #  Build base datetime from year + first DOY
    dt_base = dt.datetime(year,1,1) + dt.timedelta(days=doys[0]-1)
    one_day = dt.timedelta(days=1)
    
    # Number of days in time series
    n_days = len(doys)
    
    # Build datetime list for all days
    dtt = [dt_base + one_day*i for i in range(n_days)]

    # Loop over each da
    for i,date in enumerate(dtt):
        mon = date.month
        dayofyear = date.timetuple().tm_yday

        # Only compute targets for months May to Nov        
        if mon > 4 and mon < 12:          
            lag = 0
            
             # optional flow-dependent lag for Watt location
            if lagWatt and location==1:
                lag = int(3966.8*(35.314*FaveFlow[i])**(-0.944))
                lag = max(0,lag)
                
            # apply lagged date shift 
            dateLag = date + dt.timedelta(days=lag)

            # Convert lagged month into schedule index (May=1 to Nov=7)
            mlag = dateLag.month - 4  # index to schedule list above, accounting for lag date and jday as first index
            mlag = min(mlag,7)
            
            # Temperature index with lag, capped at end of record
            ilag = min(i + lag, n_days-1)  
            
             # day-of-year index for regression coefficients
            d = dayofyear - 1 

            # Start building row for this day's release temps
            rt = [dayofyear]
            
            # Loop over schedule rows
            for k,sched in enumerate(monthlyTempsSched):
                # Below is original code, which has some negatives that don't make a lot of sense? Copied to python exactly anyway
                # FORTRAN: ReleaseTemp(k,i)=(-((TTarg(month(k),i)-32)/1.8)+int(m)+x1(m)*aveTair(k)+x3(m)*log10(FaveFlow(k)))/(-x2(m))                
                # JYTHON:               (-1.0*((sched[mlag]-32.0)/1.8)+cf0[d] + cf1[d]*Tair[ilag] + cf3[d]*math.log10(FaveFlow[ilag]))/(-1.0*cf2[d])

                # Convert schedule target temp from F to C
                target_T = (float(sched[mlag]) - 32.0) / 1.8 

                # Forward regression formula (Folsom release temperature)
                release_T =  (-1.0*target_T + cf0[d] + cf1[d]*Tair[ilag] + cf3[d] * math.log10(FaveFlow[ilag])) / (-1.0*cf2[d])

                # TODO: insert RiverMile calc here, if location==3
                # Cap release temperature at target if flow is high or over-shoot occurs
                if release_T > target_T or FaveFlow[ilag] > 300.0:
                    release_T = target_T

                # store result for this schedule
                rt.append(release_T)
               
            # Append results for this day            
            ReleaseTemp.append(rt)    

            # save debug output (inputs + coefficients)            
            ReleaseTemp_debug.append(rt+[dayofyear,Tair[ilag],FaveFlow[ilag],cf0[d],cf1[d],cf2[d],cf3[d]])
        
        else:
            
            # Outside May to Nov: fill with -99 placeholders
            d = dayofyear - 1
            ReleaseTemp.append([dayofyear]+n99_line)
            ReleaseTemp_debug.append([dayofyear]+n99_line+[dayofyear,-99.0,-99.0,cf0[d],cf1[d],cf2[d],cf3[d]])

    # Add extra day at end (matches behavior of original FORTRAN code)
    ReleaseTemp.append([dayofyear+1]+n99_line)
    
    # Write to npt format
    with open(targt_temp_npt_filepath,'w') as fp:
        fp.write("   Temperature Target based release temps (May 1st - Nov 30) -- WTMP generated\n")
        
        # Write each row of data
        for rt_values in ReleaseTemp:
            
             # write DOY column
            fp.write("%8i,"%rt_values[0])
            
             # write release temperatures
            for rt in rt_values[1:]:
                fp.write("  %6.2f,"%rt)
                
            fp.write("\n")

    # Write corresponding debug CSV file
    with open(targt_temp_npt_filepath+'.debug.csv','w') as fp:
        
        fp.write("   Upstream Target Temp debug file -- WTMP generated\n")
        
        for rt_values in ReleaseTemp_debug:
            
            # write DOY
            fp.write("%8i,"%rt_values[0])
            
            # First part of data (all schedule temps)
            for rt in rt_values[1:-7]:
                fp.write("  %6.2f,"%rt)
                
            # Write day-of-year    
            fp.write("%8i,"%rt_values[-7])
            
            # Write regression variables and coefficients at full precision
            for rt in rt_values[-6:]:
                fp.write("  %.12e,"%rt)            
            
            fp.write("\n")

    return True


def calc_downstream_temp_W2(year,location,doys,Tair,Tmodel,FaveFlow,lagWatt=False):
    """
    Estimates downstream river temperature at the specified location using a reverse
    regression from Folsom release temperatures (from a CE-QUAL-W2 model result),
    air temperature, and flow.

    The reverse regression takes the form:
        T_downstream = cf0[d] + cf1[d]*Tair + cf2[d]*Tmodel[iBackLag] + cf3[d]*log10(flow)

    where d is the zero-based day-of-year index and iBackLag applies an optional
    backward lag to the model temperature to account for travel time.

    When flow exceeds 300 CMS, the regression is bypassed and the lagged model
    temperature is used directly as the downstream temperature.

    Inputs:
      year      -- integer; simulation year (used to build datetime objects)
      location  -- integer; downstream target location (1=Watt Ave, 3=Hazel Ave)
      doys      -- list of floats; day-of-year values for each time step
      Tair      -- list of floats; daily average air temperature (deg C)
      Tmodel    -- list of floats; Folsom daily release temperature from W2 model (deg C)
      FaveFlow  -- list of floats; daily average flow (CMS)
      lagWatt   -- bool; if True, applies a flow-dependent backward lag to Tmodel
                   when location == 1 (Watt Ave) (default False)

    Output:
      Returns a tuple (dtt, DownstreamTemp) where:
        dtt            -- list of Python datetime objects, one per day
        DownstreamTemp -- list of floats; estimated downstream temperature (deg C),
                          one per day
    """
    
    # Load interpolated daily regression coefficients
    cf0,cf1,cf2,cf3 = get_interpolated_coeffs(location)

    DownstreamTemp = []

    # Build daily datetime array
    dt_base = dt.datetime(year,1,1) + dt.timedelta(days=doys[0]-1)
    one_day = dt.timedelta(days=1)
    n_days = len(doys)
    dtt = [dt_base + one_day*i for i in range(n_days)]

    # Evaluate downstream temperature for each day
    for i,date in enumerate(dtt):
        mon = date.month
        dayofyear = date.timetuple().tm_yday 
        
        # Compute lag like upstream method
        lag = 0
        
        # optional Watt lag model
        if lagWatt and location==1:
            lag = int(3966.8*(35.314*FaveFlow[i])**(-0.944))
            lag = max(0,lag)
        
        dateLag = date + dt.timedelta(days=lag)
        
        # lagged indices
        ilag = min(i + lag, n_days-1)  # lag data index, but don't run past end of data
        d = dayofyear - 1 # day-of-year, minus 1, as index to coeffs        
        iBackLag = max(0,i-lag) # needed for reverse regression

        # below is original reverse regression code of Watt Ave:
        # avetemp1 is the daily model result temperature at Watt (Tmodel)
        # FORTRAN: WattTemp(k)=intW(m)+x1W(m)*aveTair(k)+x2W(m)*avetemp1(k-Lag)+x3W(m)*log10(FaveFlow(k))
        # JYTHON:              cf0[d] + cf1[d]*Tair[i] + cf2[d]*Tmodel[iBackLag] + cf3[d] * math.log10(FaveFlow[i])
        
        # Apply reverse regression for Watt/Hazel
        backTemp = cf0[d] + cf1[d]*Tair[i] + cf2[d]*Tmodel[iBackLag] + cf3[d] * math.log10(FaveFlow[i])

        # TODO: insert RiverMile calc here, if location==3

        # If high flows, bypass regression entirely
        if FaveFlow[i] > 300.0:
            backTemp = Tmodel[iBackLag]
            
        DownstreamTemp.append(backTemp)

    return dtt,DownstreamTemp


def update_W2_Folsom_iterative_schedule_number(run_dir,model_dir):
    """
    Updates the BestGuessTTS.csv file in the W2 Folsom model directory with the
    current ensemble number, which is read from a text file in the run directory.

    The BestGuessTTS.csv file controls which temperature schedule the W2 model
    uses during an iterative forecast run. This function overwrites the file with
    a fixed two-line format containing the current ensemble number.

    Inputs:
      run_dir   -- full path to the WAT run directory containing 'current_ensemble.txt'
      model_dir -- full path to the W2 Folsom model alternative directory containing
                   'BestGuessTTS.csv'

    Output:
      No return value. Overwrites BestGuessTTS.csv in model_dir with the format:
        Year,BG,CB
        1,<ensemble_number>,10
    """
    
    # Read current ensemble number from text file
    with open(os.path.join(run_dir,'current_ensemble.txt'),'r') as fp:
        ensemble_number = int(fp.readline().strip())

    # Update CSV file so line 2 reflects this ensemble number
    best_guess_file = os.path.join(model_dir,"BestGuessTTS.csv")
    
    with open(best_guess_file,'w') as fp:
        fp.write("Year,BG,CB\n1,%i,10\n"%ensemble_number)


def storage_to_elev(res_name,elev_stor_area,forecast_dss,storage_rec,conic=False):
    """
    Converts a storage time-series DSS record to an elevation time-series using
    the provided elevation-storage-area lookup table, and writes the elevation
    record back to the same DSS file.

    The C-part of the output DSS pathname is set to 'ELEV' and the B-part is set
    to res_name. Only linear interpolation is currently supported; conic interpolation
    will call sys.exit(-1).

    Inputs:
      res_name      -- string; reservoir name used as the B-part of the output DSS pathname
      elev_stor_area-- dict with keys 'elev', 'stor', 'area' (lists of floats); used for
                       storage-to-elevation interpolation
      forecast_dss  -- full path to the DSS file containing the storage record
      storage_rec   -- DSS pathname string of the storage record to convert
      conic         -- bool; if True, exits with an error (conic interpolation not yet
                       supported) (default False)

    Output:
      No return value. Writes the derived elevation time-series record (ft, INST-VAL)
      to forecast_dss with the same time axis as the storage record.
      Calls sys.exit(-1) if conic=True.
    """
    
    # open DSS file and read storage record
    dssFmRec = HecDss.open(forecast_dss)
    tsc = dssFmRec.get(storage_rec,True) # read ALL data in record

    elev = []
    
    # Conic interpolation not implemented; exit if requested
    if conic:
        print('Conic interpolation of elevations from storage not supported yet.')
        sys.exit(-1)
    
    else:
        
        # Loop through storage values and compute elevation by lookup
        for j in range(tsc.numberValues):
            elev.append(cbfj.linear_interpolation(elev_stor_area['stor'], elev_stor_area['elev'], tsc.values[j]))
            print('stor2elev: ',j,tsc.times[j],tsc.values[j],elev[-1])

    # Assign metadata for elevation record
    recparts = tsc.fullName.split('/')
    recparts[2] = res_name
    recparts[3] = 'ELEV'
    tsc.fullName = '/'.join(recparts)
    
    tsc.units = 'ft'
    tsc.type = 'INST-VAL'
    tsc.values = elev
    
    # write back to DSS
    dssFmRec.write(tsc)
    dssFmRec.close()

def invent_elevation(res_name,forecast_dss,storage_rec,elev_constant_ft):
    """
    Creates a constant elevation time-series DSS record by copying the time axis
    and metadata from an existing storage record and filling all values with a
    fixed elevation constant.

    Useful for reservoirs (e.g., Natoma) where elevation is approximately constant
    or where a placeholder elevation is needed for model configuration.

    The C-part of the output DSS pathname is set to 'ELEV' and the B-part is set
    to res_name.

    Inputs:
      res_name        -- string; reservoir name used as the B-part of the output DSS pathname
      forecast_dss    -- full path to the DSS file containing the source storage record
      storage_rec     -- DSS pathname string of the storage record whose time axis is copied
      elev_constant_ft-- float; constant elevation value (ft) to fill the output record

    Output:
      No return value. Writes the constant elevation time-series record (ft, INST-VAL)
      to forecast_dss with the same time axis as storage_rec.
    """
    
    # Read storage record structure to copy timing & metadata
    dssFmRec = HecDss.open(forecast_dss)
    tsc = dssFmRec.get(storage_rec,True)
    
    # Rewrite pathname for elevation
    recparts = tsc.fullName.split('/')
    recparts[2] = res_name
    recparts[3] = 'ELEV'
    
    # Fill time series with constant elevation
    tsc.fullName = '/'.join(recparts)
    tsc.units = 'ft'
    tsc.type = 'INST-VAL'
    tsc.values = [elev_constant_ft for j in range(tsc.numberValues)]
    
    dssFmRec.write(tsc)
    dssFmRec.close()

def write_forecast_elevations(currentAlternative, rtw, forecast_dss, shared_dir):
    """
    Computes and writes forecast reservoir elevation time series for Folsom Lake
    and Lake Natoma to the forecast DSS file.

    For Folsom Lake, monthly and daily storage records are converted to elevation
    using the Folsom elevation-storage-area lookup table. A forward mass-balance
    simulation (cbfj.predict_elevation) is then used to compute a daily elevation
    forecast driven by NF/SF inflows, a combined accumulation/depletion record,
    and resampled daily releases. The starting elevation is read from the monthly
    storage record for the last day of the preceding month.

    For Lake Natoma, a constant elevation (123.0 ft) is written using invent_elevation,
    since Natoma elevation is approximately constant.

    Inputs:
      currentAlternative -- WAT scripting alternative object for logging
      rtw                -- WAT run time window object providing start/end time strings
      forecast_dss       -- full path to the forecast DSS file for reading and writing
      shared_dir         -- full path to the shared directory containing the Folsom
                            elevation-storage-area CSV file ('AMR_scratch_Folsom.csv')

    Output:
      No return value. Writes the following DSS records to forecast_dss:
        - Folsom monthly elevation  (converted from storage)
        - Folsom daily elevation    (converted from storage)
        - Folsom daily release      (resampled from hourly to daily)
        - Folsom Lake ELEV-FORECAST (1-day forward mass-balance prediction)
        - Natoma ELEV               (constant 123.0 ft, monthly template)
        - Natoma ELEV               (constant 123.0 ft, daily template)
    """
    
    # Generate starting elevation timestamp (end of previous month)
    ht = HecTime(rtw.getStartTimeString())
    
    if ht.month() == 1:
        start_dt = dt.datetime(ht.year()-1,12,31)
        start_str = start_dt.strftime('%d%b%Y')+ ' 2400'
        
    else:
        start_dt = dt.datetime(ht.year(),ht.month(),1)
        start_dt = start_dt - dt.timedelta(days=1)
        start_str = start_dt.strftime('%d%b%Y')+ ' 2400'
        
    end_str = rtw.getEndTimeString()   

    # log computation window
    currentAlternative.addComputeMessage('Forecast Elevations: '+start_str+' '+end_str)

    # Load elevation-storage lookup table
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_Folsom.csv'), 'Folsom')

    # Convert monthly and daily storage to elevation
    storage_to_elev('Folsom',elev_stor_area,forecast_dss,'//FOLSOM/STORAGE//1Month/AMER_BC_SCRIPT/',conic=False)
    storage_to_elev('Folsom',elev_stor_area,forecast_dss,'/AMERICAN RIVER/FOLSOM LAKE/STORAGE-CVP//1Day/AMER_BC_SCRIPT/',conic=False)
    
    # Resample releases to daily to support elevation forecasting
    DSS_Tools.resample_dss_ts(forecast_dss,'//FOLSOM/FLOW-RELEASE//1Hour/AMER_BC_SCRIPT/',None,forecast_dss,'1DAY')
    
    # Define inflow/outflow record sets for elevation prediction
    inflow_records = ['//Folsom-NF-in/FLOW-IN//1Day/AMER_BC_SCRIPT/',
                      '//Folsom-SF-in/FLOW-IN//1Day/AMER_BC_SCRIPT/',
                      '/AMERICAN RIVER/FOLSOM LAKE/FLOW-ACC-DEP//1Day/AMER_BC_SCRIPT/']  # this actually evap, but negative already, so it goes as inflow
    outflow_records = ['//FOLSOM/FLOW-RELEASE//1Day/AMER_BC_SCRIPT/']
    
    # read starting elevation from DSS
    starting_elevation = DSS_Tools.first_value(forecast_dss,'//Folsom/ELEV//1Month/AMER_BC_SCRIPT/',start_str,end_str)
    print('starting_elevation ',starting_elevation)

    # Compute daily forecast elevation from inflows/outflows
    elev_calc_start_dt = start_dt + dt.timedelta(days=1)
    elev_calc_start_str = elev_calc_start_dt.strftime('%d%b%Y')+ ' 2400'
    
    # run reservoir routing / balance model
    cbfj.predict_elevation(currentAlternative, start_str, end_str, 'Folsom Lake', inflow_records, outflow_records, starting_elevation,
                         elev_stor_area, forecast_dss, '//Folsom Lake/ELEV-FORECAST//1DAY/AMER_BC_SCRIPT/', forecast_dss, shared_dir,
                         use_conic=False, alt_period=None, alt_period_string=None, balance_period_str='1Day')

    # Create Natoma elevation using constant value for timing
    invent_elevation('Natoma',forecast_dss,'//Folsom Lake/ELEV-FORECAST//1DAY/AMER_BC_SCRIPT/',123.0)
    invent_elevation('Natoma',forecast_dss,'//FOLSOM/STORAGE//1Month/AMER_BC_SCRIPT/',123.0)


def split_folsom_evap(forecast_dss,evap_rec):
    """
    Splits a Folsom Lake total evaporation time-series DSS record into two
    branch-specific records: a North Fork (NF) portion (2/3 of total) and a
    South Fork (SF) portion (1/3 of total), and writes both to the same DSS file.

    The split preserves the original time axis and metadata. The SF values are
    computed by subtracting the NF values from the total to avoid floating-point
    drift from independent scaling.

    Inputs:
      forecast_dss -- full path to the DSS file containing the evaporation record
      evap_rec     -- DSS pathname string of the total evaporation record to split

    Output:
      No return value. Writes two DSS records to forecast_dss:
        - 'Folsom Evap NF Split' (B-part): 66.7% of total evaporation
        - 'Folsom Evap SF Split' (B-part): 33.3% of total evaporation
          (remainder after NF subtraction)
    """
    
    # Open the DSS file containing evaporation time series
    dssFmRec = HecDss.open(forecast_dss)
    tsc = dssFmRec.get(evap_rec,True)
    
    # Parse record path and extract full evaporation time series values
    recparts = tsc.fullName.split('/')
    total_evap = tsc.values

    # 2/3 of evaporation assigned to NF branch
    recparts[2] = 'Folsom Evap NF Split'
    tsc.fullName = '/'.join(recparts)
    tsc.values = [0.667*total_evap[j] for j in range(tsc.numberValues)]
    dssFmRec.write(tsc)

    # Remaining 1/3 assigned to SF branch (using subtraction for accuracy)
    recparts[2] = 'Folsom Evap SF Split'
    tsc.fullName = '/'.join(recparts)
    sf_values =  [total_evap[j] - tsc.values[j] for j in range(tsc.numberValues)]
    tsc.values = sf_values
    dssFmRec.write(tsc)
    
    # Close DSS file after writing both split records
    dssFmRec.close()

def split_nimbus_outflow(forecast_dss,nimbus_outflow_rec):
    """
    Splits a Nimbus Dam total outflow time-series DSS record into three
    operational components - hatchery diversion, gated spillway, and uncontrolled
    spill - and writes each to the forecast DSS file.

    The allocation logic:
      - A constant hatchery demand (1.7 CMS) is subtracted from the total first.
      - The remaining gate flow is allocated to gated spillway up to a maximum
        capacity of 144.4 CMS; any excess goes to uncontrolled spill.
      - Output units are always CMS regardless of input units.

    Inputs:
      forecast_dss      -- full path to the DSS file containing the outflow record
      nimbus_outflow_rec-- DSS pathname string of the total Nimbus outflow record (CFS or CMS)

    Output:
      No return value. Writes three DSS records to forecast_dss (all in CMS):
        - /AMERICAN RIVER/LAKE NATOMA NIMBUS/FLOW-GATEDSPILLWAY//1Day/AMER_BC_SCRIPT/
        - /AMERICAN RIVER/LAKE NATOMA NIMBUS/FLOW-SPILLWAY//1Day/AMER_BC_SCRIPT/
        - /AMERICAN RIVER/LAKE NATOMA NIMBUS/FLOW-FISHHATCHERY//1Day/AMER_BC_SCRIPT/
    """
    
    # Open DSS file and read the Nimbus outflow time series
    dssFm = HecDss.open(forecast_dss)
    tsc_outflow = dssFm.get(nimbus_outflow_rec,True)
    
    # Convert units from CFS to CMS if needed
    cfs2cms = 1.0/35.314666213 if tsc_outflow.units.lower() == 'cfs' else 1.0   # need CMS

    # Define constants for hatchery flow and maximum gate capacity
    hatchery_constant = 1.7 # cms
    gate_max = 144.4 # cmm
    
    # Output category lists
    hatchery_flow = []
    spill_flow = []
    gated_spill_flow = []
      
    # Loop over each day's outflow value and allocate to components
    for vi, v in enumerate(tsc_outflow.values):
        
        # subtract hatchery demand from total flow
        gate_flow = tsc_outflow.values[vi]*cfs2cms - hatchery_constant
        
        # Hatchery flow is constant
        hatchery_flow.append(hatchery_constant)
        
         # split remaining flow between spill and gated spill
        if gate_flow > gate_max:
            spill_flow.append(gate_flow - gate_max)
            gated_spill_flow.append(gate_max)
        
        else:
            spill_flow.append(0.0)
            gated_spill_flow.append(gate_flow)

    # Set output units to CMS
    tsc_outflow.units = 'cms'
    
    # Write gated spillway output series
    tsc_outflow.fullName = '/AMERICAN RIVER/LAKE NATOMA NIMBUS/FLOW-GATEDSPILLWAY//1Day/AMER_BC_SCRIPT/'
    tsc_outflow.values = gated_spill_flow
    dssFm.write(tsc_outflow)

    # Write uncontrolled spill output series
    tsc_outflow.fullName = '/AMERICAN RIVER/LAKE NATOMA NIMBUS/FLOW-SPILLWAY//1Day/AMER_BC_SCRIPT/'
    tsc_outflow.values = spill_flow
    dssFm.write(tsc_outflow)

    # Write hatchery diversion output
    tsc_outflow.fullName = '/AMERICAN RIVER/LAKE NATOMA NIMBUS/FLOW-FISHHATCHERY//1Day/AMER_BC_SCRIPT/'
    tsc_outflow.values = hatchery_flow
    dssFm.write(tsc_outflow)

    # Close DSS file after writing all outputs
    dssFm.close()


def write_qot_7outlets_flows(forecast_dss, starttime_str, endtime_str):
    """
    Splits the daily average Folsom release flow into three CE-QUAL-W2 outlet
    components - single penstock flow, leakage, and spillway - and writes each
    to the forecast DSS file (all in CMS).

    The allocation logic:
      - Municipal pumping (FP) is subtracted from the total daily release first.
      - The remaining flow is divided between leakage (35% fraction) and three
        equal penstock flows (65% split three ways), each capped at 49.08 CMS.
      - Any flow exceeding the combined penstock + leakage maximum goes to spill.

    The hourly release record is resampled to daily before processing.

    Inputs:
      forecast_dss   -- full path to the DSS file containing release and pumping records
      starttime_str  -- start date/time string in HEC format (currently unused in reads,
                        retained for API consistency)
      endtime_str    -- end date/time string in HEC format (currently unused in reads,
                        retained for API consistency)

    Output:
      Returns True on successful completion.
      Writes three DSS records to forecast_dss (CMS, 1-day):
        - FLOW-RELEASE-SINGLEPENSTOCK  (flow through one of three equal penstocks)
        - FLOW-RELEASE-SPILL           (uncontrolled spill above turbine capacity)
        - FLOW-RELEASE-LEAKAGE         (leakage fraction of total gate flow)
    """  
    
    # Define turbine and leakage constraints
    p3_max = 49.0825 # back-calculated based on a total power outflow of 8000 cfs, per Vanessa 2025-03-06
    leakFraction = 0.35
    
    # Compute combined penstock + leakage threshold
    p_max_w_leakage = p3_max*3.0/(1.0-leakFraction)
    leakage_max = p_max_w_leakage*leakFraction
    
    # Open DSS file containing forecast data
    dssFmRec = HecDss.open(forecast_dss)
    
    # Read hourly release and convert to daily average
    tsm_flow = dssFmRec.read('//FOLSOM/FLOW-RELEASE//1HOUR/AMER_BC_SCRIPT/')#, starttime_str, endtime_str)
    tsc_flow = DSS_Tools.standardize_interval(tsm_flow,'1day').getData()
    FaveFlow = tsc_flow.values
    
    # Read municipal pumping flows
    tsc_muni_pump = dssFmRec.read('//FOLSOM/PUMPING (FP)//1Day/AMER_BC_SCRIPT/').getData()#, starttime_str, endtime_str).getData() # should end up w/ same timing as hourly above!
    muniPump = tsc_muni_pump.values
    
    # Convert to CMS if needed
    cfs2cms = 1.0
    if tsc_flow.units.lower() == 'cfs':   # need CMS
        cfs2cms = 1.0/35.314666213

    # divide up total folson flow
    # 3 penstocks are evenly split, so we will link single DSS record to three outlets
    
    # Initialize outlet component arrays
    penstock = []
    leakage = []
    spill = []
    
    # Split total flow into penstock, leakage, and spill
    for j in range(tsc_flow.numberValues):
        dam_flow = (FaveFlow[j] - muniPump[j])*cfs2cms
        
        # High-flow condition: exceed turbine capacity
        if dam_flow > p_max_w_leakage:
            
            # All penstocks at max, rest goes to spill
            leakage.append(leakage_max)
            penstock.append(p3_max)
            spill.append(dam_flow - p_max_w_leakage)
        else:
            
            # Distribute between leakage + equal 3-way penstock split
            leakage.append(dam_flow*leakFraction)
            penstock.append(dam_flow*(1.0-leakFraction)/3.0)
            spill.append(0.0)
    
    # Adjust pathname for output records    
    recparts = tsc_flow.fullName.split('/')
    recparts[5] = '1day'   # just to make sure path correct, maybe preventing DLL crash?
    tsc_flow.units = 'cms'
    
    # Write penstock record
    recparts[3] = 'FLOW-RELEASE-SINGLEPENSTOCK'
    tsc_flow.fullName = '/'.join(recparts)
    tsc_flow.values = penstock
    dssFmRec.write(tsc_flow)

    # Write spill record
    recparts[3] = 'FLOW-RELEASE-SPILL'
    tsc_flow.fullName = '/'.join(recparts)
    tsc_flow.values = spill
    dssFmRec.write(tsc_flow)

    # Write leakage record
    recparts[3] = 'FLOW-RELEASE-LEAKAGE'
    tsc_flow.fullName = '/'.join(recparts)
    tsc_flow.values = leakage
    dssFmRec.write(tsc_flow)

    # Close DSS file
    dssFmRec.close()

    return True

def study_dir_from_run_dir(run_dir):
    """
    Resolves the WAT study root directory from a given WAT run directory.

    The expected directory hierarchy is:
      study_dir / runs_dir / w2sim / run_dir

    Inputs:
      run_dir -- full path to the WAT run directory

    Output:
      Returns the full path to the study root directory (three levels up from run_dir).
    """
    
    # Strip last folder (simulation run)
    w2sim,_ = os.path.split(run_dir)
    
    # Go up another level (runs container)
    runs_dir,_ = os.path.split(w2sim)
    
    # One more level gives the parent study folder
    study_dir,_ = os.path.split(runs_dir)
    
    return study_dir

def model_dir_from_run_dir(run_dir,model_place,model_name):
    """
    Constructs the full path to a CE-QUAL-W2 model alternative directory from the
    WAT run directory and the model's place and name identifiers.

    The resulting path follows the convention:
      study_dir / cequal-w2 / model_place / model_name

    Inputs:
      run_dir     -- full path to the WAT run directory (used to resolve study_dir)
      model_place -- string; the subdirectory within cequal-w2 that organizes the model
                     (e.g., the reservoir or reach name)
      model_name  -- string; the specific CE-QUAL-W2 model alternative folder name

    Output:
      Returns the full path string to the model alternative directory.
    """
    
    # Identify the base study directory from the current run
    study_dir = study_dir_from_run_dir(run_dir)
    
    # Build path to specific CE-QUAL-W2 model directory
    model_dir = os.path.join(study_dir,'cequal-w2',model_place,model_name)
    
    return model_dir

def remove_evap_from_inflows(forecast_dss, starttime_str, endtime_str):
    """
    Partitions total Folsom Lake evaporation proportionally between the North Fork
    (NF) and South Fork (SF) inflow records, subtracting the evaporation share from
    each inflow and writing the adjusted records to the forecast DSS file.

    The evaporation (stored as a negative FLOW-ACC-DEP record) is allocated in
    proportion to the relative magnitude of NF and SF inflows at each time step.
    The resulting adjusted NF and SF inflows are written under new B-part names.

    Inputs:
      forecast_dss  -- full path to the forecast DSS file for reading and writing
      starttime_str -- start date/time string in HEC format for the read window
      endtime_str   -- end date/time string in HEC format for the read window

    Output:
      No return value. Writes two DSS records to forecast_dss:
        - //Folsom-NF-in-minus-evap/FLOW-IN//1Day/AMER_BC_SCRIPT/
        - //Folsom-SF-in-minus-evap/FLOW-IN//1Day/AMER_BC_SCRIPT/
    """
    
    # Open DSS file containing inflow and evaporation records
    dssFmRec = HecDss.open(forecast_dss)
    
    # Read NF and SF inflows
    tsc_nf = dssFmRec.read('//Folsom-NF-in/FLOW-IN//1Day/AMER_BC_SCRIPT/', starttime_str, endtime_str).getData()
    tsc_sf = dssFmRec.read('//Folsom-SF-in/FLOW-IN//1Day/AMER_BC_SCRIPT/', starttime_str, endtime_str).getData()
    
    # Read evaporation (negative values)
    tsc_evap = dssFmRec.read('/AMERICAN RIVER/FOLSOM LAKE/FLOW-ACC-DEP//1Day/AMER_BC_SCRIPT/', starttime_str, endtime_str).getData() # should be negative of evap

    nf_less_evap = []
    sf_less_evap = []

    # Partition evaporation proportional to NF/SF inflows
    for j in range(tsc_nf.numberValues):
        in_sum = tsc_nf.values[j]+tsc_sf.values[j]
        
        # Weighted subtraction of evaporation
        nf_less_evap.append( tsc_nf.values[j] + tsc_evap.values[j]*tsc_nf.values[j]/in_sum )
        sf_less_evap.append( tsc_sf.values[j] + tsc_evap.values[j]*tsc_nf.values[j]/in_sum )

    # Write modified NF record
    tsc_nf.fullName = '//Folsom-NF-in-minus-evap/FLOW-IN//1Day/AMER_BC_SCRIPT/'
    tsc_nf.values = nf_less_evap
    dssFmRec.write(tsc_nf)
    
    # Write modified SF record
    tsc_sf.fullName = '//Folsom-SF-in-minus-evap/FLOW-IN//1Day/AMER_BC_SCRIPT/'
    tsc_sf.values = sf_less_evap
    dssFmRec.write(tsc_sf)

    dssFmRec.close()

def subtract_muni_pump(forecast_dss):
    """
    Computes the net Folsom Lake release to Natoma by subtracting the daily
    municipal pumping (PUMPING FP) from the daily total release, and writes the
    result to the forecast DSS file.

    The hourly release record is first resampled to daily resolution before
    subtraction. The output record represents the flow actually delivered downstream
    to Lake Natoma.

    Inputs:
      forecast_dss -- full path to the forecast DSS file for reading and writing

    Output:
      No return value. Writes one DSS record to forecast_dss:
        - //FOLSOM/FLOW-RELEASE-TO-NATOMA//1Day/AMER_BC_SCRIPT/
    """
    
    # Resample hourly release to daily before subtraction
    DSS_Tools.resample_dss_ts(forecast_dss,'//FOLSOM/FLOW-RELEASE//1Hour/AMER_BC_SCRIPT/',None,forecast_dss,'1DAY')
    
    # Open DSS and read daily release + pumping
    dssFmRec = HecDss.open(forecast_dss)
    tsc_outflow = dssFmRec.get('//FOLSOM/FLOW-RELEASE//1Day/AMER_BC_SCRIPT/',True)
    tsc_pump = dssFmRec.get('//FOLSOM/PUMPING (FP)//1Day/AMER_BC_SCRIPT/',True)
    
    out_to_natoma = []
    
    # Subtract pumping from total flow to get actual downstream release
    for j in range(tsc_outflow.numberValues):
        out_to_natoma.append( tsc_outflow.values[j] - tsc_pump.values[j])
    
    # Update path for new flow record
    tsc_outflow.fullName = '//FOLSOM/FLOW-RELEASE-TO-NATOMA//1Day/AMER_BC_SCRIPT/'
    tsc_outflow.values = out_to_natoma
    
    # Write record
    dssFmRec.write(tsc_outflow)
    dssFmRec.close()


def load_tt_data(forecast_dss, starttime_str, endtime_str):
    """
    Reads the daily Nimbus flow and hourly Fair Oaks air temperature from the
    forecast DSS file, converts units to CMS and deg C respectively, and computes
    the day-of-year vector for use in downstream temperature regression calculations.

    The hourly air temperature is resampled to daily average before being returned.

    Inputs:
      forecast_dss  -- full path to the forecast DSS file for reading
      starttime_str -- start date/time string in HEC format for the read window
      endtime_str   -- end date/time string in HEC format for the read window

    Output:
      Returns a tuple (doys, FaveFlow, Tair) where:
        doys     -- list of floats; decimal day-of-year values for each daily time step
        FaveFlow -- list of floats; daily average Nimbus flow (CMS)
        Tair     -- list of floats; daily average Fair Oaks air temperature (deg C)
    """
    
    # Open DSS and read Nimbus actual flow record
    dssFmRec = HecDss.open(forecast_dss)
    tsm_flow = dssFmRec.read('/AMERICAN RIVER/LAKE NATOMA/FLOW-NIMBUS ACTUAL//1Day/AMER_BC_SCRIPT/', starttime_str, endtime_str)
    
    # Read Fair Oaks air temperature
    tsm_at = dssFmRec.read('/MR Am.-Natoma Lake/Fair Oaks/Temp-Air//1Hour/251.40.53.1.1/', starttime_str, endtime_str)
    dssFmRec.close()

    # Convert hourly air temp to daily average
    tsc_flow = tsm_flow.getData()
    tsc_at = DSS_Tools.standardize_interval(tsm_at,'1day').getData()

    FaveFlow = tsc_flow.values
    
    # Convert flow units from CFS to CMS if needed
    if tsc_flow.units.lower() == 'cfs':   # need CMS
        for j in range(tsc_flow.numberValues):
            FaveFlow[j] = FaveFlow[j] / 35.314666213

    # Convert air temp F to C if needed
    Tair = tsc_at.values
    
    # Convert air temperature from Fahrenheit to Celsius if needed
    if tsc_at.units.lower() == 'f':   # need C
        for j in range(tsc_at.numberValues):
            Tair[j] = (Tair[j] - 32.0)* 5.0 / 9.0
    
    # Compute day-of-year vector       
    doys = DSS_Tools.jday_from_tsc(tsc_flow)  
    
    return doys,FaveFlow,Tair


def get_downstream_loc(forecastDSS):
    """
    Reads the downstream temperature control location identifier from a DSS record
    and returns it as an integer.

    The location integer is stored as a text value in a special INTEGER-type DSS record
    and is used to select the appropriate regression coefficient set in
    get_interpolated_coeffs and related functions.

    Inputs:
      forecastDSS -- full path to the forecast DSS file containing the location record

    Output:
      Returns an integer identifying the downstream temperature control location:
        1 = Watt Avenue
        3 = Hazel Avenue
      Prints the raw record value and parsed integer for diagnostic purposes.
    """
    
    # Open DSS and read stored integer specifying downstream location
    dssFm = HecDss.open(forecastDSS)        
    tsc = dssFm.get('//DOWNSTREAM_CONTROL_LOC///INTEGER/AMER_TARGET_TEMP/', True) # this should be passed in a linked record at some point
    
    # Parse integer from text field
    loc = int(str(tsc.getText()).strip())
    print('Downstream Loc: ',str(tsc),loc)
    
    dssFm.close()
    return loc

def remove_folsom_lower_river_use(forecast_dss,lro_use_rec):
    """
    Disables the Folsom Lower River Outlet usage record by overwriting all values
    in the specified DSS record with -1.0, signaling to the downstream model that
    this outlet is not in service.

    Inputs:
      forecast_dss -- full path to the forecast DSS file for reading and writing
      lro_use_rec  -- DSS pathname string of the lower river outlet usage record to disable

    Output:
      No return value. Overwrites all values in lro_use_rec with -1.0 in forecast_dss.
    """
    
    # Overwrite record with placeholder values (disabling usage)
    dssFm = HecDss.open(forecast_dss)
    tsc = dssFm.get(lro_use_rec,True)
    
    # Replace all values with -1.0 (disabling)
    values = []
    
    for i in range(tsc.numberValues):
        values.append(-1.0)
        
    tsc.values = values
    dssFm.write(tsc)
    dssFm.close()

def modify_w2_selective_start_date(rtw,w2_sel_filepath):

    """
    Modifies the CE-QUAL-W2 selective withdrawal control file (w2_selective.npt)
    in-place to update the TEND and TSTR day-of-year fields for the current
    forecast start date.

    Lines 10 to 12 (1-indexed) contain TEND fields (index 6) that are updated to
    start_doy + 1. Lines 13 to 15 contain TSTR fields (index 5) that are updated to
    start_doy + 2. All other lines are copied unchanged. The start DOY is clamped
    to a minimum of 149 (May 29) so that the model always has valid initial conditions
    for the temperature management season.

    The expected file format is documented in the original docstring below.
    Fields are comma-delimited; lines are rewritten in-place using fileinput.

    Inputs:
      rtw              -- WAT run time window object providing the start time string
      w2_sel_filepath  -- full path to the w2_selective.npt file to be modified

    Output:
      Returns True on successful completion.
      Modifies w2_sel_filepath in-place.
      Prints a confirmation message with the modified DOY.
    """
    
    # Extract start date as DOY
    starttime_str = rtw.getStartTimeString()
    start_doy = DSS_Tools.hectime_to_julian(HecTime(starttime_str))
    
    # Clamp to minimum allowed DOY (149=May 29)
    start_doy = max(149,start_doy) # default is to set these lines to 150/151


    print('Attempting to modify w2_selective: %s'%w2_sel_filepath)
    
    # Process file in-place, modifying specific lines
    for line in fileinput.input(w2_sel_filepath, inplace=True):
        lineno = fileinput.filelineno() # 1 index (not python zero index)
        
        # Lines 10 to 12: Modify TEND field (index 6)
        if lineno >= 10 and lineno <= 12:
            tokens = line.rstrip().split(',') # need that rstrip()!
            tokens[6] = "%i"%(start_doy+1)
            print("%s" % ",".join(tokens))
        
        # Lines 13 to 15: Modify TSTR field (index 5)
        elif lineno >= 13 and lineno <=15:
            tokens = line.rstrip().split(',') # need that rstrip()!
            tokens[5] = "%i"%(start_doy+2)
            print("%s" % ",".join(tokens))
        
        # Other lines copied unchanged
        else:
            print("%s" % line.rstrip())

    print("Modifed w2_selective.ntp for Folsom auto gate start after start of doy: %i"%start_doy)
    return True


def update_W2_Folsom_iterative_restart_date_and_shutters(rtw,model_dir,w2_elevs):    

    """
    Updates two CE-QUAL-W2 model input files - w2_con.csv and folsom_in.npt -
    in-place to configure the restart date and initial shutter (gate) elevations
    for the current forecast run.

    In w2_con.csv:
      - Line 403: the restart DOY is written (May 1 / DOY 120 if before that date,
        otherwise start DOY + 1)
      - Lines 166, 167, 168: the initial elevation (column 0) for each of the three
        shutters is set from w2_elevs

    In folsom_in.npt:
      - Line 4: the restart DOY is written

    Inputs:
      rtw       -- WAT run time window object providing the start time string
      model_dir -- full path to the W2 Folsom model alternative directory containing
                   w2_con.csv and folsom_in.npt
      w2_elevs  -- list of three floats; the starting CE-QUAL-W2 gate elevation (m)
                   for each of the three shutters (see snap_inital_gate_elev_from_folsom_elev)

    Output:
      No return value. Modifies w2_con.csv and folsom_in.npt in-place in model_dir.
    """
        
    # Open the DSS file containing evaporation time series
    starttime_str = rtw.getStartTimeString()
    doy = DSS_Tools.hectime_to_julian(HecTime(starttime_str))
    
    # Split full pathname for later rewriting
    restart_doy = 120 if doy < 120 else doy + 1 # write restart file either May 1 or start day + 1

    # fileinput should backup and write over file, with inplace=True

    # replace doy in line 403 (index 402) in w2_con.csv
    w2_con = os.path.join(model_dir,"w2_con.csv")
    for line in fileinput.input(w2_con, inplace=True):
        lineno = fileinput.filelineno() # 1 index (not python zero index)
        if lineno==403: 
            print("%i,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,"%restart_doy)
        elif lineno==166:
            tokens = line.rstrip().split(',')
            tokens[0] = str(w2_elevs[0])
            print("%s" % ",".join(tokens))
        elif lineno==167:
            tokens = line.rstrip().split(',')
            tokens[0] = str(w2_elevs[1])
            print("%s" % ",".join(tokens))
        elif lineno==168:
            tokens = line.rstrip().split(',')
            tokens[0] = str(w2_elevs[2])
            print("%s" % ",".join(tokens))        
        else:
            print("%s" % line.rstrip())

    # replace doy in line 4 in folsom_in.npt
    folsom_in = os.path.join(model_dir,"folsom_in.npt")
    for line in fileinput.input(folsom_in, inplace=True):
        if fileinput.filelineno()==4:
            print("%i,,,,,,,,,,,,,,,,,,"%restart_doy)
        else:
            print("%s" % line.rstrip())


def snap_inital_gate_elev_from_folsom_elev(elev_in):
    """
    Reads the initial forecast shutter (gate) positions from DSS for each of the
    Folsom Dam selective withdrawal shutters, converts the ResSim elevations to the
    corresponding CE-QUAL-W2 gate elevations, and returns the values for writing
    into the W2 control file.

    For simulations starting at or before day 120 (approximately May 1), all shutters
    are set to the penstock default elevation (93.57 m) because the temperature
    management season has not yet started.

    For each shutter, the DSS pathname and file are resolved via
    DSS_Tools.getDataLocationDSSInfo, the elevation time series is read, and the
    second value (index 1) is used to avoid potential issues with the first record.
    The ResSim elevation is then mapped to a W2 gate elevation via
    snap_inital_gate_elev_from_folsom_elev.

    Inputs:
      rtw                -- WAT run time window object providing start/end time strings
      currentAlternative -- WAT scripting alternative object for resolving linked DSS paths
      computeOptions     -- WAT compute options object providing the compute DSS filename
      shutter_objs       -- list of WAT DataLocation objects; one per Folsom shutter,
                            ordered to match w2_con.csv rows 166, 167, 168

    Output:
      Returns a list of three floats; the W2 starting gate elevation (m) for each
      shutter, suitable for passing directly to
      update_W2_Folsom_iterative_restart_date_and_shutters.
      Returns [93.57, 93.57, 93.57] if the simulation start is on or before DOY 120.
    """
    
    # Highest gate configuration (all shutters installed / "all-in")
    if elev_in >= 395.:
        return 124.94 
        
    # One shutter removed ("one-out")    
    elif elev_in >= 355.:
        return 113.05 
    
    # Historical support for a deganged-middle shutter configuration.
    # This option is currently disabled because it is not intended to be
    # supported by this workflow.
    #elif elev_in >= 346.:
    #    return 349.0

    # Intermediate gate configuration.
    elif elev_in >= 320.:
        return 105.13  
        
    # Lower gate configuration corresponding to penstock operations.    
    elif elev_in >= 300.:
        return 93.57 
    
    # If elevation is below all defined thresholds, return the input value
    # unchanged as a fallback.
    else:
        return elev_in

def get_initial_shutter_positions(rtw, currentAlternative, computeOptions, shutter_objs):
    """
    Reads the initial forecast shutter (gate) positions from DSS for each of the
    Folsom Dam selective withdrawal shutters, converts the ResSim elevations to the
    corresponding CE-QUAL-W2 gate elevations, and returns the values for writing
    into the W2 control file.

    For simulations starting at or before day 120 (approximately May 1), all shutters
    are set to the penstock default elevation (93.57 m) because the temperature
    management season has not yet started.

    For each shutter, the DSS pathname and file are resolved via
    DSS_Tools.getDataLocationDSSInfo, the elevation time series is read, and the
    second value (index 1) is used to avoid potential issues with the first record.
    The ResSim elevation is then mapped to a W2 gate elevation via
    snap_inital_gate_elev_from_folsom_elev.

    Inputs:
      rtw                -- WAT run time window object providing start/end time strings
      currentAlternative -- WAT scripting alternative object for resolving linked DSS paths
      computeOptions     -- WAT compute options object providing the compute DSS filename
      shutter_objs       -- list of WAT DataLocation objects; one per Folsom shutter,
                            ordered to match w2_con.csv rows 166, 167, 168

    Output:
      Returns a list of three floats; the W2 starting gate elevation (m) for each
      shutter, suitable for passing directly to
      update_W2_Folsom_iterative_restart_date_and_shutters.
      Returns [93.57, 93.57, 93.57] if the simulation start is on or before DOY 120.
    """
    
    # Get simulation start and end times from the run-time window.
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()
    
    # Convert simulation start time to Julian day-of-year.
    doy = DSS_Tools.hectime_to_julian(HecTime(starttime_str))

    # W2 Folsom expects all shutters to start at penstock elevations before
    # the temperature management season begins (day 120).
    if doy <= 120:
        return [93.57,93.57,93.57]

    # Store converted W2 shutter elevations.
    w2_elevs = []
    
    # Process each shutter object independently.
    for shutter in shutter_objs:
        
        # Determine the DSS pathname and DSS file containing the shutter data.
        shutter_dss_rec,shutter_dss_filepath = DSS_Tools.getDataLocationDSSInfo(shutter, currentAlternative, computeOptions)
        
        # Retrieve the shutter elevation time series from DSS.
        shutter_elevs = DSS_Tools.data_from_dss(shutter_dss_filepath,shutter_dss_rec,starttime_str,endtime_str)
        
        # Use the first available forecast value (index 1) rather than the
        # initial record because the start record can occasionally be missing
        # or unreliable.
        #
        # Convert the ResSim shutter elevation to the corresponding W2 gate
        # elevation and store the result.
        w2_elevs.append(snap_inital_gate_elev_from_folsom_elev(shutter_elevs[1])) # use first position in case there are issues with start record, as is sometime the case
    
    # Return the list of W2 starting elevations for all shutters.
    return w2_elevs

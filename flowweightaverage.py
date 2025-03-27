from hec.heclib.dss import HecDss
import hec.hecmath.TimeSeriesMath as tsmath
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE
import math,sys,datetime
import DSS_Tools
reload(DSS_Tools)

import tz_offset
reload(tz_offset)

#version 2.1
#modified 11-30-2022 by Scott Burdick-Yahya

def organizeLocations(currentAlternative, locations):
    locations_list = []
    if len(locations) % 2 != 0:
        currentAlternative.addComputeMessage("Uneven amount of Flow/Temp pairings. Check inputs.")
        sys.exit(1)
    for li, location in enumerate(locations):
        tspath =str(currentAlternative.loadTimeSeries(location))
        tspath = DSS_Tools.fixInputLocationFpart(currentAlternative, tspath)
        if li % 2 == 0: 
            current_pair = [tspath]
        else:
            current_pair.append(tspath)
            locations_list.append(current_pair)
    return locations_list


def flow_in_cfs(units,flows):
    if units.lower()=='cfs':
        return flows
    elif units.lower()=='cms':
        values_converted = []
        for f in flows:
            values_converted.append(f * 35.314666213)
        return values_converted
    else:
        print('FWA2: flow units not known:',units)
        sys.exit(-1)

def temperature_in_C(units,temps):
    if units.lower()=='c' or units.lower()=='deg c':
        return temps
    elif units.lower()=='f' or units.lower()=='deg f':
        values_converted = []
        for t in temps:
            values_converted.append((t - 32.0)*5.0/9.0)
        return values_converted
    else:
        print('FWA2: temperature units not known:',units)
        sys.exit(-1)

def FWA2(currentAlt, dssFile, timewindow, DSSPaths_list, outputname, cfs_limit=None, bad_data_fill_tempC=None, 
         last_override=False,return_tsc=False):
    '''Made a new flow-weighted average temperature function; other one was producing weirdness '''
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)

    flow_total = []
    flowtemp_total = []
    n_pairs = []

    flow_limit = 0.0 if cfs_limit is None else cfs_limit
    fill_value = UNDEFINED_DOUBLE if bad_data_fill_tempC is None else bad_data_fill_tempC
    
    for dspi, dsspaths in enumerate(DSSPaths_list):
        flow_dss_path = dsspaths[0]
        temp_dss_path = dsspaths[1]
        currentAlt.addComputeMessage(str(flow_dss_path))
        print('FWA2 Reading:',flow_dss_path)
        tsc_flow = dssFm.read(flow_dss_path, starttime_str, endtime_str, False).getData()
        flows = flow_in_cfs(tsc_flow.units,tsc_flow.values)
        print('FWA2 Reading:',temp_dss_path)
        last_rec_valid = False  # test to see if we can use override values
        try:
            tsc_temp = dssFm.read(temp_dss_path, starttime_str, endtime_str, False).getData()
            temps = temperature_in_C(tsc_temp.units,tsc_temp.values)
            print('tscf',tsc_flow.values[0])
            print(flows[0])
            print('tsct',tsc_temp.values[0])
            print(temps[0])
    
            # use type of 1st temp record
            if dspi==0:
                nrecs = len(flows)
                temp_type = tsc_temp.type
    
            if len(flows) != nrecs or len(temps) != nrecs:
                currentAlt.addComputeMessage("FWA2: record lengths do not match!")
                print("FWA2: record lengths do not match!",nrecs,len(flows),len(temps))
                sys.exit(-1)
    
            for i in range(nrecs):
                if dspi==0:
                    n_pairs.append(0) # init counter for number of flow/temp pairs in weighted average
                    flow_total.append(0.0)
                    flowtemp_total.append(0.0)
                # perform a lot of checks on data
                #print(i,flows[i],temps[i])
                if not math.isnan(flows[i]) and not math.isnan(temps[i]):
                    if flows[i] > flow_limit and flows[i] < 9.0e6: # could lower upper limit to something relevant to watershed
                        if temps[i] >= 0.0 and temps[i] <= 80.0:
                            # passed the data checks
                            
                            n_pairs[i] += 1
                            flow_total[i] += flows[i]
                            flowtemp_total[i] += flows[i]*temps[i]
    
                            #print(dspi,i,n_pairs[i],flows[i],temps[i],flow_total[i],flowtemp_total[i])
            last_rec_valid = True
        except:
            currentAlt.addComputeMessage('FWA2: data not addeded for record: '+temp_dss_path)
            last_rec_valid = False
        
    fwat = []
    print('nrecs:',nrecs)
    for i in range(nrecs):
        if n_pairs[i] > 0:
            fwat.append(flowtemp_total[i]/flow_total[i])        
        else:
            fwat.append(fill_value)
        if last_override and last_rec_valid:
            if flows[i] > flow_limit and flows[i] < 9.0e6: # could lower upper limit to something relevant to watershed
                if temps[i] >= 0.0 and temps[i] <= 80.0:
                    fwat[i] = temps[i]
        #print(i,fwat[i])

    # use last temp container to write
    tsc_temp.type = temp_type
    tsc_temp.fullName = outputname
    tsc_temp.values = fwat
    dssFm.write(tsc_temp)
    dssFm.close()
    if return_tsc:
        return tsc_temp
    else:
        return 0

def FWA(currentAlt, dssFile, timewindow, DSSPaths_list, outputname, cfs_limit=None):
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)
    dss_data = {}
    for dspi, dsspaths in enumerate(DSSPaths_list):
        flow_dss_path = dsspaths[0]
        currentAlt.addComputeMessage(str(flow_dss_path))
        flowTS = dssFm.read(flow_dss_path, starttime_str, endtime_str, False)

#        readabledates = []
#        for i in range(len(flowTS.getData().times)):
#            hecdate = flowTS.getData().getHecTime(i).dateAndTime(4) #delete later
#            readabledates.append(hecdate)

        
        flowTS = flowTS.getData()
        hecstarttimes = flowTS.times
        
        flow_units = flowTS.units
        if flow_units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')
            flowvals = []
            for flow in flowTS.values:
                flowvals.append(flow * 35.314666213)
            dss_data[dspi] = {'flow': flowvals} #start the dict
        else:
            dss_data[dspi] = {'flow': flowTS.values} #start the dict

        
        if cfs_limit != None:
            for fi, flow in enumerate(dss_data[dspi]['flow']):
                if flow < cfs_limit:
#                    print('Flow of {0} removed for being under limit: {1}'.format(flow, cfs_limit)) #noisy..
                    dss_data[dspi]['flow'][fi] = 0
                

        temp_dss_path = dsspaths[1]
        currentAlt.addComputeMessage(str(temp_dss_path))
        TempTS = dssFm.read(temp_dss_path, starttime_str, endtime_str, False)
        TempTS = TempTS.getData()
        tempunits = TempTS.units
        temptype = TempTS.type
        dss_data[dspi]['temp'] = TempTS.values

#    print(readabledates)
        
    for dspi in dss_data.keys():
        flowtemps = []
        offset = 0
        
        for i, flow in enumerate(dss_data[dspi]['flow']):
            
            temp = dss_data[dspi]['temp'][i]
            flowtemps.append(flow*temp)
            
        dss_data[dspi]['flowtemp'] = flowtemps
        
    total_flows = []
    dspi = dss_data.keys()[0]
    for i, flow in enumerate(dss_data[dspi]['flow']):
        temptotalflow = flow
        for key in dss_data.keys():
            if key != dspi:
                temptotalflow += dss_data[key]['flow'][i]
        total_flows.append(temptotalflow)
        
    total_flowtemp = []
    dspi = dss_data.keys()[0]
    for i, flowtemp in enumerate(dss_data[dspi]['flowtemp']):
        temptotalflowtemp = flowtemp
        for key in dss_data.keys():
            if key != dspi:             
                if not math.isnan(dss_data[key]['flowtemp'][i]):
                    if math.isnan(temptotalflowtemp):
                        temptotalflowtemp = dss_data[key]['flowtemp'][i]
                    else:
                        temptotalflowtemp += dss_data[key]['flowtemp'][i]
        total_flowtemp.append(temptotalflowtemp)
    
    FW_Avg_vals = []
    for i, flow in enumerate(total_flows):
        flowtemp = total_flowtemp[i]
        if flow == 0:
            FW_Avg_vals.append(UNDEFINED_DOUBLE)
        else:
            FW_Avg_vals.append(flowtemp / flow)
    
    tsc = TimeSeriesContainer()
    tsc.times = hecstarttimes
    tsc.fullName = outputname
    tsc.values = FW_Avg_vals
    tsc.startTime = hecstarttimes[0]
    tsc.units = tempunits
    tsc.type = temptype
    tsc.endTime = hecstarttimes[-1]
    tsc.numberValues = len(FW_Avg_vals)
    tsc.startHecTime = timewindow.getStartTime()
    tsc.endHecTime = timewindow.getEndTime()
    dssFm.write(tsc)
    dssFm.close()
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(FW_Avg_vals)))
    return 0

def FWA_Daily(currentAlt, dssFile, timewindow, DSSPaths_list, outputname, cfs_limit=None,delay_days=0):
    starttime_str = timewindow.getStartTimeString()
    if delay_days > 0:
        dt_start = DSS_Tools.hec_str_time_to_dt(starttime_str) + datetime.timedelta(days=delay_days)
        starttime_str = dt_start.strftime('%d%b%Y %H%M')      
    
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)
    dss_data = {}
    for dspi, dsspaths in enumerate(DSSPaths_list):
        flow_dss_path = dsspaths[0]
        currentAlt.addComputeMessage(str(flow_dss_path))
        flowTS = dssFm.read(flow_dss_path, starttime_str, endtime_str, False)
        
        readabledates = []
        for i in range(len(flowTS.getData().times)):
            hecdate = flowTS.getData().getHecTime(i).dateAndTime(4) #delete later
            readabledates.append(hecdate)
        
        
        flowTS = flowTS.getData()
        hecstarttimes = flowTS.times
        daystart_idx = [ti for ti, t in enumerate(readabledates) if t.split(':')[0][-2:] == '24']
        daily_times = [hecstarttimes[t] for t in daystart_idx]
        daily_times_readable = [readabledates[t] for t in daystart_idx]
        #dayend_idx = daystart_idx + 24
        
        flow_units = flowTS.units
        if flow_units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')
            flowvals = []
            for flow in flowTS.values:
                flowvals.append(flow * 35.314666213)
            dss_data[dspi] = {'flow': flowvals} #start the dict
        else:
            dss_data[dspi] = {'flow': flowTS.values} #start the dict
        
        
        if cfs_limit != None:
            for fi, flow in enumerate(dss_data[dspi]['flow']):
                if flow < cfs_limit:
#                   print('Flow of {0} removed for being under limit: {1}'.format(flow, cfs_limit)) #noisy..
                    dss_data[dspi]['flow'][fi] = 0

        temp_dss_path = dsspaths[1]
        currentAlt.addComputeMessage(str(temp_dss_path))
        TempTS = dssFm.read(temp_dss_path, starttime_str, endtime_str, False)
        TempTS = TempTS.getData()
        tempunits = TempTS.units
        dss_data[dspi]['temp'] = TempTS.values
        
        flowtemps = []
        total_flows = []
        offset = 0
        
        for i, flow in enumerate(dss_data[dspi]['flow']): #thats everyday bro            
            temp = dss_data[dspi]['temp'][i]
            if i < 10:
                print(i,'Flow-temp pair:',flow,temp)
            flowtemps.append(flow*temp)
        
        FWTemps_daily = []
        flowsums_daily = []
        for ti in daystart_idx:
           dei = ti - 23
           if dei < 0:
                dei = 0
           daily_flowsums = sum(dss_data[dspi]['flow'][dei:ti+1])
           flowsums_daily.append(daily_flowsums)
           flowtemp_div_flows = []
           for flowtemp in flowtemps[dei:ti]:
                if daily_flowsums == 0:
                    flowtemp_div_flows.append(MISSING_DOUBLE)
                else:
                    flowtemp_div_flows.append(flowtemp/daily_flowsums)
           FWTemps_daily.append(sum(flowtemp_div_flows))
           
        dss_data[dspi]['FWTemps_daily'] = FWTemps_daily
        dss_data[dspi]['DailyFlowsums'] = flowsums_daily

    total_Daily_flowsums = []
    dspi = dss_data.keys()[0]
    dailyflowsums = dss_data[dspi]['DailyFlowsums']
    for i, dfs in enumerate(dailyflowsums):
        temptotalflow = dfs
        for key in dss_data.keys():
            if key != dspi:
                temptotalflow += dss_data[key]['DailyFlowsums'][i]
        total_Daily_flowsums.append(temptotalflow)

    total_FWA_FlowTemps = []
    dspi = dss_data.keys()[0]
    FWTemps = dss_data[dspi]['FWTemps_daily']
    for i, fwtd in enumerate(FWTemps):
        dflow_sum = dss_data[dspi]['DailyFlowsums'][i]
        FW_Flowtemps = fwtd * dflow_sum
        temp_Fwflowtemps_summed = FW_Flowtemps           
        for key in dss_data.keys():
            if key != dspi:
                dflow_sum = dss_data[key]['FWTemps_daily'][i]
                dFWTemp = dss_data[key]['DailyFlowsums'][i]
                temp_Fwflowtemps_summed += dflow_sum * dFWTemp
        total_FWA_FlowTemps.append(temp_Fwflowtemps_summed)
        
    FW_Avg_vals = []
    print('len total_FWA_FlowTemps:',len(total_FWA_FlowTemps))
    print('len total_Daily_flowsums:',len(total_Daily_flowsums))
    for i, tfwaft in enumerate(total_FWA_FlowTemps):
        total_Daily_flowsum = total_Daily_flowsums[i]
        if i < 10:
            print(i,'total_Daily_flowsum:',total_Daily_flowsum)
            print(i,'tfwaft:',tfwaft)
        if total_Daily_flowsum < 24.0:
            FW_Avg_vals.append(UNDEFINED_DOUBLE)
            #FW_Avg_vals.append(-5)
        else:
            FW_Avg_vals.append(tfwaft/total_Daily_flowsum)
    
    tsc = TimeSeriesContainer()
    tsc.times = daily_times
    tsc.fullName = outputname
    tsc.values = FW_Avg_vals
    tsc.startTime = daily_times[0]
    tsc.units = tempunits
    tsc.endTime = daily_times[-1]
    tsc.numberValues = len(FW_Avg_vals)
#    tsc.startHecTime = timewindow.getStartTime()
#    tsc.endHecTime = timewindow.getEndTime()
    dssFm.write(tsc)
    dssFm.close()
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(FW_Avg_vals)))
    return 0

def F_to_C(t,is_in_F):
    if is_in_F:
        return (t-32.)*5./9.
    else:
        return t

def FWA2_Daily(currentAlt, dssFile, timewindow, DSSPaths_list, outputname, cfs_limit=-1,delay_days=0,delay_hours=0):
    '''FWA_Daily is not working so I made this one'''
    
    starttime_str = timewindow.getStartTimeString()
    if delay_days > 0:
        dt_start = DSS_Tools.hec_str_time_to_dt(starttime_str) + datetime.timedelta(days=delay_days,hours=delay_hours)
        starttime_str = dt_start.strftime('%d%b%Y %H%M')      
    
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)

    for dspi, dsspaths in enumerate(DSSPaths_list):
        flow_dss_path = dsspaths[0]
        currentAlt.addComputeMessage(str(flow_dss_path))        
        tsc_flow = dssFm.read(flow_dss_path, starttime_str, endtime_str, False).getData()
        if dspi==0:
            # get datetimes for filtering by day
            dtt = DSS_Tools.hectime_to_datetime(tsc_flow)
        
        to_cfs = 1.0
        if tsc_flow.units.lower() == 'cms':
            to_cfs = 35.314666213

        temp_dss_path = dsspaths[1]
        currentAlt.addComputeMessage(str(temp_dss_path))
        tsc_temp = dssFm.read(temp_dss_path, starttime_str, endtime_str, False).getData()
        if dspi==0:    
            # make daily tsc container from input, arrays once we know length
            out_tsc = tsmath(tsc_temp).transformTimeSeries('1DAY',"","AVE").getData()
            nDaily = len(out_tsc.values)
            flow = [0.0 for i in range(nDaily)]
            flowtemp = [0.0 for i in range(nDaily)]
            fwa = [UNDEFINED_DOUBLE for i in range(nDaily)]
            
        to_C = False
        if tsc_temp.units.lower() == 'f' or tsc_temp.units.lower() == 'degf':
            to_C = True

        # iteraite over small-than-daily entries and put in sequenial day arrays
        di = 0
        dt0 = dtt[0] + tz_offset.timedelta
        last_day = dt0.day
        for i,(f,t) in enumerate(zip(tsc_flow.values,tsc_temp.values)):
            # increment day di when day changes
            dttz = dtt[i] + tz_offset.timedelta
            if i < 50:
                print(i,di,last_day,dttz.strftime('%d%b%Y %H%M'))
            if dttz.day != last_day:
                last_day = dttz.day
                di += 1
                if di >= nDaily:
                    break
            f1 = f*to_cfs
            if i < 50:
                print(i,di,last_day,dttz.strftime('%d%b%Y %H%M'),f1,F_to_C(t,to_C))
            if f1 > cfs_limit and t >= 0.0 and t < 120.:  # some basic filters, sometimes there is an undefined value in at start/end
                flow[di] += f1
                flowtemp[di] += f1*F_to_C(t,to_C)

    # find temp
    for i in range(nDaily):
        if i < 4:
            print(i,flowtemp[i],flow[i])
        if flow[i] > 0.0:
            fwa[i] = flowtemp[i]/flow[i]
    
    out_tsc.fullName = outputname
    out_tsc.units = 'C'
    out_tsc.values = fwa
    dssFm.write(out_tsc)
    dssFm.close()
        
    return 0
           




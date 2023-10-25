from hec.heclib.dss import HecDss
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE
import math

#version 2.1
#modified 02-16-2023 by Scott Burdick-Yahya

def organizeLocations(currentAlternative, locations):
    locations_list = []

    if len(locations) % 2 != 0:
        currentAlternative.addComputeMessage("Uneven amount of Flow/Temp pairings. Check inputs.")
        sys.exit(1)
    for li, location in enumerate(locations):
        tspath =str(currentAlternative.loadTimeSeries(location))
        tspath = fixInputLocationFpart(currentAlternative, tspath)
        if li % 2 == 0: 
            current_pair = [tspath]
        else:
            current_pair.append(tspath)
            locations_list.append(current_pair)
    return locations_list

def fixInputLocationFpart(currentAlternative, tspath):
    new_fpart_start = ':'.join(currentAlternative.getInputFPart().split(':')[:-1])
    tspath = tspath.split('/')
    fpart = tspath[6]
    fpart_split = fpart.split(':')
    new_fpart = new_fpart_start + ':' + fpart_split[-1]
    tspath[6] = new_fpart
    tspath = '/'.join(tspath)
    return tspath

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
            FW_Avg_vals.append(MISSING_DOUBLE)
        else:
            FW_Avg_vals.append(flowtemp / flow)
    
    tsc = TimeSeriesContainer()
    tsc.times = hecstarttimes
    tsc.fullName = outputname
    tsc.values = FW_Avg_vals
    tsc.startTime = hecstarttimes[0]
    tsc.units = tempunits
    tsc.endTime = hecstarttimes[-1]
    tsc.numberValues = len(FW_Avg_vals)
    tsc.startHecTime = timewindow.getStartTime()
    tsc.endHecTime = timewindow.getEndTime()
    dssFm.write(tsc)
    dssFm.close()
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(FW_Avg_vals)))
    return 0

def FWA_Daily(currentAlt, dssFile, timewindow, DSSPaths_list, outputname, cfs_limit=None):
    starttime_str = timewindow.getStartTimeString()
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
    for i, tfwaft in enumerate(total_FWA_FlowTemps):
        total_Daily_flowsum = total_Daily_flowsums[i]
        if total_Daily_flowsum == 0:
            FW_Avg_vals.append(MISSING_DOUBLE)
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
           

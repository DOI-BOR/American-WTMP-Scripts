'''
Created on 8/7/2023
@author: Scott Burdick-Yahya
@organization: Resource Management Associates
@contact: scott@rmanet.com
@note:
'''

import create_balance_flow_jython as cbfj
reload(cbfj)
from com.rma.model import Project
import os
import Simple_DSS_Functions as sdf
reload(sdf)

def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()

    # Folsom Inputs **********************************************************************
    # ******* Use same time resolution as ResSim hydro model time step ************
    # Flows are assumed to be period averaged
    # Evap assumed to be period accumulated length (e.g., ft)
    # Stage assumed to be instantaneous values
    
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    balance_period_str = currentAlternative.getTimeStep()
    shared_dir = os.path.join(project_dir, 'shared')

    DMS_hydro_dss_file = os.path.join(shared_dir, "DMS_AmericanHydroTS.dss")
    output_dss_file = os.path.join(shared_dir,'DMS_American_ResSim_Pre-Process.dss')
    fallback_dss_file = os.path.join(shared_dir,'WTMP_American_Historical.dss')

    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR Am.-Folsom Lake/NF American River-Flow/Flow//1Day/250.400.125.1.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR Am.-Folsom Lake/SF American River-Flow/Flow//1Day/250.402.125.1.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(output_dss_file,'/MR Am.-Folsom Lake/MormonR_NewcastlePP_Sum/Flow//1Day/ResSim_PreProcess/',rtw,output_dss_file,'1HOUR')

    inflow_records = ['::'.join([output_dss_file,'/MR Am.-Folsom Lake/NF American River-Flow/Flow//1Hour/250.400.125.1.1/']),
                      '::'.join([output_dss_file,'/MR Am.-Folsom Lake/SF American River-Flow/Flow//1Hour/250.402.125.1.1/']),
                      '::'.join([output_dss_file,'/MR Am.-Folsom Lake/MormonR_NewcastlePP_Sum/Flow//1Hour/ResSim_PreProcess/']),]

    outflow_records = ['/MR Am.-Folsom Lake/FOL-Generation Release U1/Flow//1Hour/250.3.125.4.1/',
                       '/MR Am.-Folsom Lake/FOL-Generation Release U2/Flow//1Hour/250.3.125.5.1/',
                       '/MR Am.-Folsom Lake/FOL-Generation Release U3/Flow//1Hour/250.3.125.6.1/',
                       '/MR Am.-Folsom Lake/FOL-Pumping Plant Release/Flow//1Hour/250.3.125.3.1/',
                       '/MR Am.-Folsom Lake/FOL-Spill Release/Flow//1Hour/250.3.125.7.1/',
                       '::'.join([output_dss_file,'/MR Am.-Folsom Lake/Upper_River_Outlets_Sum_min4/Flow//1Hour/ResSim_PreProcess/']),
                       '::'.join([output_dss_file,'/MR Am.-Folsom Lake/Upper_River_Outlets_Sum_min4/Flow//1Hour/ResSim_PreProcess/']),
                       '::'.join([os.path.join(shared_dir,'Folsom_balance_6.dss'),'//EID/FLOW/*/1HOUR/USGS-CARDNO-MERGED/']),]

    stage_record = '/MR Am.-Folsom Lake/FOL-Elevation/Elev//1Hour/250.3.145.1.1/'
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])

    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_Folsom.csv'), 'Folsom') #TODO: check this

    use_conic = True
    write_evap = False
    write_storage = False

    evap_dss_record_name = "/FOLSOM LAKE/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    storage_dss_record_name = "/FOLSOM LAKE/STORAGE/FLOW//1HOUR/DERIVED/"
    output_dss_record_name = "/FOLSOM LAKE/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    if use_conic:
        output_dss_record_name = "/FOLSOM LAKE/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
        if 'ZEROS' in evap_record:
            output_dss_record_name = "/FOLSOM LAKE/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"

    cbfj.create_balance_flows(currentAlternative, rtw, 'Folsom', inflow_records, outflow_records, stage_record, evap_record,
                                elev_stor_area, DMS_hydro_dss_file, output_dss_record_name, output_dss_file, shared_dir,
                                evap_dss_record_name=evap_dss_record_name, storage_dss_record_name=storage_dss_record_name,
                                balance_period_str=balance_period_str, use_conic=use_conic, write_evap=write_evap, write_storage=write_storage,
                                alt_period=1440, alt_period_string='1Day')


    # Natoma Inputs **********************************************************************
    # ******* Use same time resolution as ResSim hydro model time step ************
    # Flows are assumed to be period averaged
    # Evap assumed to be period accumulated length (e.g., ft)
    # Stage assumed to be instantaneous values

    inflow_records = ['/MR Am.-Folsom Lake/FOL-Generation Release U1/Flow//1Hour/250.3.125.4.1/',
                       '/MR Am.-Folsom Lake/FOL-Generation Release U2/Flow//1Hour/250.3.125.5.1/',
                       '/MR Am.-Folsom Lake/FOL-Generation Release U3/Flow//1Hour/250.3.125.6.1/',
                       '/MR Am.-Folsom Lake/FOL-Spill Release/Flow//1Hour/250.3.125.7.1/',
                       '::'.join([output_dss_file,'/MR Am.-Folsom Lake/Upper_River_Outlets_Sum_min4/Flow//1Hour/ResSim_PreProcess/']),
                       '::'.join([output_dss_file,'/MR Am.-Folsom Lake/Upper_River_Outlets_Sum_min4/Flow//1Hour/ResSim_PreProcess/']),]

    outflow_records = ['/MR Am.-Natoma Lake/NAT-Fish Hatchery Flow/Flow//1Hour/251.4.125.26.1/',
                       '::'.join([output_dss_file,'/MR Am.-Natoma Lake/NAT-Gen Release Sum/Flow//1Hour/ResSim_PreProcess/']),
                       '/MR Am.-Natoma Lake/NAT-South Canal Diversion Flw/Flow//1Hour/251.4.125.1.1/',
                       '/MR Am.-Natoma Lake/NAT-Spill Release/Flow//1Hour/251.4.125.6.1/']

    stage_record = '/MR Am.-Natoma Lake/NAT-Elevation/Elev//1Hour/251.4.145.1.1/'
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])

    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_Natoma.csv'), 'Natoma') #TODO: check this

    use_conic = True
    write_evap = False
    write_storage = False

    evap_dss_record_name = "/LAKE NATOMA/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    storage_dss_record_name = "/LAKE NATOMA/STORAGE/FLOW//1HOUR/DERIVED/"
    output_dss_record_name = "/LAKE NATOMA/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    if use_conic:
        output_dss_record_name = "/LAKE NATOMA/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
        if 'ZEROS' in evap_record:
            output_dss_record_name = "/LAKE NATOMA/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"

    cbfj.create_balance_flows(currentAlternative, rtw, 'Natoma', inflow_records, outflow_records, stage_record, evap_record,
                                elev_stor_area, DMS_hydro_dss_file, output_dss_record_name, output_dss_file, shared_dir,
                                evap_dss_record_name=evap_dss_record_name, storage_dss_record_name=storage_dss_record_name,
                                balance_period_str=balance_period_str, use_conic=use_conic, write_evap=write_evap, write_storage=write_storage,
                                alt_period=60*3, alt_period_string='3Hour')


    #######################################################################################
    # TODO: Calculate River balances
    #######################################################################################
    # Trinity River: Limekiln Gulch

    # Trinity River: Douglas City

    # Trinity River: Junction City

    # Clear Creek at South Fork junction (IGO)

    # Sacramento River at Bend Bridge

    return True

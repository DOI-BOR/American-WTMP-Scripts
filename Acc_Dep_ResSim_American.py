

import create_balance_flow_jython as cbfj
reload(cbfj)
from com.rma.model import Project
import os
import Simple_DSS_Functions as sdf
reload(sdf)

def computeAlternative(currentAlternative, computeOptions):
    # Log start of alternative computation.
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')
    
    # Retrieve DSS file and runtime window from compute options.
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
   
    # Determine balance computation timestep.
    balance_period_str = currentAlternative.getTimeStep()
    
    # Shared folder contains DSS files and lookup tables.
    shared_dir = os.path.join(project_dir, 'shared')

    # Define input hydro DSS file.
    DMS_hydro_dss_file = os.path.join(shared_dir, "DMS_AmericanHydroTS.dss")
    
    # Define preprocessing output DSS file.
    output_dss_file = os.path.join(shared_dir,'DMS_American_Pre-Process.dss')
    
    # Historical DSS file used as fallback source.
    fallback_dss_file = os.path.join(shared_dir,'WTMP_American_Historical.dss')

    # Resample Placerville daily flow to hourly resolution.
    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR Am.-Folsom Lake/11444500 Placerville-Daily Flow/Flow//1Day/250.201.125.1.1/',rtw,output_dss_file,'1HOUR')
    
    # Resample Mormon Ravine/Newcastle inflow data.
    sdf.resample_dss_ts(output_dss_file,'/MR Am.-Folsom Lake/MormonR_NewcastlePP_Sum/Flow//1Day/ResSim_PreProcess/',rtw,output_dss_file,'1HOUR')
    
    # Resample North Arm inflow data.
    sdf.resample_dss_ts(output_dss_file,'/MR Am.-Folsom Lake/North Arm/Flow//1Day/ResSim_PreProcess/',rtw,output_dss_file,'1HOUR')

    # Define inflow records contributing to Folsom Lake balance.
    inflow_records = ['::'.join([output_dss_file,'/MR Am.-Folsom Lake/North Arm/Flow//1Hour/ResSim_PreProcess/']),
                      '::'.join([output_dss_file,'/MR Am.-Folsom Lake/11444500 Placerville-Daily Flow/Flow//1Houray/250.201.125.1.1/']),
                      '::'.join([output_dss_file,'/MR Am.-Folsom Lake/MormonR_NewcastlePP_Sum/Flow//1Hour/ResSim_PreProcess/']),]

    # Define all releases and diversions leaving Folsom Lake.
    outflow_records = ['/MR Am.-Folsom Lake/FOL-Generation Release U1/Flow//1Hour/250.3.125.4.1/',
                       '/MR Am.-Folsom Lake/FOL-Generation Release U2/Flow//1Hour/250.3.125.5.1/',
                       '/MR Am.-Folsom Lake/FOL-Generation Release U3/Flow//1Hour/250.3.125.6.1/',
                       '/MR Am.-Folsom Lake/FOL-Pumping Plant Release/Flow//1Hour/250.3.125.3.1/',
                       '/MR Am.-Folsom Lake/FOL-Spill Release/Flow//1Hour/250.3.125.7.1/',
                       '::'.join([output_dss_file,'/MR Am.-Folsom Lake/Upper_River_Outlets_Sum_min4/Flow//1Hour/ResSim_PreProcess/']),
                       '::'.join([output_dss_file,'/MR Am.-Folsom Lake/Lower_River_Outlets_Sum_min4/Flow//1Hour/ResSim_PreProcess/']),
                       '::'.join([os.path.join(shared_dir,'Folsom_balance_6.dss'),'//EID/FLOW/*/1HOUR/USGS-CARDNO-MERGED/']),]

    # Reservoir elevation record used to estimate storage.
    stage_record = '/MR Am.-Folsom Lake/FOL-Elevation/Elev//1Hour/250.3.145.1.1/'
    
    # Currently using zero evaporation series.
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])

    # Read elevation-storage-area relationship for Folsom.
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_Folsom.csv'), 'Folsom') #TODO: check this

    # Enable conic interpolation for storage calculations.
    use_conic = True
    
    # Disable writing intermediate evaporation output.
    write_evap = False
    
    # Disable writing intermediate storage output.
    write_storage = False

    # Define DSS record names for optional intermediate outputs.
    evap_dss_record_name = "/FOLSOM LAKE/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    storage_dss_record_name = "/FOLSOM LAKE/STORAGE/FLOW//1HOUR/DERIVED/"
    
    # Default balance flow output record.
    output_dss_record_name = "/FOLSOM LAKE/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    
    # Modify output DSS pathname to indicate conic interpolation method.
    if use_conic:
        output_dss_record_name = "/FOLSOM LAKE/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
        
        # Further identify output as a no-evaporation balance calculation.
        if 'ZEROS' in evap_record:
            output_dss_record_name = "/FOLSOM LAKE/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"
    
    # Generate Folsom Lake balance flows using inflows, outflows, 
    # stage data, and elevation-storage relationships.
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

    # Natoma inflows are primarily releases from Folsom Reservoir.
    inflow_records = ['/MR Am.-Folsom Lake/FOL-Generation Release U1/Flow//1Hour/250.3.125.4.1/',
                       '/MR Am.-Folsom Lake/FOL-Generation Release U2/Flow//1Hour/250.3.125.5.1/',
                       '/MR Am.-Folsom Lake/FOL-Generation Release U3/Flow//1Hour/250.3.125.6.1/',
                       '/MR Am.-Folsom Lake/FOL-Spill Release/Flow//1Hour/250.3.125.7.1/',
                       '::'.join([output_dss_file,'/MR Am.-Folsom Lake/Upper_River_Outlets_Sum_min4/Flow//1Hour/ResSim_PreProcess/']),
                       '::'.join([output_dss_file,'/MR Am.-Folsom Lake/Lower_River_Outlets_Sum_min4/Flow//1Hour/ResSim_PreProcess/']),]

    # Define all known releases and diversions leaving Lake Natoma.
    outflow_records = ['/MR Am.-Natoma Lake/NAT-Fish Hatchery Flow/Flow//1Hour/251.4.125.26.1/',
                       '::'.join([output_dss_file,'/MR Am.-Natoma Lake/NAT-Gen Release Sum/Flow//1Hour/ResSim_PreProcess/']),
                       '/MR Am.-Natoma Lake/NAT-South Canal Diversion Flw/Flow//1Hour/251.4.125.1.1/',
                       '/MR Am.-Natoma Lake/NAT-Spill Release/Flow//1Hour/251.4.125.6.1/']

    # Lake Natoma elevation record used for storage calculations.
    stage_record = '/MR Am.-Natoma Lake/NAT-Elevation/Elev//1Hour/251.4.145.1.1/'
    
    # Use zero evaporation series for current balance calculations.
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])

    # Read Natoma elevation-storage-area lookup table.
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_Natoma.csv'), 'Natoma') #TODO: check this

    # Enable conic interpolation for storage estimation.
    use_conic = True
    
    # Disable writing evaporation time series.
    write_evap = False
    
    # Disable writing storage time series.
    write_storage = False

    # DSS pathname for optional evaporation output.
    evap_dss_record_name = "/LAKE NATOMA/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    
    # DSS pathname for optional storage output.
    storage_dss_record_name = "/LAKE NATOMA/STORAGE/FLOW//1HOUR/DERIVED/"
    
    # Default DSS pathname for balance flow output.
    output_dss_record_name = "/LAKE NATOMA/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    
    # Modify output pathname to reflect conic interpolation method.
    if use_conic:
        output_dss_record_name = "/LAKE NATOMA/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
        
        # Add no-evaporation designation when evaporation is disabled.
        if 'ZEROS' in evap_record:
            output_dss_record_name = "/LAKE NATOMA/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"

    # Generate Lake Natoma balance flows. 
    # Balance is computed using observed stage changes,
    # inflow records, outflow records, and storage relationships.
    cbfj.create_balance_flows(currentAlternative, rtw, 'Natoma', inflow_records, outflow_records, stage_record, evap_record,
                                elev_stor_area, DMS_hydro_dss_file, output_dss_record_name, output_dss_file, shared_dir,
                                evap_dss_record_name=evap_dss_record_name, storage_dss_record_name=storage_dss_record_name,
                                balance_period_str=balance_period_str, use_conic=use_conic, write_evap=write_evap, write_storage=write_storage,
                                alt_period=60*3, alt_period_string='3Hour')


    #######################################################################################
    # TODO: Calculate River balances
    #  
    # These locations are planned for implementation using the same 
    # mass-balance methodology applied to Folsom and Natoma. 
    # # Required inputs will likely include: 
    # - Upstream inflows 
    # - Tributary inflows 
    # - Diversions 
    # - Measured downstream flows 
    # - Optional storage or routing adjustments
     #######################################################################################
    
    # Trinity River: Limekiln Gulch

    # Trinity River: Douglas City

    # Trinity River: Junction City

    # Clear Creek at South Fork junction (IGO)

    # Sacramento River at Bend Bridge

    return True

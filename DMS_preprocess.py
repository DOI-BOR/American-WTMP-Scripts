
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

units_need_fixing = ['tenths','m/s','deg',] #'radians',]

def fix_DMS_types_units(dss_file):
    '''This method was implemented to change data types to PER-AVER that are not coming from the DMS that way'''
    dss = HecDss.open(dss_file)
    recs = dss.getPathnameList()
    for r in recs:
        tsm = dss.read(r)
        rlow = r.lower()
        if "/flow" in rlow or "/1day/" in rlow:
            tsm.setType('PER-AVER')
            dss.write(tsm)
        
        if tsm.getUnits() in units_need_fixing:
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

def DMS_fix_units_types(hydro_dss,met_dss_file):
    fix_DMS_types_units(hydro_dss)
    fix_DMS_types_units(met_dss_file)

def compute_folsom_flows(currentAlternative, rtw, hydro_dss, output_dss_file):
    # Add Mormon Ravine and Newcastle PP flows
    inflow_records = ['/MR Am.-Folsom Lake/11425416 Newcastle PP-Daily Flow/Flow//1Day/250.114.125.1.1/',
                      '/MR Am.-Folsom Lake/11433930 Mormon Ravine-Daily Flow/Flow//1Day/250.115.125.1.1/']
    DSS_Tools.add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Am.-Folsom Lake/MormonR_NewcastlePP_Sum/Flow//1Day/ResSim_PreProcess/', output_dss_file)

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

    output_dss_file = os.path.join(shared_dir,'DMS_American_ResSim_Pre-Process.dss') 

    hydro_dss = os.path.join(shared_dir, 'DMS_AmericanHydroTS.dss')
    fix_DMS_types_units(hydro_dss)
    met_dss_file = os.path.join(shared_dir,'DMS_AmericanMet.dss')
    fix_DMS_types_units(met_dss_file)

    splice_lewiston_met_data(currentAlternative, rtw, met_dss_file, output_dss_file, months=[1,2,3,12])
    compute_5Res_outflows(currentAlternative, rtw, hydro_dss, output_dss_file)


def preprocess_ResSim_American(currentAlternative, computeOptions):
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    balance_period = currentAlternative.getTimeStep()
    shared_dir = os.path.join(project_dir, 'shared')

    output_dss_file = os.path.join(shared_dir,'DMS_American_ResSim_Pre-Process.dss')

    hydro_dss = os.path.join(shared_dir, 'DMS_AmericanHydroTS.dss')
    fix_DMS_types_units(hydro_dss)
    met_dss_file = os.path.join(shared_dir,'DMS_AmericanMet.dss')
    fix_DMS_types_units(met_dss_file)

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


                        
    # if template IDs exist still, remove them
    #DSS_Tools.strip_templateID_and_rename_records(hydro_dss,currentAlternative)
    #DSS_Tools.strip_templateID_and_rename_records(met_dss_file,currentAlternative)

    compute_folsom_flows(currentAlternative, rtw, hydro_dss, output_dss_file)
    compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file)

    return True


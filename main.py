# ==============================================
# Automated System-wide Strength Evaluation Tool (ASSET) - Main Function
# Contributors: Pranav Sharma, Bin Wang, Leonardo Rese, Shahil Shah
# Last modified: 12/12/24
# Sharma, P., Rese, L., Wang, B., Vyakaranam, B., & Shah, S. (2023). Grid Strength Analysis for Integrating 30 GW of Offshore Wind Generation by 2030 in the U.S. Eastern Interconnection: Preprint. Paper presented at 22nd Wind and Solar Integration Workshop, Copenhagen, Denmark. https://www.nrel.gov/docs/fy24osti/87392.pdf

# ==============================================
# Copyright (c) 2025 Alliance for Sustainable Energy, LLC and the University of Texas at San Antonio
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ==============================================

import os, sys
import math
import time
from datetime import datetime

# =============================================================================================
# Environment Requirements:
# - Python 2.7 32-bit: https://www.python.org/download/releases/2.7/
# - IDE: PyCharm
# - Package: `xlrd` version 1.2.0 (Install: `pip install xlrd==1.2.0`)
# - PSS/E Software:
#   - PSS/E v34 Free Version (50 Buses):
#     https://psspy.org/psse-help-forum/question/6626/where-can-i-download-psse-version-34-for-students/
# =============================================================================================

# =============================================================================================
# Setting the PSS/E Installation Folder:
# Note: Update the paths below to match your PSS/E version.
# =============================================================================================
# If you know the installation folder, set it directly:
# Example:
# pssbindir = r"C:\Program Files (x86)\PTI\PSSE33\PSSBIN"
# pssepydir = r"""C:\Program Files\PTI\PSSE35\35.4\PSSPY39"""
# Uncomment and update for PSS/E v34:
# pssbindir = r"C:\Program Files (x86)\PTI\PSSE34\PSSBIN"
# pssepydir = r"""C:\Program Files (x86)\PTI\PSSE34\PSSPY27"""

# Alternatively, use the `pssepath` package for auto-detection:
# Install: `pip install pssepath`
import pssepath
pssepath.add_pssepath()

# =============================================================================================
# Load and Initialize the PSS/E API:
# Note: Ensure the correct version of PSS/E is loaded.
# =============================================================================================
# Uncomment the appropriate version:
# import pssexplore34
import psse34
import psspy, excelpy, dyntools, redirect, pssarrays

# Import constants for PSS/E
from psspy import _i, _f, _s

# Redirect PSS/E output to Python
redirect.psse2py()

# =============================================================================================

# Import other required packages
# Note: Uncomment xlrd if needed for Excel file handling
# import xlrd
import numpy as np
import assetlib
import pandas as pd

# Define the working directory
working_dir = os.getcwd()

# =============================================================================================
# Input Parameters (GUI under development):
# =============================================================================================
GUIinput_powerflowfile = working_dir + "\input\savnw.sav"  # Power flow file
GUIinput_poidata = working_dir + "\input\poidata.csv"  # POI (Point of Interconnection) data file
GUIinput_outputfolder = working_dir + "\output2"  # Folder for output files
GUIinput_DiscAllOwu = 1  # Include fault current contribution from IBR units? (0: No, 1: Yes)
GUIinput_SimMode = 2  # Simulation mode (0: Critical N-1/N-2, 1: Contingency scan, 2: SCRIF Computation)
GUIinput_K = 10  # Fault current analysis depth (branches K-level away from tested POI)
GUIinput_contfolder = working_dir + "\ctg"  # Folder for contingency data

# =============================================================================================
# Configuration Parameters:
# =============================================================================================
busNlimit = 50  # Maximum number of buses for initialization
ctg_num = 5  # Number of contingencies to consider


# =============================================================================================
# Main Function
# =============================================================================================

def main():
    """Main function to execute the ASSET tool."""

    # Create output folder
    out_dir = GUIinput_outputfolder
    try:
        os.mkdir(out_dir)
    except OSError:
        print("Output folder already exists.")

    # Create temporary folder for intermediate files
    tmp_dir = os.path.join(working_dir, "temp")
    try:
        os.mkdir(tmp_dir)
    except OSError:
        print("Temp folder already exists.")

    # Load power flow file and POI data
    PFcase = GUIinput_powerflowfile
    POIdatafile = GUIinput_poidata
    POIdata = pd.read_csv(POIdatafile)

    # Extract relevant data from POI file
    poi_list = POIdata['POI bus #'].values
    OSWid = POIdata['OW id'].tolist()
    poi_pmax = POIdata['MW capacity'].values

    # Initialize PSS/E
    psspy.psseinit(busNlimit)
    psspy.progress_output(6, "", [0, 0])
    psspy.alert_output(6, "", [0, 0])
    psspy.prompt_output(6, "", [0, 0])
    ierr = psspy.case(PFcase)
    assert ierr == 0, "SAV file cannot be opened."

    # =============================================================================================
    # Simulation Mode: Critical N-1/N-2
    if GUIinput_SimMode == 0:
        K = GUIinput_K
        Z9999_flag = GUIinput_DiscAllOwu

        SCMVA_POI, sum_brch, SCR_POI = [], [], []

        for i, POIi in enumerate(poi_list):
            print("Processing POI %d/%d" % (i + 1, len(poi_list)))

            # Reinitialize PSS/E for each POI
            psspy.psseinit(busNlimit)
            ierr = psspy.case(PFcase)
            assert ierr == 0, "SAV file cannot be opened."

            # Calculate short-circuit currents at K levels
            SCC_POI, SCC_kl, SCC_kl_array = assetlib.SC_k_lvl(psspy, POIi, OSWid[i], K, poi_list, Z9999_flag)

            # Identify top 2 critical branches
            SCC_ranking = np.flipud(np.argsort(SCC_kl_array[:, 2]))
            tp1, tp2 = assetlib.IdTop2Brch(PFcase, busNlimit, SCC_kl, SCC_kl_array, SCC_ranking)
            SCC_POIi_N1, SCC_POIi_N2 = assetlib.CalcSccN12(PFcase, busNlimit, POIi, OSWid[i], SCC_kl, SCC_kl_array, tp1, tp2, poi_list, Z9999_flag)

            # Append results
            SCMVA_POI.append([int(SCC_POI[0]), SCC_POI[1], SCC_POIi_N1, SCC_POIi_N2])
            sum_brch.append([int(SCC_POI[0]), int(SCC_kl_array[tp1[0]][0]), int(SCC_kl_array[tp1[0]][1]), SCC_kl[tp1[0]][2],
                             int(SCC_kl_array[tp2[0]][0]), int(SCC_kl_array[tp2[0]][1]), SCC_kl[tp2[0]][2]])
            SCR_POI.append([int(SCC_POI[0]), SCC_POI[1] / poi_pmax[i], SCC_POIi_N1 / poi_pmax[i], SCC_POIi_N2 / poi_pmax[i]])
       
       # Save results to CSV files
        pd.DataFrame(SCMVA_POI, columns=['Bus number of POI', 'SCMVA(N-0)', 'SCMVA(N-1)', 'SCMVA(N-2)']).to_csv(
            os.path.join(out_dir, "result_SCMVA.csv"), index=False)
        pd.DataFrame(sum_brch, columns=['Bus number of POI', 'From bus - Branch 1', 'To bus - Branch 1', 'ID - Branch 1',
                                        'From bus - Branch 2', 'To bus - Branch 2', 'ID - Branch 2']).to_csv(
            os.path.join(out_dir, "result_brch.csv"), index=False)
        pd.DataFrame(SCR_POI, columns=['Bus number of POI', 'SCR(N-0)', 'SCR(N-1)', 'SCR(N-2)']).to_csv(
            os.path.join(out_dir, "result_SCR.csv"), index=False)

    # =============================================================================================
    # Simulation Mode: Contingency Scan
    if GUIinput_SimMode == 1:
        now = datetime.now()
        res_dir = working_dir + '\\Results_' + now.strftime("%Y-%m-%d_%H-%M-%S")

        ctg_branch = pd.read_csv(GUIinput_contfolder + "\\CTGs_branch.csv")
        ctg_3wind = pd.read_csv(GUIinput_contfolder + "\\CTGs_3wind.csv")
        ctg_gentrip = pd.read_csv(GUIinput_contfolder + "\\CTGs_gentrip.csv")
        # ctg_swshunt = pd.read_csv(GUIinput_contfolder + "\\CTGs_swshunt.csv")
        ctg_discbus = pd.read_csv(GUIinput_contfolder + "\\CTGs_discbus.csv")
        ctg_clobranch = pd.read_csv(GUIinput_contfolder + "\\CTGs_clobranch.csv")
        # ctg_blckDC = pd.read_csv(GUIinput_contfolder + "\\CTGs_blckDC.csv")

        res_df = pd.DataFrame({"POI": poi_list})
        print('\n\n\nStarting SCR Computation...')
        start_time = time.time()

        for ctgi in range(0, ctg_num + 1, 1):
            print('Contingency: ' + str(ctgi))
            # start_time = time.time()
            psspy.psseinit(busNlimit)
            psspy.progress_output(6, "", [0, 0])
            psspy.alert_output(6, "", [0, 0])
            psspy.prompt_output(6, "", [0, 0])
            ierr = psspy.case(PFcase)
            assert ierr == 0, 'SAV file cannot be opened'

            ierr = psspy.short_circuit_coordinates(ival=1) # 0-rec, 1-polar
            ierr = psspy.short_circuit_units(ival=1) # 0-pu, 1-physical units
            # print("--- %s seconds ---" % (time.time() - start_time))

            ik = []
            scmva = []
            scr = []
            if ctgi == 0:
                ierr = psspy.bsys(sid=1, numbus=len(poi_list), buses=poi_list.tolist())
                results = pssarrays.iecs_currents(sid=1,
                                                  all=0,
                                                  flt3ph=1,
                                                  fltlg=0,
                                                  fltllg=0,
                                                  fltll=0,
                                                  fltloc=0,
                                                  linout=0,
                                                  linend=0,
                                                  tpunty=0,
                                                  lnchrg=1,
                                                  shntop=1,
                                                  dcload=0,
                                                  zcorec=0,
                                                  optnftrc=0,
                                                  loadop=1,
                                                  genxop=0,
                                                  brktime=0.08333)
                for k in range(0, len(poi_list)):
                    ik.append(abs(results.flt3ph[k].ia1))
                    ierr, nomv = psspy.busdat(ibus=poi_list[k], string='BASE')
                    scmva.append(abs(results.flt3ph[k].ia1) * nomv * math.sqrt(3) / 1000)
                    scr.append(abs(results.flt3ph[k].ia1) * nomv * math.sqrt(3) / 1000 / poi_pmax[k])
                res_df['SCMVA_' + str(ctgi)] = scmva
                res_df['SCR_' + str(ctgi)] = scr
            else:
                # Trip: AC lines and two-winding transformers
                tmp = ctg_branch[ctg_branch['CTG num'] == ctgi]
                tmp = tmp.reset_index(drop=True)
                if len(tmp) > 0:
                    for idx_tmp in range(len(tmp)):
                        from_bus = tmp.loc[idx_tmp, 'From bus']
                        to_bus = tmp.loc[idx_tmp, 'To bus']
                        br_id = tmp.loc[idx_tmp, 'ID']
                        ierr, v1 = psspy.busdat(ibus=from_bus, string='BASE')
                        ierr, v2 = psspy.busdat(ibus=to_bus, string='BASE')
                        if v1 == v2:
                            ierr = psspy.branch_chng_3(ibus=from_bus,
                                                       jbus=to_bus,
                                                       ckt=str(br_id),
                                                       intgar=[0, _i, _i, _i, _i, _i],
                                                       realar=[_f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f],
                                                       ratings=[_f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f],
                                                       namear=_s)
                        else:
                            ierr, realaro = psspy.two_winding_chng_5(ibus=from_bus,
                                                                     jbus=to_bus,
                                                                     ckt=str(br_id),
                                                                     intgar=[0, _i, _i, _i, _i, _i, _i, _i, _i, _i, _i,
                                                                             _i, _i, _i, _i],
                                                                     realari=[_f, _f, _f, _f, _f, _f, _f, _f, _f, _f,
                                                                              _f, _f, _f, _f, _f, _f, _f, _f, _f, _f,
                                                                              _f],
                                                                     ratings=[_f, _f, _f, _f, _f, _f, _f, _f, _f, _f,
                                                                              _f, _f],
                                                                     namear=_s,
                                                                     vgrpar=_s)
                del tmp

                # Close: AC lines and two-winding transformers
                tmp = ctg_clobranch[ctg_clobranch['CTG num'] == ctgi]
                tmp = tmp.reset_index(drop=True)
                if len(tmp) > 0:
                    for idx_tmp in range(len(tmp)):
                        from_bus = tmp.loc[idx_tmp, 'From bus']
                        to_bus = tmp.loc[idx_tmp, 'To bus']
                        br_id = tmp.loc[idx_tmp, 'ID']
                        ierr, v1 = psspy.busdat(ibus=from_bus, string='BASE')
                        ierr, v2 = psspy.busdat(ibus=to_bus, string='BASE')
                        if v1 == v2:
                            ierr = psspy.branch_chng_3(ibus=from_bus,
                                                       jbus=to_bus,
                                                       ckt=str(br_id),
                                                       intgar=[1, _i, _i, _i, _i, _i],
                                                       realar=[_f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f],
                                                       ratings=[_f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f],
                                                       namear=_s)
                        else:
                            ierr, realaro = psspy.two_winding_chng_5(ibus=from_bus,
                                                                     jbus=to_bus,
                                                                     ckt=str(br_id),
                                                                     intgar=[1, _i, _i, _i, _i, _i, _i, _i, _i, _i, _i,
                                                                             _i, _i, _i, _i],
                                                                     realari=[_f, _f, _f, _f, _f, _f, _f, _f, _f, _f,
                                                                              _f, _f, _f, _f, _f, _f, _f, _f, _f, _f,
                                                                              _f],
                                                                     ratings=[_f, _f, _f, _f, _f, _f, _f, _f, _f, _f,
                                                                              _f, _f],
                                                                     namear=_s,
                                                                     vgrpar=_s)
                del tmp

                # Trip: three-winding transformers
                tmp = ctg_3wind[ctg_3wind['CTG num'] == ctgi]
                tmp = tmp.reset_index(drop=True)
                if len(tmp) > 0:
                    for idx_tmp in range(len(tmp)):
                        bus1 = tmp.loc[idx_tmp, 'From bus']
                        bus2 = tmp.loc[idx_tmp, 'To bus 1']
                        bus3 = tmp.loc[idx_tmp, 'To bus 2']
                        br_id = tmp.loc[idx_tmp, 'ID']
                        ierr, realaro = psspy.three_wnd_imped_chng_4(ibus=bus1,
                                                                     jbus=bus2,
                                                                     kbus=bus3,
                                                                     ckt=str(br_id),
                                                                     intgar=[_i, _i, _i, _i, _i, _i, _i, 0, _i, _i, _i,
                                                                             _i, _i],
                                                                     realari=[_f, _f, _f, _f, _f, _f, _f, _f, _f, _f,
                                                                              _f, _f, _f, _f, _f, _f, _f],
                                                                     namear=_s,
                                                                     vgrpar=_s)
                del tmp

                # Disconnect bus
                tmp = ctg_discbus[ctg_discbus['CTG num'] == ctgi]
                tmp = tmp.reset_index(drop=True)
                if len(tmp) > 0:
                    for idx_tmp in range(len(tmp)):
                        bus = tmp.loc[idx_tmp, 'Bus']
                        ierr = psspy.bus_chng_4(ibus=bus,
                                                inode=_i,
                                                intgar=[4, _i, _i, _i],
                                                realar=[_f, _f, _f, _f, _f, _f, _f],
                                                name=_s)
                del tmp

                # Disconnect generator
                tmp = ctg_gentrip[ctg_gentrip['CTG num'] == ctgi]
                tmp = tmp.reset_index(drop=True)
                if len(tmp) > 0:
                    for idx_tmp in range(len(tmp)):
                        bus = tmp.loc[idx_tmp, 'Bus']
                        geid = tmp.loc[idx_tmp, 'ID']
                        ierr = psspy.machine_chng_2(ibus=bus,
                                                    id=str(geid),
                                                    intgar=[0, _i, _i, _i, _i, _i],
                                                    realar=[_f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f,
                                                            _f, _f],
                                                    )
                                                    # name=_s)
                del tmp

                treeobj = psspy.treedat(sizeislands=_i)
                a = treeobj.get('island_busnum')
                islanded_buses = [element for innerList in a for element in innerList]
                ierr = psspy.island()

                poi_list_new = poi_list.copy()
                for idx2, bus_tmp in enumerate(poi_list):
                    if (bus_tmp in islanded_buses):
                        poi_list_new.remove(bus_tmp)

                ierr = psspy.bsys(sid=1, numbus=len(poi_list_new), buses=poi_list_new.tolist())
                results = pssarrays.iecs_currents(sid=1,
                                                  all=0,
                                                  flt3ph=1,
                                                  fltlg=0,
                                                  fltllg=0,
                                                  fltll=0,
                                                  fltloc=0,
                                                  linout=0,
                                                  linend=0,
                                                  tpunty=0,
                                                  lnchrg=1,
                                                  shntop=1,
                                                  dcload=0,
                                                  zcorec=0,
                                                  optnftrc=0,
                                                  loadop=1,
                                                  genxop=0,
                                                  brktime=0.08333)
                flag = False
                for n1, poi_old in enumerate(poi_list):
                    for n2, poi_new in enumerate(poi_list_new):
                        if poi_old == poi_new:
                            flag = True
                            idxicc = n2
                    if flag == True:
                        ik.append(abs(results.flt3ph[idxicc].ia1))
                        ierr, nomv = psspy.busdat(ibus=poi_old, string='BASE')
                        scmva.append(abs(results.flt3ph[idxicc].ia1) * nomv * math.sqrt(3) / 1000)
                        scr.append(abs(results.flt3ph[idxicc].ia1) * nomv * math.sqrt(3) / 1000 / poi_pmax[idxicc])
                        flag = False
                    else:
                        ik.append(0.0)
                        scmva.append(0.0)
                        scr.append(0.0)
                res_df = pd.concat((res_df, pd.DataFrame({'SCMVA_' + str(ctgi): scmva})), axis=1)
                res_df = pd.concat((res_df, pd.DataFrame({'SCR_' + str(ctgi): scr})), axis=1)
                # res_df['SCMVA_'+ str(ctgi)] = scmva
                # res_df['SCR_'+ str(ctgi)] = scr

        res_df.to_csv(GUIinput_outputfolder + '\\results_CTG_SCAN_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv', index=False)
        print("Contingency Scan Completed.")
    # =============================================================================================
    # Simulation Mode: SCRIF Computation
    elif GUIinput_SimMode == 2:    
        print("SCRIF mode selected")
        K = GUIinput_K
        Z9999_flag = GUIinput_DiscAllOwu
        SCMVA_POI = []
        V_base = []
        # Step 1: Get base case SCMVA_POI and V_base for all POIs
        for i in range(len(poi_list)):
            POIi = poi_list[i]
            print("Computing base case for POI %d/%d" % (i + 1, len(poi_list)))
            psspy.psseinit(busNlimit)
            ierr = psspy.case(PFcase)
            if ierr != 0:
                raise Exception("SAV file cannot be opened.")

            # Get base SCC
            SCC_POI, SCC_kl, SCC_kl_array = assetlib.SC_k_lvl(psspy, POIi, OSWid[i], K, poi_list, Z9999_flag)
            SCMVA_POI.append([int(SCC_POI[0]), SCC_POI[1]])

            # Get base voltage
            ierr, v = psspy.busdat(POIi, 'PU')
            if ierr != 0:
                raise Exception("Failed to get voltage for bus %d" % POIi)
            V_base.append(abs(v))
        # Step 2: Loop over POI buses to compute SCRIF
        SCRIF_POI = []
        V_base = np.array(V_base)
        for i in range(len(poi_list)):
            POIi = poi_list[i]
            print("Computing SCRIF for POI %d/%d" % (i + 1, len(poi_list)))
            psspy.psseinit(busNlimit)
            ierr = psspy.case(PFcase)
            if ierr != 0:
                raise Exception("SAV file cannot be opened.")
            # Add negative load
            psspy.load_data_5(POIi, "TP", [_i, _i, _i, _i, _i, _i, _i], [_f, -100.0, _f, _f, _f, _f, _f, _f])
            psspy.fnsl([0, 0, 0, 0, 0, 1, 0, 0])  # Power flow
            # Get voltage at all POI buses
            V_del = []
            for POIj in poi_list:
                ierr, v = psspy.busdat(POIj, 'PU')
                if ierr != 0:
                    raise Exception("Failed to get voltage for bus %d" % POIj)
                V_del.append(abs(v))
            V_del = np.array(V_del)
            # Compute Impact Factor
            delta_i = abs(V_del[i] - V_base[i])
            Imp_Factor = []
            for j in range(len(poi_list)):
                if i == j:
                    Imp_Factor.append(1.0)
                else:
                    delta_j = abs(V_del[j] - V_base[j])
                    IF = 0.0
                    if delta_i > 1e-7:
                        IF = delta_j / delta_i
                        if IF < 1e-7:
                            IF = 0.0
                    Imp_Factor.append(IF)
            # Purge load and rerun PF
            psspy.purgload(POIi, "TP")
            psspy.fnsl([0, 0, 0, 0, 0, 1, 0, 0])
            # Compute SCRIF
            denom = 0.0
            for j in range(len(poi_list)):
                denom += Imp_Factor[j] * poi_pmax[j]
            scrif = SCMVA_POI[i][1] / denom if denom > 0 else 0
            SCRIF_POI.append([SCMVA_POI[i][0], scrif])
        # Save results
        pd.DataFrame(SCMVA_POI, columns=['Bus number of POI', 'SCMVA(N-0)']).to_csv(
            os.path.join(out_dir, "result_SCMVA_SCRIF.csv"), index=False)
        pd.DataFrame(SCRIF_POI, columns=['Bus number of POI', 'SCRIF']).to_csv(
            os.path.join(out_dir, "result_SCRIF.csv"), index=False) 
    print("Simulation complete.")
# Main function execution
if __name__ == "__main__":
    main() 
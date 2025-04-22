# ==============================================
# Automated System-wide Strength Evaluation Tool (ASSET) - Main Function
# Contributors: Pranav Sharma, Bin Wang, Leonardo Rese, Shahil Shah
# Last modified: 12/12/24
# Sharma, P., Rese, L., Wang, B., Vyakaranam, B., & Shah, S. (2023). Grid Strength Analysis for Integrating 30 GW of Offshore Wind Generation by 2030 in the U.S. Eastern Interconnection: Preprint. Paper presented at 22nd Wind and Solar Integration Workshop, Copenhagen, Denmark. https://www.nrel.gov/docs/fy24osti/87392.pdf

# ==============================================
# Copyright (c) 2024 Alliance for Sustainable Energy, LLC and the University of Texas at San Antonio
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

import sys, os
import numpy as np


# =============================================================================================
# Set PSS/E installation folder
# Remark: change to match the required PSS/E version
# =============================================================================================
# If the users know the exact installation folder, they can directly set it up:
# pssbindir = r"C:\Program Files (x86)\PTI\PSSE33\PSSBIN"
# pssepydir=(r"""C:\Program Files\PTI\PSSE35\35.4\PSSPY39""")

# If the users do not know the installation folder, they can use pssepath package for auto-setup
# package pssepath: https://pypi.org/project/pssepath/
import pssepath
pssepath.add_pssepath()

# =============================================================================================
# Load and initialize PSS/E API
# Remark: change to match the required PSS/E version
# =============================================================================================
import psse34
import psspy,excelpy,dyntools,redirect,pssarrays

from psspy import _i
from psspy import _f
from psspy import _s
redirect.psse2py()
# =============================================================================================


working_dir = os.getcwd()




def SC_k_lvl(psspy, busN, OSW_id, k, all_POIs, Z9999_flag):
	tmp_dir = working_dir+"\\temp"
	SCC_kl = []

	# Run SC analysis
	# set OSW Z=9999 for one or all OW plants
	if Z9999_flag == 1:
		psspy.machine_chng_2(busN, OSW_id, [psspy._i, psspy._i, psspy._i,
											psspy._i, psspy._i, psspy._i],
							 [psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, 9999.0,
							  psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f])
	else:
		for busi in all_POIs:
			psspy.machine_chng_2(busi, OSW_id, [psspy._i, psspy._i, psspy._i,
												psspy._i, psspy._i, psspy._i],
								 [psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
								  9999.0, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
								  psspy._f])

	psspy.short_circuit_coordinates(1)
	ierr = psspy.short_circuit_units(ival=1)
	psspy.progress_output(6,"",[0,0])
	psspy.bsys(1, 0, [0.0, 0.0], 0, [], 2, [busN], 0, [], 0, [])
	# psspy.iecs_4(1, 0, [1, 0, 0, 0, 3, 3, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0], [0.08333, 1.1], "", "", "")

	temp = sys.stdout	# store original stdout object for later

	sys.stdout = open(tmp_dir+"\\POI_" + str(busN) + "_lvl_" + str(k) + ".txt","w")
	sid = 1
	all = 0

	flt3ph = 1
	fltlg = 0
	fltllg = 0
	fltll = 0
	rptop = 3

	# k
	fltloc = 0
	linout = 0
	linend = 0
	tpunty = 0

	lnchrg = 1
	shntop = 1
	dcload = 0
	zcorec = 0
	cfactor = 0 # optnftrc in iecs_current()

	loadop = 1 # orig 0
	genxop = 0


	brktime = 0.08333
	psspy.iecs_4(sid, all, [flt3ph, fltlg, fltllg, fltll, rptop, k, fltloc, linout, linend, tpunty, lnchrg, shntop, dcload, zcorec, cfactor, loadop, genxop], [brktime, 1.1], "", "", "")
	sys.stdout.close()

	sys.stdout = temp	# restore print commands to interactive prompt

	f = open(tmp_dir+"\\POI_" + str(busN) + "_lvl_" + str(k) + ".txt",'r')
	# f = open('POI_126287_lvl_3.txt', 'r')
	alllines = f.readlines()
	size = len(alllines)
	i = 0
	strbus = ' X------------ BUS ------------X '
	strbus_len = len(strbus)
	stratbus = ' AT BUS'
	stratbus_len = len(stratbus)
	strfrom = ' X----------- FROM ------------X '
	strfrom_len = len(strfrom)
	strend = ' --------------------------------------'
	strend_len = len(strend)
	str3PH = '3PH'
	sc_k = []
	while i <= size-1:
		aa = alllines[i]
		size_aa = len(aa)
		if size_aa > strbus_len and aa[0:strbus_len] == strbus:
			bb = alllines[i+1]
			cc = bb.split()
			ind = cc.index(str3PH)
			c_mag = cc[ind+1]
			c_ang = cc[ind+2]
			bus_num = cc[0]
			bus_ind = int(bus_num)
			sc_mag = float(c_mag)
			if c_ang.find('-') == -1:
				ind_minus = None
			else:
				ind_minus = c_ang.index('-')
			if ind_minus != None:
				if ind_minus>1:
					c_ang = c_ang[0:ind_minus-1]
			sc_ang = float(c_ang)
			sc_k_0 = np.array([bus_ind,sc_mag])
		if size_aa > stratbus_len and aa[0:stratbus_len] == stratbus:  # go to AT BUS
			bb = alllines[i]
			cc = bb.split()
			ind = 2
			bus_to = cc[ind]
			bus_to_ind = int(bus_to)
		if size_aa > strfrom_len and aa[0:strfrom_len] == strfrom:
			checkedallfrombus = 0
			i = i + 1
			while checkedallfrombus == 0:
				bb = alllines[i]
				if len(bb)>strend_len and bb[0:strend_len] == strend:
					checkedallfrombus = 1
					continue
				else:
					cc = bb.split()
					bus_from = cc[0]
					flag = bus_from.isdigit()
					if len(cc)>1 and flag:
						ind = bb.find(']')
						cc = bb[ind + 1:].split()
						line_id = cc[1]
						ele_cc = cc[3]
						if ele_cc.find('-') == -1:
							ind_minus = None
						else:
							ind_minus = ele_cc.index('-')
						if ind_minus != None:
							if ind_minus > 1:
								ele_cc = ele_cc[0:ind_minus - 1]
						ele_sc = float(ele_cc)
						if int(bus_from)<=int(bus_to): # make sure the same branch is only included once
							SCC_kl.append([int(bus_from), int(bus_to), line_id, ele_sc])
				i = i + 1
		i = i + 1

	leng = len(SCC_kl)
	SCC_kl_array = np.zeros([leng, 3])
	for i in range(leng):
		SCC_kl_array[i][0] = float(SCC_kl[i][0])
		SCC_kl_array[i][1] = float(SCC_kl[i][1])
		SCC_kl_array[i][2] = float(SCC_kl[i][3])
		pass
	return sc_k_0, SCC_kl, SCC_kl_array


def IdTop2Brch(PFcase, busNlimit, SCC_kl, SCC_kl_array, SCC_ranking):
	tp1 = []
	tp2 = []
	for i in range(len(SCC_ranking)):
		psspy.psseinit(busNlimit)
		psspy.case(PFcase)

		Line_RX = np.asarray(psspy.abrncplx(-1, 1, 1, 1, 1, ['RX'])[1][0])
		Line_chg = np.asarray(psspy.abrnreal(-1, 1, 1, 1, 1, ['CHARGING'])[1][0])
		Line_from = np.asarray(psspy.abrnint(-1, 1, 1, 1, 1, ['FROMNUMBER'])[1][0])
		Line_to = np.asarray(psspy.abrnint(-1, 1, 1, 1, 1, ['TONUMBER'])[1][0])
		Line_id = psspy.abrnchar(-1, 0, 0, 1, 1, ['ID'])[1][0]


		if tp1:
			ierr = psspy.branch_chng_3(int(SCC_kl_array[tp1[0]][0]), int(SCC_kl_array[tp1[0]][1]),
								SCC_kl[tp1[0]][2], [0, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i],
								[psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
								 psspy._f, psspy._f, psspy._f, psspy._f],
								[psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
								 psspy._f, psspy._f, psspy._f, psspy._f], "")
			if ierr != 0:
				psspy.two_winding_chng_5(int(SCC_kl_array[tp1[0]][0]), int(SCC_kl_array[tp1[0]][1]),
										 SCC_kl[tp1[0]][2], [0, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i,
										  psspy._i, psspy._i, psspy._i, 0, psspy._i, psspy._i, psspy._i],
										 [psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
										  psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
										  psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f],
										 [psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
										  psspy._f, psspy._f, psspy._f, psspy._f, psspy._f], "", "")

		ierr = psspy.branch_chng_3(int(SCC_kl_array[SCC_ranking[i]][0]), int(SCC_kl_array[SCC_ranking[i]][1]),
							SCC_kl[SCC_ranking[i]][2], [0, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i],
							[psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
							 psspy._f, psspy._f, psspy._f, psspy._f],
							[psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
							 psspy._f, psspy._f, psspy._f, psspy._f], "")
		if ierr != 0:
			psspy.two_winding_chng_5(int(SCC_kl_array[SCC_ranking[i]][0]), int(SCC_kl_array[SCC_ranking[i]][1]),
							SCC_kl[SCC_ranking[i]][2], [0, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i,
									  psspy._i, psspy._i, psspy._i, 0, psspy._i, psspy._i, psspy._i],
									 [psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
									  psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
									  psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f],
									 [psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
									  psspy._f, psspy._f, psspy._f, psspy._f, psspy._f], "", "")

		ierr, buses = psspy.tree(1, -1)  # len(buses)=0 if there is no island, otherwise there is island
		psspy.tree(2, -1)  # need to call tree for 2nd time as required by this API, so program can move on

		if buses == 0:
			if not tp1:
				if not tp2:
					temp1 = np.absolute(Line_from - int(SCC_kl_array[SCC_ranking[i]][1])) + np.absolute(
						Line_to - int(SCC_kl_array[SCC_ranking[i]][0]))
					temp2 = np.absolute(Line_from - int(SCC_kl_array[SCC_ranking[i]][0])) + np.absolute(
						Line_to - int(SCC_kl_array[SCC_ranking[i]][1]))
					temp3 = np.multiply(np.sign(temp1), np.sign(temp2))
					idx = np.argmin(temp3)
					Z = np.absolute(Line_RX[idx])
					if Z > 1e-4:
						tp1.append(SCC_ranking[i])
			elif tp1:
				if tp2:
					break

				temp1 = np.absolute(Line_from - int(SCC_kl_array[SCC_ranking[i]][1])) + np.absolute(
					Line_to - int(SCC_kl_array[SCC_ranking[i]][0]))
				temp2 = np.absolute(Line_from - int(SCC_kl_array[SCC_ranking[i]][0])) + np.absolute(
					Line_to - int(SCC_kl_array[SCC_ranking[i]][1]))
				temp3 = np.multiply(np.sign(temp1), np.sign(temp2))
				idx = np.argmin(temp3)
				Z = np.absolute(Line_RX[idx])
				# if Z > 0.005:
				if Z > 1e-4:
					tp2.append(SCC_ranking[i])
				else:
					continue
	return tp1, tp2



def CalcSccN12(PFcase, busNlimit, POIi, OSW_id, SCC_kl, SCC_kl_array, tp1, tp2, all_POIs, Z9999_flag):
	# SCC at N-1
	psspy.psseinit(busNlimit)
	psspy.case(PFcase)
	psspy.short_circuit_coordinates(1)
	ierr = psspy.short_circuit_units(ival=1)
	psspy.progress_output(6,"",[0,0])

	ierr = psspy.branch_chng_3(int(SCC_kl_array[tp1[0]][0]), int(SCC_kl_array[tp1[0]][1]),
						SCC_kl[tp1[0]][2], [0, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i],
						[psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
						 psspy._f, psspy._f, psspy._f, psspy._f],
						[psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
						 psspy._f, psspy._f, psspy._f, psspy._f], "")
	if ierr!=0:
		psspy.two_winding_chng_5(int(SCC_kl_array[tp1[0]][0]), int(SCC_kl_array[tp1[0]][1]),
						SCC_kl[tp1[0]][2],
								 [0, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i,
								  psspy._i, psspy._i, 0, psspy._i, psspy._i, psspy._i],
								 [psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
								  psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
								  psspy._f, psspy._f, psspy._f, psspy._f, psspy._f],
								 [psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
								  psspy._f, psspy._f, psspy._f, psspy._f], "", "")
	SCC_POIi_N1, temp, temp = SC_k_lvl(psspy, POIi, OSW_id, 0, all_POIs, Z9999_flag)

	# SCC at N-2
	psspy.short_circuit_coordinates(1)
	psspy.psseinit(busNlimit)
	psspy.case(PFcase)
	psspy.short_circuit_coordinates(1)
	ierr = psspy.short_circuit_units(ival=1)
	psspy.progress_output(6,"",[0,0])
	ierr = psspy.branch_chng_3(int(SCC_kl_array[tp1[0]][0]), int(SCC_kl_array[tp1[0]][1]),
						SCC_kl[tp1[0]][2], [0, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i],
						[psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
						 psspy._f, psspy._f, psspy._f, psspy._f],
						[psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
						 psspy._f, psspy._f, psspy._f, psspy._f], "")
	if ierr!=0:
		psspy.two_winding_chng_5(int(SCC_kl_array[tp1[0]][0]), int(SCC_kl_array[tp1[0]][1]),
						SCC_kl[tp1[0]][2],
								 [0, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i,
								  psspy._i, psspy._i, 0, psspy._i, psspy._i, psspy._i],
								 [psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
								  psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
								  psspy._f, psspy._f, psspy._f, psspy._f, psspy._f],
								 [psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
								  psspy._f, psspy._f, psspy._f, psspy._f], "", "")

	ierr = psspy.branch_chng_3(int(SCC_kl_array[tp2[0]][0]), int(SCC_kl_array[tp2[0]][1]),
						SCC_kl[tp2[0]][2], [0, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i],
						[psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
						 psspy._f, psspy._f, psspy._f, psspy._f],
						[psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
						 psspy._f, psspy._f, psspy._f, psspy._f], "")
	if ierr!=0:
		psspy.two_winding_chng_5(int(SCC_kl_array[tp2[0]][0]), int(SCC_kl_array[tp2[0]][1]),
						SCC_kl[tp2[0]][2],
								 [0, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i,
								  psspy._i, psspy._i, 0, psspy._i, psspy._i, psspy._i],
								 [psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
								  psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
								  psspy._f, psspy._f, psspy._f, psspy._f, psspy._f],
								 [psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
								  psspy._f, psspy._f, psspy._f, psspy._f], "", "")

	SCC_POIi_N2, temp, temp = SC_k_lvl(psspy, POIi, OSW_id, 0, all_POIs, Z9999_flag)
	return SCC_POIi_N1[1], SCC_POIi_N2[1]
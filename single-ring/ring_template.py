#!/usr/bin/env python
#encoding=utf-8

"""
This python code for generating initial configuration of a single ring.
[USAGE] python ConfigGenerator.py [Input_Parameter_File] [Output_Config_File] [LAMMPS MOLECULE FILE] [check=on/off]
"""

import os, sys, math
import numpy as np
global l


if sys.argv[1] == "-h" or sys.argv[1] == "h" or sys.argv[1] == "help" or sys.argv[1] == "-help":
	print (__doc__)
	sys.exit()
	
FORMAT = "%.4f"	

try:
	ParameterIn = sys.argv[1]
	Outputfile = sys.argv[2]
	MoleFile = sys.argv[3]
	check = sys.argv[4]
except:
	print (__doc__)
	sys.exit()

if __name__ == '__main__':
	if len(sys.argv)==1:
		print (__doc__)
		sys.exit()
	else:
		In = open(ParameterIn,'r')
		Out = open(Outputfile,'w')

	print ("\n")
	print("#### [1] reading configuration parameters ####")
	
	lattPara=[]
	mass=[]
	angle=[]
	bond_1=[]
	pair_1=[]
	angle_1=[]
	dihedral_1=[]
	
	line = In.readline()
	while line:
		if line.startswith("Diameter"):
			strs = line.strip().split()
			Diameter = float(strs[1])
		if line.startswith("beeds_number"):
			strs = line.strip().split()
			beeds_number=float(strs[1])
		if line.startswith("Atom_Type"):
			strs = line.strip().split()
			Atomtype=float(strs[1])
			atom_mass=float(strs[2])
			mass.append(float(strs[2]))
		if line.startswith("Angle_Type"):
			strs = line.strip().split()
			AngleType = strs[1]
		if line.startswith("angle_parameter"):
			strs = line.strip().split()
			angle_1.append(strs[1])
			angle_1.append(int(strs[2]))
			angle_1.append(float(strs[3]))
		if line.startswith("Pair"):
			strs = line.strip().split()
			pair_1.append(strs[1])
			pair_1.append(float(strs[2]))
			pair_1.append(float(strs[3]))
			pair_1.append(float(strs[4]))
		if line.startswith("Bond_Type"):
			strs = line.strip().split()
			BondType = strs[1]
		if line.startswith("bond_parameter"):
			strs=line.strip().split()
			bond_1.append(strs[1])
			bond_1.append(int(strs[2]))
			bond_1.append(float(strs[3]))
		if line.startswith("Atom_Type_mass"):
			strs = line.strip().split()
			mass = strs[2]

		line = In.readline()


	Out.write("LAMMPS data file for hard-sphere chain with chiral SW interaction\n\n")
	Out.write(str(int(beeds_number))+" atoms\n")
	
	Out.write(str(int(beeds_number)*int(mass))+" mass\n")
	
	Out.write(str(0.0)+"	"+str(0.0)+"	"+str(0.0)+" com\n")
	
	Out.write(str(int(beeds_number))+" bonds\n")
	Out.write(str(int(beeds_number))+" angles\n")
	
	Out.write("\n")
	
	Out.write("Coords\n\n")
	i=0
	while i <beeds_number:
		Out.write(str(int(i+1))+"		"+str(0.0)+" "+str(Diameter/2/math.sin(math.pi/beeds_number)*math.sin(2*math.pi/beeds_number*(i+1)))+" "+str(Diameter/2/math.sin(math.pi/beeds_number)*math.cos(2*math.pi/beeds_number*(i+1)))+"\n")
		i=i+1
	Out.write("\n\n")
	
	Out.write("Types\n\n")	
	i=1 
	while i<beeds_number+1:
		Out.write(str(i)+"		"+str(1)+"\n")
		i=i+1
	Out.write("\n")
	
	Out.write("Bonds\n\n")
	i=1
	while i<beeds_number:
		Out.write(str(i)+"		"+str(1)+"		"+str(i)+"		"+str(i+1)+"\n")
		i=i+1
	if i==beeds_number:
		Out.write(str(int(beeds_number))+"		"+str(1)+"		"+str(int(beeds_number))+"		"+str(1)+"\n")
	Out.write("\n")
	
	Out.write("Angles\n\n")
	i=1
	while i<beeds_number-1:
		Out.write(str(i)+"		"+str(1)+"		"+str(i)+"		"+str(i+1)+"		"+str(i+2)+"\n")
		i=i+1
	if i==beeds_number-1:
		Out.write(str(int(beeds_number-1))+"		"+str(1)+"		"+str(int(beeds_number-1))+"		"+str(int(beeds_number))+"		"+str(1)+"\n")
		i=i+1
	if i==beeds_number:
		Out.write(str(int(beeds_number))+"		"+str(1)+"		"+str(int(beeds_number))+"		"+str(1)+"		"+str(2)+"\n")
		i=i+1
	Out.write("\n")
	
	# Out.write("Shake Flags\n\n")
	# i=1
	# while i<beeds_number+1:
		# Out.write(str(i)+"		"+str(1)+"\n")
		# i=i+1
	# Out.write("\n")
	
	# Out.write("Shake Atoms\n\n")
	# i=1
	# while i<beeds_number+1:
		# Out.write(str(i)+"		"+str(1)+"\n")
		# i=i+1
	# Out.write("\n")
	
	# Out.write("Shake Bonds Types\n\n")
	# i=1
	# while i<beeds_number+1:
		# Out.write(str(i)+"		"+str(1)+"		"+str(1)+"		"+str(1)+"\n")
		# i=i+1
	# Out.write("\n")
	
	# Out.write("Special Bond Counts\n\n")
	# i=1
	# while i<beeds_number+1:
		# Out.write(str(i)+"		"+str(2)+"		"+str(2)+"		"+str(2)+"\n")
		# i=i+1
	# Out.write("\n")
	
	# Out.write("Special Bonds\n\n")
	# i=1
	# if i==1:
		# Out.write(str(1)+"		"+str(2)+"		"+str(28)+"\n")
		# i=i+1
	# while i<beeds_number:
		# Out.write(str(i)+"		"+str(i-1)+"		"+str(i+1)+"\n")
		# i=i+1
	# if i==beeds_number:
		# Out.write(str(28)+"		"+str(27)+"		"+str(1)+"\n")
	# Out.write("\n")
	
	Out.write("Diameters\n\n")
	i=1
	while i<beeds_number+1:
		Out.write(str(i)+"		"+str(Diameter)+"\n")
		i=i+1
	Out.write("\n")
	
	Out.write("Masses\n\n")
	i=1
	while i<beeds_number+1:
		Out.write(str(i)+"		"+str(mass)+"\n")
		i=i+1
	Out.write("\n")
	
	


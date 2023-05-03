#环境配制等前期工作
import os, sys, math
import numpy as np

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

#charge?
#分子坐标 n表示每个圆环中的beeds数  这里我们考虑在z平面扩展成圆环 sig为原子直径
def Molecule(Coor,sig,n):
	i=1
	MoleXYZ=[]
	Atomtype=1
	Charge=0.0
	while i<=2.0*n:
		if i<=n:
			x=Coor[0]
			y=Coor[1]+sig/2/math.sin(math.pi/n)*math.sin(2*math.pi/n*i)
			z=Coor[2]+sig/2/math.sin(math.pi/n)*math.cos(2*math.pi/n*i)
		else:
			x=Coor[0]+sig/2/math.sin(math.pi/n)*math.sin(2*math.pi/n*i)
			y=Coor[1]
			z=Coor[2]+sig/math.sin(math.pi/n)/2+sig/2/math.sin(math.pi/n)*math.cos(2*math.pi/n*i)
		MoleXYZ.append((i,Atomtype,Charge,float(x),float(y),float(z)))
		i+=1
	return MoleXYZ

# 晶格构建 nac为配位数 eta（packing fraction）
# 晶格构建 nac为配位数 eta（packing fraction）
def Lattice(nac,ncx,ncy,ncz,Vmol,eta,lds):
	global nPart

	nPart = nac*ncx*ncy*ncz # nPart is number of molecules
	print ("#### Coordinates of ", nPart, "molecules created by this script ####")

	root2 = 2.0**(0.5)
	root3 = 3.0**(0.5)
	# if nac == 4:
		# zstep = 2.0*(root2/root3+lds) #fcc 
	# elif nac == 6:
		# zstep = 3.0*(root2/root3+lds) #hpc

	# ux = [5.0, 15.0, 15.0, 5.0, 15.0, 5.0]
	# uz = []; uy = []
	
	# uz.append(lds/2.0 + 1.0/root2/root3)
	# uz.append(lds/2.0 + 1.0/root2/root3)
	# uz.append(3.0*lds/2.0 + root3/root2)
	# uz.append(3.0*lds/2.0 + root3/root2)
	# uz.append(5.0*lds/2.0 + 5.0/root2/root3)
	# uz.append(5.0*lds/2.0 + 5.0/root2/root3)
	
	# uy.append(3.0*root3/4.0)
	# uy.append(3.0*root3/4.0)
	# uy.append(3.0*root3/4.0)
	# uy.append(3.0*root3/4.0)
	# uy.append(3.0*root3/4.0)
	# uy.append(3.0*root3/4.0)
	
	


	# lbox = []
	# lbox.append(1.0)
	# lbox.append(1.0*root3*float(ncy)/float(ncx))
	# lbox.append(1.0*zstep*float(ncz)/float(ncx))
	root2 = 2.0**(0.5)
	root3 = 3.0**(0.5)

	ux = [0.5]
	uz = []; uy = []
	
	uz.append(0.5)

	uy.append(0.5)
	
	lbox = []
	lbox.append(1.0)
	lbox.append(1.0)
	lbox.append(1.0)
	l = (nPart*Vmol/lbox[1]/lbox[2]/eta)**(1.0/3.0)
	print ("#### The dimensions of the lattice created are (x/y/z):", lbox[0]*l,lbox[1]*l,lbox[2]*l)

	print ("#### creating fcc/hcp lattice for your system, please wait.... ####")
	i=0
	rPart = []
	
	for j in range(nac):
		for xx in range(ncx):
			for yy in range(ncy):
				for zz in range(ncz):
					i+=1
					#1.02 is scaling factor assigned to Y-axis
					# rPart.append(((float(xx)*20+ux[j])/float(ncx+1)*l-lbox[0]*0.5, \
								  # (root3*float(yy)*10+uy[j])/float(ncx+1)*l-lbox[1]*0.5,\
								  # (zstep*float(zz)+uz[j])/float(ncx+1)*l-lbox[2]*0.5  \
								# ))
					rPart.append((float(xx)*l+l*ux[j], \
								  float(yy)*l+l*uy[j],\
								  float(zz)*l+l*uz[j] \
								))
					if len(rPart)%500 == 0:
						print ("####", i,"-th Molecule, Lattice Coordinates are", rPart[len(rPart)-1], "####")
						print( "####", FORMAT%(float(len(rPart))/float(nPart)*100), "% is finished, please waiting.... ####")

	return rPart

#参数导入
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

	line = In.readline()
	while line:
		if line.startswith("Diameter"):
			strs = line.strip().split()
			sig = float(strs[1])
		if line.startswith("beeds_number"):
			strs = line.strip().split()
			beeds_number=float(strs[1])
		if line.startswith("Atom_Type"):
			strs = line.strip().split()
			Atomtype=float(strs[1])
			atom_mass=float(strs[2])
			mass.append(float(strs[2]))
		if line.startswith("Packing_fraction"):
			strs = line.strip().split()
			eta = float(strs[1])
		if line.startswith("Lattice_Parameter"):
			strs = line.strip().split()
			lattPara.append(int(strs[1]))
			lattPara.append(int(strs[2]))
			lattPara.append(int(strs[3]))
			lattPara.append(int(strs[4]))
		if line.startswith("Lattice_Relax"):
			strs = line.strip().split()
			lattRex = float(strs[1])
		if line.startswith("Accuracy"):
			strs = line.strip().split()
			Nacc = int(strs[1])
		if line.startswith("rcut"):
			strs = line.strip().split()
			rcut = float(strs[1])
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
		line = In.readline()


	Vmol = 4.0/3.0*(sig/2.0)**3*math.pi*beeds_number*2

	iMole = 0
	nMole = beeds_number*2 #atoms per molecule

	MoleXYZ0 = []; MoleXYZ = []
	ConfigXYZ = []; ConfigX = []; ConfigY = []; ConfigZ = []
	Coor=[]

	print ("#### [3] Creating lattice for the system ####")
	# calculating the location of the 1st atom of a molecule in fcc/hcp lattice
	MoleXYZ0 = Lattice(lattPara[0],lattPara[1],lattPara[2],lattPara[3],Vmol,eta,sig*beeds_number)

	MolNum=nPart
	print("#### [4] Creating Coordinates for all atoms in the system ####")
	while iMole<MolNum:
		Coor.append(MoleXYZ0[iMole][0])
		Coor.append(MoleXYZ0[iMole][1])
		Coor.append(MoleXYZ0[iMole][2])
		MoleXYZ=Molecule(Coor,sig,beeds_number)
		for iAtom in range(int(nMole)):
			ConfigXYZ.append((MoleXYZ[iAtom][0]+iMole*nMole,\
							iAtom+1,\
							MoleXYZ[iAtom][1],\
							MoleXYZ[iAtom][2],\
							MoleXYZ[iAtom][3],\
							MoleXYZ[iAtom][4],\
							MoleXYZ[iAtom][5]))
			ConfigX.append(MoleXYZ[iAtom][3])
			ConfigY.append(MoleXYZ[iAtom][4])
			ConfigZ.append(MoleXYZ[iAtom][5])	
		iMole =iMole+1
		Coor=[]
		
	xlo = FORMAT%(min(ConfigX));	xhi = FORMAT%(max(ConfigX))
	ylo = FORMAT%(min(ConfigY));	yhi = FORMAT%(max(ConfigY))	
	zlo = FORMAT%(min(ConfigZ));	zhi = FORMAT%(max(ConfigZ))
	
	i=0; shiftXYZ=[]
	while i<nMole*MolNum:
		shiftXYZ.append((float(ConfigXYZ[i][4])-float(xlo),\
						 float(ConfigXYZ[i][5])-float(ylo),\
						 float(ConfigXYZ[i][6])-float(zlo)))
		i+=1
	xlo=float(xlo)-float(xlo); xhi=float(xhi)-float(xlo)
	ylo=float(ylo)-float(ylo); yhi=float(yhi)-float(ylo)
	zlo=float(zlo)-float(zlo); zhi=float(zhi)-float(zlo)
	max=max(xhi,yhi,zhi)
	print ("#### The size of simulation box is ", "[X:",xlo,"|",xhi,"]",\
			"[Y:",ylo,"|",yhi,"]", "[Z:",zlo,"|",zhi, "] ####")

			
	Bond=[];Angle=[]
	nBond = (beeds_number)*MolNum*2
	nAngle = MolNum*beeds_number*2
	
	
	print ("#### [5] Creating tables for Bonds and Angles ####")
	i=1;iMole=0;iAtom=0   #i 表示bond number
	while iMole < MolNum:
		for iAtom in range (0,int(nMole)):
			BondType=1
			if iAtom==int(nMole)/2-1 :
				Bond.append((int(i),BondType,int(iAtom+1+iMole*nMole),int(1+iMole*nMole)))   
				i+=1
			elif iAtom==int(nMole)-1:
				Bond.append((int(i),BondType,int(iAtom+1+iMole*nMole),int(int(nMole)/2+1+iMole*nMole)))
				i+=1
			elif iAtom<int(nMole)/2-1 or (iAtom<int(nMole)-1 and iAtom>int(nMole)/2-1) :
				Bond.append((int(i),BondType,int(iAtom+1+iMole*nMole),int(iAtom+2+iMole*nMole)))
				i+=1	
		iAtom=0
		iMole+=1
		
		
	i=1;iMole=0;iAtom=0
	while iMole <MolNum:
		for iAtom in range(0,int(nMole)):
			AngleType=1
			if iAtom==0:
				Angle.append((int(i),AngleType,int(int(nMole)/2+iMole*nMole),int(iAtom+iMole*nMole+1),int(iAtom+iMole*nMole+2)))
				i+=1
			elif iAtom==int(nMole)/2-1 :
				Angle.append((int(i),AngleType,int(iAtom+iMole*nMole),int(iAtom+iMole*nMole+1),int(iMole*nMole+1)))
				i+=1
			elif iAtom==int(nMole)/2:
				Angle.append((int(i),AngleType,int(int(nMole)+iMole*nMole),int(iAtom+iMole*nMole+1),int(iAtom+iMole*nMole+2)))
				i+=1
			elif iAtom==int(nMole)-1:
				Angle.append((int(i),AngleType,int(iAtom+iMole*nMole),int(iAtom+iMole*nMole+1),int(int(nMole)/2+1+iMole*nMole)))
				i+=1
			else:
				Angle.append((int(i),AngleType,int(iAtom+iMole*nMole),int(iAtom+iMole*nMole+1),int(iAtom+iMole*nMole+2)))
				i+=1
		iAtom=0
		iMole+=1
		
	print ("#### [6] Writing to LAMMPS DATA FILE ####.")
	Out.write("LAMMPS data file for hard-sphere chain with chiral SW interaction\n\n")
	Out.write(str(int(MolNum*nMole))+" atoms\n")
	Out.write(str(int(nBond))+" bonds\n")
	Out.write(str(int(nAngle))+" angles\n")
	Out.write("0 impropers\n")
	Out.write("\n")

	Out.write(str(1)+"  atom types\n")
	Out.write(str(BondType)+"  bond types\n")
	Out.write(str(AngleType)+"  angle types\n")
	Out.write("0  improper types\n")
	Out.write("\n")	

	Out.write(str(xlo) + "	" + str(int(max)) +"		xlo xhi\n")
	Out.write(str(ylo) + "	" + str(int(max)) +"		ylo yhi\n")		
	Out.write(str(zlo) + "	" + str(int(max)) +"		zlo zhi\n")
	Out.write("\n")
	
	Out.write("Masses\n\n")
	i=0
	while i<1:
		Out.write(str(i+1)+"	"+str(mass[i])+" \n" )
		i+=1
	Out.write("\n")
	
	Out.write("Bond Coeffs\n\n")
	Out.write(str(1)+"	"+str(bond_1[1])+"	"+str(bond_1[2])+" \n")
	Out.write("\n")
	
	Out.write("Angle Coeffs\n\n")
	Out.write(str(1)+"	"+str(angle_1[1])+"	"+str(angle_1[2])+" \n")
	Out.write("\n")
	
	Out.write("Atoms\n\n")
	i=0 	
	while i<nMole*MolNum:
		#ConfigXYZ[i][1]=int(ConfigXYZ[i][1])-28
		Out.write(str(int(ConfigXYZ[i][0]))+"		"+str(int(i//beeds_number+1))+"		"+\
				str(ConfigXYZ[i][2])+"		"+str(ConfigXYZ[i][3])+"		"+\
				str(FORMAT%(shiftXYZ[i][0]))+"		"+\
				str(FORMAT%(shiftXYZ[i][1]))+"		"+\
				str(FORMAT%(shiftXYZ[i][2]))+"\n")
#				str(FORMAT%(ConfigXYZ[i][4]))+"		"+\
#				str(FORMAT%(ConfigXYZ[i][5]))+"		"+str(FORMAT%(ConfigXYZ[i][6]))+"	\n")
		
		i+=1
	Out.write("\n")
	
	Out.write("Bonds\n\n")
	i=0
	while i<nBond:
		Out.write("		"+str(Bond[i][0])+"		"+str(Bond[i][1])+"		"+str(Bond[i][2])+"		"+\
				str(Bond[i][3])+"	\n")
		i+=1
	Out.write("\n")
	
	Out.write("Angles\n\n")
	i=0
	while i<nAngle:
		Out.write("		"+str(Angle[i][0])+"		"+str(Angle[i][1])+"		"+str(Angle[i][2])+"		"+\
				str(Angle[i][3])+"		"+str(Angle[i][4])+"	\n")
		i+=1
	Out.write("\n")	
	
	print ("#### [END] LAMMPS DATA FILE is successfully created. [END] ###")
	
# *********************************************************************************
	
	print ("#### Creating LAMMPS MOLECULE FILE ")
	MoleOut = open(MoleFile, 'w')
	MoleOut.write("# LAMMPS MOLECULE FILE\n\n")
	MoleOut.write(str(int(nMole))+"		atoms\n")
	MoleOut.write(str(int(nBond/MolNum))+"		bonds\n")
	MoleOut.write(str(int(nAngle/MolNum))+"		angles\n")	
	MoleOut.write(str(0)+"		impropers\n")
	MoleOut.write(str(beeds_number*atom_mass)+"		mass\n")

	Coor=[0,0,0]
	MoleCoor=Molecule(Coor,sig,beeds_number)
	
	tx=0.0; ty=0.0; tz=0.0
	tmass=beeds_number*atom_mass
	for aa in MoleCoor:
		tx=tx+float(aa[2]*atom_mass)
		ty=ty+float(aa[3]*atom_mass)
		tz=tz+float(aa[4]*atom_mass)
	tx = FORMAT % (float(tx)/tmass)
	ty = FORMAT % (float(ty)/tmass)
	tz = FORMAT % (float(tz)/tmass)	
	MoleOut.write(str(tx)+"		"+str(ty)+"		"+str(tz)+"		com\n")
		
	Ixx=0.0; Iyy=0.0; Izz=0.0
	Ixy=0.0; Ixz=0.0; Iyz=0.0
	for aa in MoleCoor:
		Ixx = Ixx + atom_mass*(float(aa[4])*float(aa[4])+float(aa[3])*float(aa[3]))
		Iyy = Iyy + atom_mass*(float(aa[2])*float(aa[2])+float(aa[4])*float(aa[4]))
		Izz = Izz + atom_mass*(float(aa[3])*float(aa[3])+float(aa[2])*float(aa[2]))
		Ixy = Ixy + -1.0*atom_mass*(float(aa[2])*float(aa[3]))
		Ixz = Ixz + -1.0*atom_mass*(float(aa[2])*float(aa[4]))
		Iyz = Iyz + -1.0*atom_mass*(float(aa[4])*float(aa[3]))
	MoleOut.write(str(FORMAT % (Ixx))+"	"+str(FORMAT % (Iyy))+"	"+str(FORMAT % (Izz))+\
				"	"+str(FORMAT % (Ixy))+"	"+str(FORMAT % (Ixz))+"	"+str(FORMAT % (Iyz))+" inertia\n")
	
	
	MoleOut.write("\n")
	MoleOut.write("Coords\n")
	for aa in MoleCoor:
		MoleOut.write(str(aa[0])+"	"+str(FORMAT%(aa[2]))+"		"+str(FORMAT%(aa[3]))+"		"+str(FORMAT%(aa[4]))+"\n")
	
	MoleOut.write("\n")
	MoleOut.write("Diameters\n")
	for aa in MoleCoor:
		MoleOut.write(str(aa[0])+"	"+str(FORMAT%(sig))+"\n")
	
	MoleOut.write("\n")
	MoleOut.write("Masses\n")
	for aa in MoleCoor:
		MoleOut.write(str(aa[0])+"	"+str(FORMAT%(atom_mass))+"\n")
	
	MoleOut.write("\n")
	MoleOut.write("Bonds\n")
	i=0
	while i<nBond/MolNum:
		MoleOut.write("		"+str(Bond[i][0])+"		"+str(Bond[i][1])+"		"+str(Bond[i][2])+"		"+\
				str(Bond[i][3])+"	\n")
		i+=1
	
	MoleOut.write("\n")
	MoleOut.write("Angles\n")		
	i=0
	while i<nAngle/MolNum:
		MoleOut.write("		"+str(Angle[i][0])+"		"+str(Angle[i][1])+"		"+str(Angle[i][2])+"		"+\
				str(Angle[i][3])+"		"+str(Angle[i][4])+"	\n")
		i+=1	
		
#	Calculating distance between a atom pair
	if check == "on" :
		print ("#### Calculating distance between a atom pair ####")
		rij=[];	iAtom=0; jAtom=0; i=0
		for iAtom in range(int(MolNum*nMole)):
			for jAtom in range(int(MolNum*nMole)): 
				if iAtom != jAtom and ConfigXYZ[iAtom][1]!=ConfigXYZ[jAtom][1] and (iAtom%2!=0 and jAtom%2!=0):
					rx=(ConfigXYZ[iAtom][4]-ConfigXYZ[jAtom][4])**2.0
					ry=(ConfigXYZ[iAtom][5]-ConfigXYZ[jAtom][5])**2.0
					rz=(ConfigXYZ[iAtom][6]-ConfigXYZ[jAtom][6])**2.0
					if (rx+ry+rz)**0.5<1.0:
						rij.append((iAtom,jAtom)); i+=1
						print (i,"	", iAtom, jAtom, "in bad configuration", FORMAT%(rx) ,FORMAT%(ry), FORMAT%(rz))
				
				
		
		
		
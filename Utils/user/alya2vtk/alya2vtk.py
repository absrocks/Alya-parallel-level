#!/usr/bin/env python

import sys

def usage():
	print "                      ___    __           ___ _    ____  __     "
	print "                     /   |  / /_  ______ |__ \ |  / / /_/ /__   "
	print "                    / /| | / / / / / __ `/_/ / | / / __/ //_/   "
	print "                   / ___ |/ / /_/ / /_/ / __/| |/ / /_/ ,<      "   
	print "                  /_/  |_/_/\__, /\__,_/____/|___/\__/_/|_|     "  
	print "                           /____/                               "
	print " ---------------------------------------------------------------------------------"
	print "| Mesh convertor from Alya format to VTK (ASCII) format                           |"
	print "| Whit this tool you don't need to execute Alya to visualize the mesh in Paraview |"
	print "| So far it only works for 3D meshes with tetrahedral elements                    |"
	print "| But I have plans to extend it for 2D and different type of elements             |"
	print "| Usage: $./alya2vtk.py input_file or $python alya2vtk.py input_file              |"
	print "| The input_file must include the keywords in the first column and                |"
	print "| the name of each file in the second one. Check the examples!!!                  |"
	print "| Example:                                                                        |"
	print "| nodes_file nodes.dom.dat (mandatory)                                            |"
	print "| elements_file elem.dom.dat (mandatory)                                          |"
	print "| boundaries_file bounda.dom.dat                                                  |"
	print "| codes_file codes.dom.dat                                                        |"
	print "| materials_file mats.dat                                                         |"
	print "| vectors_file fibers.dat                                                         |"
	print "| outfile_name meshinvtk.vtk                                                      |"
	print "| STL (if you want to plot STL's, but is in testing mode)                         |"
	print "| (As you can realize the data must be splitted in different files)               |"
	print "| The keywords are nodes_file, elements_file, boundaries_file,                    |"
	print "| codes_file, materials_file and output_file (and STL, but is in beta testing)    |" 
	print "| Additional info: data must be sorted and in the case of materials,              |"
	print "| each cell has to have its own material (i.e. zero if none)                      |"
	print "| When using the boundary file, it result will only be visible using              |"
	print "| clips in paraview, because this option allows you to visualize the              |"
	print "| boundary definition in a different color. So if everything went OK,             |" 
	print "| you won't be able to see the diference of colours unless you clip the volume.   |"
	print "| One last thing: codes_file allows you to plot the boundary sets. i.e. you can   |"
	print "| see using colours the boundary set you have defined.                            |"
	print "| One more time: check the examples!!!!!                                          |"
	print "| If you have any problem using this converter,                                   |"
	print "| or if you really really need a 2D version or support for other type of element  |"
	print "| please send me an email: matias.rivero@bsc.es or call me!                       |"
	print " ---------------------------------------------------------------------------------"
	sys.exit()
	#TODO list: 
	#1- Corregir el problema con el encabezado ON_BOUNDARIES en codes_file. Actualmente no funciona con el encabezado.
	#2- Generalizarlo a 2D y para cualquier tipo de elemento (tanto en 3D como 2D).

nodesfile_present = False
elemfile_present = False
boundfile_present = False
codesfile_present = False
matfile_present = False
vecfile_present = False
outfile_present = False
STL_present = False #maneje para que raul pueda sacar STLs
main_list = []
main_list2 = []

if (len(sys.argv) < 2):
	usage()
elif (sys.argv[1] == '-v') | (sys.argv[1] == '--version') | (sys.argv[1] == '-version'):
	usage()

try:
	input_file = open(sys.argv[1],'r')
except:
	print "Couldn't open input_file. Something is wrong..."
	sys.exit()
for line in input_file:
	line = line.strip()
	if line.startswith('nodes_file'):
		line = line.split()
		nodesfile_name = line[1]
		nodesfile_present = True
		try:
			tmp = open(nodesfile_name,'r') 
		except:
			print "Couldn't open nodes_file. Check name..."
			sys.exit()
		first_line = tmp.readline()
		first_line = first_line.strip()
		if first_line.startswith('COORD'): #check if the first line is the header
			nodesfile_header = True
			offset_nodes = 1
		else: 
			nodesfile_header = False
			offset_nodes = 0
		tmp.close()
	elif line.startswith('elements_file'):
		line = line.split()
		elemfile_name = line[1]
		elemfile_present = True
		try:
			tmp = open(elemfile_name,'r')
		except:
			print "Couldn't open elements_file. Check name..."
			sys.exit()
		first_line = tmp.readline()
		first_line = first_line.strip()
		if first_line.startswith('ELEME'):
			elemfile_header = True
			offset_eleme = 1
		else: 
			elemfile_header = False
			offset_eleme = 0			
		tmp.close()
	elif line.startswith('boundaries_file'):
		line = line.split()
		boundfile_name = line[1]
		boundfile_present = True
		try:
			tmp = open(boundfile_name,'r')
		except:
			print "Couldn't open boundaries_file. Check name..."
			sys.exit()
		first_line = tmp.readline()
		first_line = first_line.strip()
		if first_line.startswith('BOUNDA'):
			boundfile_header = True
			offset_bounda = 1			
		else: 
			boundfile_header = False
			offset_bounda = 0		
		tmp.close()
	elif line.startswith('codes_file'):
                if not boundfile_present:
                	print "No boundary file present - must be included if you want to plot codes"
                        sys.exit()
		line = line.split()
		codesfile_name = line[1]
		codesfile_present = True
		try:
			tmp = open(codesfile_name,'r')
		except:
			print "Couldn't open codes_file. Check name..."
			sys.exit()
		first_line = tmp.readline()
		first_line = first_line.strip()
		if first_line.startswith('ON_BOUNDA'):
			codesfile_header = True
			offset_codes = 1			
		else: 
			codesfile_header = False
			offset_codes = 0		
		tmp.close()
	elif line.startswith('materials_file'):
		line = line.split()
		matfile_name = line[1]
		matfile_present = True
		try:
			tmp = open(matfile_name,'r')
		except:
			print "Couldn't open materials_file. Check name..."
			sys.exit()
		first_line = tmp.readline()
		first_line = first_line.strip()
		if first_line.startswith('MATER'):
			matfile_header = True
			offset_mater = 1		
		else: 
			matfile_header = False
			offset_mater = 0		
		tmp.close()
	elif line.startswith('vectors_file'):
		line = line.split()
		vecfile_name = line[1]
		vecfile_present = True
		try:
			tmp = open(vecfile_name,'r')
		except:
			print "Couldn't open vectors_file. Check name..."
			sys.exit()
		first_line = tmp.readline()
		first_line = first_line.strip()
		if first_line.startswith('FIELD'):
			vecfile_header = True
			offset_vector = 1		
		else: 
			vecfile_header = False
			offset_vector = 0		
		tmp.close()
	elif line.startswith('output_file'):
		line = line.split()
		outfile_name = line[1]
		outfile_present = True			
	elif line.startswith('STL'):
		STL_present = True
input_file.close()	


if not nodesfile_present:
	print "Nodes file name is missing, bye..."
	sys.exit()	
if not elemfile_present:
	print "Elements file name is missing, bye..."
	sys.exit()
if outfile_present:
	vtkfile = open(outfile_name,'w')
else:
	print "Output file name is missing, bye..."
	sys.exit()	
if not boundfile_present:
	print "No boundaries file present"
if not codesfile_present:
	print "No codes file present"
if not matfile_present:	
	print "No materials file present"
if not vecfile_present:	
	print "No vectors file present"


#file de coordenadas
print "computing coordinates file..."
nodesfile = open(nodesfile_name,'r')
num_nodes = 0
for line in nodesfile:
	line = line.strip()
	main_list.append(line.split())
	num_nodes = num_nodes + 1
nodesfile.close()
if nodesfile_header:
	num_nodes = num_nodes - 2

vtkfile.write('# vtk DataFile Version 2.0\n')
vtkfile.write('Generated with Alya2Vtk\n')
vtkfile.write('ASCII\n')
vtkfile.write('DATASET UNSTRUCTURED_GRID\n')	
vtkfile.write('POINTS')
vtkfile.write(' ')
vtkfile.write(str(num_nodes))
vtkfile.write(' ')
vtkfile.write('FLOAT \n')
for i in range(num_nodes):
	vtkfile.write(str(main_list[i+offset_nodes][1]))
	vtkfile.write(' ')
	vtkfile.write(str(main_list[i+offset_nodes][2]))
	vtkfile.write(' ')	
	vtkfile.write(str(main_list[i+offset_nodes][3]))
	vtkfile.write('\n')	

del main_list[:]
print "done!"


#file de elementos
print "computing elements file..."
elemfile = open(elemfile_name,'r')
num_elem = 0
for line in elemfile:
	line = line.strip()
	main_list.append(line.split())
	num_elem = num_elem + 1
elemfile.close()
if elemfile_header:
	num_elem = num_elem - 2

vtkfile.write('CELLS')
if not STL_present:
	vtkfile.write(' ')
	vtkfile.write(str(num_elem))
	vtkfile.write(' ')
	vtkfile.write(str(num_elem*5))
	vtkfile.write('\n')
	for i in range(num_elem):
		vtkfile.write('4')
		vtkfile.write(' ')
		vtkfile.write(str(int(main_list[i+offset_eleme][1])-1))
		vtkfile.write(' ')
		vtkfile.write(str(int(main_list[i+offset_eleme][2])-1))
		vtkfile.write(' ')	
		vtkfile.write(str(int(main_list[i+offset_eleme][3])-1))
		vtkfile.write(' ')	
		vtkfile.write(str(int(main_list[i+offset_eleme][4])-1))
		vtkfile.write('\n')	
	vtkfile.write('CELL_TYPES')
	vtkfile.write(' ')
	vtkfile.write(str(num_elem))
	vtkfile.write('\n')
	for i in range(num_elem):
		vtkfile.write('10 \n')
if STL_present:
	vtkfile.write(' ')
	vtkfile.write(str(num_elem))
	vtkfile.write(' ')
	vtkfile.write(str(num_elem*4))
	vtkfile.write('\n')
	for i in range(num_elem):
		vtkfile.write('3')
		vtkfile.write(' ')
		vtkfile.write(str(int(main_list[i+offset_eleme][1])-1))
		vtkfile.write(' ')
		vtkfile.write(str(int(main_list[i+offset_eleme][2])-1))
		vtkfile.write(' ')	
		vtkfile.write(str(int(main_list[i+offset_eleme][3])-1))
		vtkfile.write('\n')	
	vtkfile.write('CELL_TYPES')
	vtkfile.write(' ')
	vtkfile.write(str(num_elem))
	vtkfile.write('\n')
	for i in range(num_elem):
		vtkfile.write('5 \n')

del main_list[:]
print "done!"


#INPUTS OPTATIVOS
#POINT DATA
#file de boundaries
if boundfile_present:
        print "computing boundaries file..."
	boundfile = open(boundfile_name,'r')
	num_bound = 0
	for line in boundfile:
		line = line.strip()
		main_list.append(line.split())
		num_bound = num_bound + 1
	boundfile.close()
	if boundfile_header:
		num_bound = num_bound - 2
	tmp = [0 for i in range(num_nodes)]
	for i in range(num_bound):
		tmp[int(main_list[i+offset_bounda][1])-1] = 1
		tmp[int(main_list[i+offset_bounda][2])-1] = 1
		tmp[int(main_list[i+offset_bounda][3])-1] = 1

	vtkfile.write('POINT_DATA')
	vtkfile.write(' ')
	vtkfile.write(str(num_nodes))
	vtkfile.write('\n')
	vtkfile.write('SCALARS')
	vtkfile.write(' ')
	vtkfile.write('Boundaries')
	vtkfile.write(' ')
	vtkfile.write('INT\n')
	vtkfile.write('LOOKUP_TABLE')
	vtkfile.write(' ')
	vtkfile.write('default\n')
	for i in range(num_nodes):
		vtkfile.write(str(tmp[i]))
		vtkfile.write('\n')
	print "done!"


#file de codes
if codesfile_present:
	print "computing codes file..."
	codfile = open(codesfile_name,'r')
	for line in codfile:
		line = line.strip()
		main_list2.append(line.split())
	tmp = [0 for i in range(num_nodes)]

	codfile.seek(0)
	
	i = 0
	for line in codfile:
		tmp[int(main_list[int(main_list2[i+offset_codes][0])-1+offset_bounda][1])-1] = int(main_list2[i+offset_codes][1])
		tmp[int(main_list[int(main_list2[i+offset_codes][0])-1+offset_bounda][2])-1] = int(main_list2[i+offset_codes][1])
		tmp[int(main_list[int(main_list2[i+offset_codes][0])-1+offset_bounda][3])-1] = int(main_list2[i+offset_codes][1])
		i = i + 1
	codfile.close()

	vtkfile.write('SCALARS')
	vtkfile.write(' ')
	vtkfile.write('Codes')
	vtkfile.write(' ')
	vtkfile.write('INT\n')
	vtkfile.write('LOOKUP_TABLE')
	vtkfile.write(' ')
	vtkfile.write('default\n')
	for i in range(num_nodes):
		vtkfile.write(str(tmp[i]))
		vtkfile.write('\n')
	print "done!"

	del main_list2[:] #delete main 
del main_list[:]


#file de vectores
if vecfile_present:
        print "computing vectors file..."
	vecfile = open(vecfile_name,'r')
	for line in vecfile:
		line = line.strip()
		main_list.append(line.split())
	vecfile.close()

	if not boundfile_present:
		vtkfile.write('POINT_DATA')
		vtkfile.write(' ')
		vtkfile.write(str(num_nodes))
		vtkfile.write('\n')
	vtkfile.write('VECTORS')
	vtkfile.write(' ')
	vtkfile.write('Vectors')
	vtkfile.write(' ')
	vtkfile.write('FLOAT\n')
	for i in range(num_nodes):
		vtkfile.write(str(main_list[i+offset_vector][1]))
		vtkfile.write(' ')
		vtkfile.write(str(main_list[i+offset_vector][2]))
		vtkfile.write(' ')
		vtkfile.write(str(main_list[i+offset_vector][3]))
		vtkfile.write(' ')
		vtkfile.write('\n')

	del main_list[:]
	print "done!"


#file de materiales
#CELL DATA
if matfile_present:
        print "computing materials file..."
	matfile = open(matfile_name,'r')
	for line in matfile:
		line = line.strip()
		main_list.append(line.split())
	matfile.close()

	vtkfile.write('CELL_DATA')
	vtkfile.write(' ')
	vtkfile.write(str(num_elem))
	vtkfile.write('\n')
	vtkfile.write('SCALARS')
	vtkfile.write(' ')
	vtkfile.write('Material')
	vtkfile.write(' ')
	vtkfile.write('INT\n')
	vtkfile.write('LOOKUP_TABLE')
	vtkfile.write(' ')
	vtkfile.write('default\n')
	for i in range(num_elem):
		vtkfile.write(str(main_list[i+offset_mater][1]))
		vtkfile.write('\n')

	del main_list[:]
	print "done!"

vtkfile.close()

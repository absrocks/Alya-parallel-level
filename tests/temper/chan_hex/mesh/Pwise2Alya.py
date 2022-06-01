# Loads files from PointWise to Alya format, including case of 2nd order

import os
import sys
import numpy as np
import array as ar
import copy

#os.environ['PYTHONINSPECT'] = 'TRUE'
caseName = sys.argv[1]

# read linear mesh file (gmsh format 2.2.08)

linearMeshName = sys.argv[2]
linearMeshFile = '{}.msh'.format(linearMeshName)
f = open(linearMeshFile,'r')

line = True
inNODE = False
inELEM = False
inPhNa = False

# Read line by line
while line:
    line = f.readline()
    #Detect end of file
    if "$EndFile" in line:
            print("End of linear mesh file!")
            line = False
    else:
        # Read nodes section
        if "Nodes" in line:
            if "NumNodes" in line:
                inNODE = False
                line = True
            elif "EndNodes" in line:
                print("Read linear mesh Nodes!")
                inNODE = False
                line = True
            else:
                inNODE = True
                print("Reading Nodes...")
                line = True
        # Read Elements section
        elif "Elements" in line:
            if "NumElements" in line:
                inELEM = False
                line = True
            elif "EndElements" in line:
                print("Read linear mesh Elements!")
                inELEM = False
                line = True
            else:
                inELEM = True
                print("Reading Elements...")
                line = True
        # Read Physical names section
        elif "PhysicalNames" in line:
            if "EndPhysicalNames" in line:
                print("Read physical groups!")
                inPhNa = False
                line = True
            else:
                inPhNa = True
                print("Reading physical groups...")
                line = True
        # Start operations and data storage
        else:
            # Nodal info storage
            if inNODE:
                stripline = line.strip().split()
                s = len(stripline)
                if s == 1:
                    numNodes = int(stripline[0])
                    linearCoor = np.zeros((numNodes,3))
                    print("Nodes on linear mesh := ",numNodes)
                    i = 0
                elif s == 3:
                    linearCoor[i,0] = stripline[1] # x
                    linearCoor[i,1] = stripline[2] # y
                    numDim = 2
                    i = i+1
                elif s == 4:
                    linearCoor[i,0] = stripline[1] # x
                    linearCoor[i,1] = stripline[2] # y
                    linearCoor[i,2] = stripline[3] # z
                    numDim = 3
                    i = i+1
            # Elemental info storage
            elif inELEM:
                stripline = line.strip().split()
                s = len(stripline)
                if s == 1:
                    numElem = int(stripline[0])
                    linearMeshTable = np.zeros((numElem,28))
                    linearElemTypes = np.zeros(numElem)
                    print("All elements on mesh := ",numElem)
                    i = 0
                else:
                    linearElemTypes[i] = stripline[1] # Linear element type (gmsh)
                    linearMeshTable[i,0] = stripline[3] # Physical group
                    linearMeshTable[i,1:s-4] = stripline[5:s] # Nodal indexes
                    i = i+1
            # Groups storage
            elif inPhNa:
                stripline = line.strip().split()
                s = len(stripline)
                if s == 1:
                    numBoun = int(stripline[0])
                    linearBoun = np.zeros((numBoun,2))
                    print("Total groups := ",numBoun)
                    print("Total boundaries := ",numBoun-1)
                    i = 0
                else:
                    linearBoun[i,0] = stripline[0] # Dimension
                    linearBoun[i,1] = stripline[1] # Group ID
                    i = i+1
f.close()

# Remove boundary elements from mesh table

linearMeshTable.astype(int)
linearBoun.astype(int)
counter = 0
for i in range(numElem):
    if linearMeshTable[i,0] == linearBoun[0,1]:
        counter = counter+1
numIntElem = counter
print("Total interior elements := ",counter)
elements = copy.deepcopy(linearMeshTable[0:counter,1:])
elements.astype(int)
boundaryGroups = copy.deepcopy(linearMeshTable[counter:,0])
boundaries = copy.deepcopy(linearMeshTable[counter:,1:])
numBounElem = boundaries.shape[0]
print("Number of elements on boundaries is = ",numBounElem)

# Read high-order file if it exists

if len(sys.argv) == 4:
    HoMeshName = sys.argv[3]
    HoMeshFile = '{}.txt'.format(HoMeshName)
    f = open(HoMeshFile,'r')
    line = True
    inNODE = False
    inELEM = False
    inELNO = False
    # Start reading additional file
    while line:
        line = f.readline()
        #Detect end of file
        if "$EndFile" in line:
            print("End of HO mesh file!")
            line = False
        else:
            # Read nodes section
            if "Coord" in line:
                if "EndCoord" in line:
                    print("Read HO Nodes!")
                    inNODE = False
                    line = True
                else:
                    inNODE = True
                    print("Reading HO Nodes...")
                    line = True
            # Read Elements section
            elif "Elements" in line:
                if "EndElements" in line:
                    print("Read HO mesh table!")
                    inELEM = False
                    line = True
                else:
                    inELEM = True
                    print("Reading HO mesh table...")
                    line = True
            # Read ElementsNodes section
            elif "ElemNodes" in line:
                if "EndElemNodes" in line:
                    print("Read HO elements/nodes total!")
                    inELNO = False
                    line = True
                else:
                    inELNO = True
                    print("Reading HO elements/nodes total...")
                    line = True    
            # Start operations and data storage
            else:
                # HO elements/nodes info storage
                if inELNO:
                    stripline = line.strip().split()
                    numNodesHO = int(stripline[0])
                    print("Number of HO nodes on mesh := ",numNodesHO)
                    numElemHO = int(stripline[1])
                    hoCoor = np.zeros((numNodesHO,3))
                    hoMeshTable = np.zeros((numElemHO,27))
                    i = 0
                    j = 0
                # Nodal info storage
                elif inNODE:
                    stripline = line.strip().split()
                    s = len(stripline)
                    if s == 2:
                        hoCoor[i,0] = stripline[0] # x
                        hoCoor[i,1] = stripline[1] # y
                        numDimHO = 2
                        i = i+1
                    elif s == 3:
                        hoCoor[i,0] = stripline[0] # x
                        hoCoor[i,1] = stripline[1] # y
                        hoCoor[i,2] = stripline[2] # z
                        numDimHO = 3
                        i = i+1
                # HOMT storage
                elif inELEM:
                    stripline = line.strip().split()
                    s = len(stripline)
                    hoMeshTable[j,0:s] = stripline[0:s] # Nodal indexes
                    j = j+1
    f.close()

    # Start reordering interior connectivity
    # Nodes are organized as [corners|edges|faces|volumes]

    newMeshTable = np.zeros((numElemHO,27))
    for ielem in range(numElemHO):
        # Interior TRI element type
        if linearElemTypes[ielem] == 2:
            newMeshTable[ielem,0] = hoMeshTable[ielem,0]
            elements[ielem,0] = newMeshTable[ielem,0]
            newMeshTable[ielem,1] = hoMeshTable[ielem,2]
            elements[ielem,1] = newMeshTable[ielem,1]
            newMeshTable[ielem,2] = hoMeshTable[ielem,5]
            elements[ielem,2] = newMeshTable[ielem,2]
            newMeshTable[ielem,3] = hoMeshTable[ielem,1]
            newMeshTable[ielem,4] = hoMeshTable[ielem,4]
            newMeshTable[ielem,5] = hoMeshTable[ielem,3]
        # Interior QUA element type
        elif linearElemTypes[ielem] == 3:
            newMeshTable[ielem,0] = hoMeshTable[ielem,0]
            elements[ielem,0] = newMeshTable[ielem,0]
            newMeshTable[ielem,1] = hoMeshTable[ielem,2]
            elements[ielem,1] = newMeshTable[ielem,1]
            newMeshTable[ielem,2] = hoMeshTable[ielem,8]
            elements[ielem,2] = newMeshTable[ielem,2]
            newMeshTable[ielem,3] = hoMeshTable[ielem,6]
            elements[ielem,3] = newMeshTable[ielem,3]
            newMeshTable[ielem,4] = hoMeshTable[ielem,1]
            newMeshTable[ielem,5] = hoMeshTable[ielem,5]
            newMeshTable[ielem,6] = hoMeshTable[ielem,7]
            newMeshTable[ielem,7] = hoMeshTable[ielem,3]
            newMeshTable[ielem,8] = hoMeshTable[ielem,4]
        # Interior TET element type
        elif linearElemTypes[ielem] == 4:
            newMeshTable[ielem,0] = hoMeshTable[ielem,0]
            elements[ielem,0] = newMeshTable[ielem,0]
            newMeshTable[ielem,1] = hoMeshTable[ielem,2]
            elements[ielem,1] = newMeshTable[ielem,1]
            newMeshTable[ielem,2] = hoMeshTable[ielem,5]
            elements[ielem,2] = newMeshTable[ielem,2]
            newMeshTable[ielem,3] = hoMeshTable[ielem,9]
            elements[ielem,3] = newMeshTable[ielem,3]
            newMeshTable[ielem,4] = hoMeshTable[ielem,1]
            newMeshTable[ielem,5] = hoMeshTable[ielem,4]
            newMeshTable[ielem,6] = hoMeshTable[ielem,3]
            newMeshTable[ielem,7] = hoMeshTable[ielem,6]
            newMeshTable[ielem,8] = hoMeshTable[ielem,7]
            newMeshTable[ielem,9] = hoMeshTable[ielem,8]
        # Interior HEX element type
        elif linearElemTypes[ielem] == 5:
            print()
            newMeshTable[ielem,0] = hoMeshTable[ielem,0]
            elements[ielem,0] = newMeshTable[ielem,0]
            newMeshTable[ielem,1] = hoMeshTable[ielem,2]
            elements[ielem,1] = newMeshTable[ielem,1]
            newMeshTable[ielem,2] = hoMeshTable[ielem,5]
            elements[ielem,2] = newMeshTable[ielem,2]
            newMeshTable[ielem,3] = hoMeshTable[ielem,9]
            elements[ielem,3] = newMeshTable[ielem,3]
            newMeshTable[ielem,4] = hoMeshTable[ielem,1]
            elements[ielem,4] = newMeshTable[ielem,4]
            newMeshTable[ielem,5] = hoMeshTable[ielem,4]
            elements[ielem,5] = newMeshTable[ielem,5]
            newMeshTable[ielem,6] = hoMeshTable[ielem,3]
            elements[ielem,6] = newMeshTable[ielem,6]
            newMeshTable[ielem,7] = hoMeshTable[ielem,6]
            elements[ielem,7] = newMeshTable[ielem,7]
            newMeshTable[ielem,8] = hoMeshTable[ielem,7]
            newMeshTable[ielem,9] = hoMeshTable[ielem,8]
            newMeshTable[ielem,10] = hoMeshTable[ielem,1]
            newMeshTable[ielem,11] = hoMeshTable[ielem,1]
            newMeshTable[ielem,12] = hoMeshTable[ielem,4]
            newMeshTable[ielem,13] = hoMeshTable[ielem,3]
            newMeshTable[ielem,14] = hoMeshTable[ielem,6]
            newMeshTable[ielem,15] = hoMeshTable[ielem,7]
            newMeshTable[ielem,16] = hoMeshTable[ielem,8]
            newMeshTable[ielem,17] = hoMeshTable[ielem,4]
            newMeshTable[ielem,18] = hoMeshTable[ielem,3]
            newMeshTable[ielem,19] = hoMeshTable[ielem,6]
            newMeshTable[ielem,20] = hoMeshTable[ielem,7]
            newMeshTable[ielem,21] = hoMeshTable[ielem,8]
            newMeshTable[ielem,22] = hoMeshTable[ielem,8]
            newMeshTable[ielem,23] = hoMeshTable[ielem,4]
            newMeshTable[ielem,24] = hoMeshTable[ielem,3]
            newMeshTable[ielem,25] = hoMeshTable[ielem,6]
            newMeshTable[ielem,26] = hoMeshTable[ielem,7]
            print("Do not trust this yet!!!!")
        else:
            print("Element type not yet defined!!!")
            break

    # Start reordering boundary nodes
    # Check which boundary belongs to which element

    numBounElem = boundaries.shape[0]
    pairing = np.zeros((numBounElem,2))
    newBoundaries = np.zeros((numBounElem,27))
    for iboun in range(numBounElem):
        bounElem = []
        newBounElem = []
        for inodb in range(27):
            if boundaries[iboun,inodb] == 0:
                sizeBoun = len(bounElem)
                break
            else:
                bounElem.append(boundaries[iboun,inodb].astype(int))
        for ielem in range(numElemHO):
            oldElement = linearMeshTable[ielem,1:]
            counter = 0
            match = False
            for inodb in bounElem:
                if inodb in oldElement:
                    counter += 1
                if counter == sizeBoun:
                    match = True
                    break
            if match:
                pairing[iboun,0] = iboun
                pairing[iboun,1] = ielem
                break
        indexElem = pairing[iboun,1]
        oldElement = linearMeshTable[ielem,1:]
        newElement = newMeshTable[ielem,:]
        # Renumber the original boundary
        bounNodeIndex = []
        for inodb in bounElem:
            if inodb in oldElement:
                indexNodes = oldElement.tolist().index(inodb)
                bounNodeIndex.append(indexNodes)
                newNodeNum = newElement[indexNodes].astype(int)
                newBounElem.append(newNodeNum)
        newBoundaries[iboun,0:sizeBoun] = newBounElem
        # Add new nodes
        if linearElemTypes[ielem] == 4: # TET
            if 0 in bounNodeIndex and 1 in bounNodeIndex and 2 in bounNodeIndex:
                bounNodeIndex.append(4)
                bounNodeIndex.append(5)
                bounNodeIndex.append(6)
            elif 0 in bounNodeIndex and 2 in bounNodeIndex and 3 in bounNodeIndex:
                bounNodeIndex.append(6)
                bounNodeIndex.append(9)
                bounNodeIndex.append(7)
            elif 0 in bounNodeIndex and 1 in bounNodeIndex and 3 in bounNodeIndex:
                bounNodeIndex.append(7)
                bounNodeIndex.append(8)
                bounNodeIndex.append(4)
            elif 1 in bounNodeIndex and 2 in bounNodeIndex and 3 in bounNodeIndex:
                bounNodeIndex.append(5)
                bounNodeIndex.append(9)
                bounNodeIndex.append(8)
            newSizeBoun = len(bounNodeIndex)
            newBoundaries[iboun,3:newSizeBoun] = newMeshTable[ielem,bounNodeIndex[3:]]
        else:
            print("Element not coded!")

    # Start writing files

    hoDimsFile = 'ho_{}.dims.dat'.format(caseName)
    hoGeoDatFile = 'ho_{}.geo.dat'.format(caseName)
    hoFixBouFile = 'ho_{}.fix.bou'.format(caseName)

    # *.dims.dat

    f=open(hoDimsFile,'w')
    f.write(' '.join(map(str,["NODAL_POINTS ",numNodesHO]))+'\n')
    f.write(' '.join(map(str,["ELEMENTS     ",numElemHO]))+'\n')
    f.write(' '.join(map(str,["BOUNDARIES   ",numBounElem])))
    f.close()

    # *.fix.bou

    f=open(hoFixBouFile,'w')
    for iboun in range(numBounElem):
        if iboun == numBounElem-1:
            f.write(' '.join(map(str,[iboun+1,linearMeshTable[numElemHO+iboun,0].astype(int)])))
        else:
            f.write(' '.join(map(str,[iboun+1,linearMeshTable[numElemHO+iboun,0].astype(int)]))+'\n')
    f.close()

    # *.geo.dat

    f=open(hoGeoDatFile,'w')
    line = True
    while line:
        f.write("NODES_PER_ELEMENT"+'\n')
        for ielem in range(numElemHO):
            counter = 0
            for inode in newMeshTable[ielem,:].astype(int):
                if inode == 0:
                    break
                else:
                    counter += 1
            f.write(' '.join(map(str,[ielem+1,counter]))+'\n')
        f.write("END_NODES_PER_ELEMENT"+'\n')
        f.write("ELEMENTS"+'\n')
        for ielem in range(numElemHO):
            counter = 0
            for inode in newMeshTable[ielem,:].astype(int):
                if inode == 0:
                    break
                else:
                    counter += 1
            element = newMeshTable[ielem,0:counter].astype(int)
            element = np.insert(element,0,ielem+1)
            string = str(element).strip('[]')
            f.write(string+'\n')
        f.write("END_ELEMENTS"+'\n')
        f.write("COORDINATES"+'\n')
        for inode in range(numNodesHO):
            xyz = hoCoor[inode,0:numDim]
            xyz = np.insert(xyz,0,(inode+1))
            xyz = xyz.tolist()
            xyz[0] = int(xyz[0])
            if numDim == 2:
                f.write(str(xyz[0])+' '+str(xyz[1])+' '+str(xyz[2])+'\n')
            elif numDim == 3:
                f.write(str(xyz[0])+' '+str(xyz[1])+' '+str(xyz[2])+' '+str(xyz[3])+'\n')
        f.write("END_COORDINATES"+'\n')
        f.write("BOUNDARIES"+'\n')
        for iboun in range(numBounElem):
            counter = 0
            for inodb in newBoundaries[iboun,:].astype(int):
                if inodb == 0:
                    break
                else:
                    counter += 1
            boundary = newBoundaries[iboun,0:counter].astype(int)
            boundary = np.insert(boundary,0,iboun+1)
            string = str(boundary).strip('[]')
            f.write(string+'\n')
        f.write("END_BOUNDARIES"+'\n')
        f.write("SKEW_SYSTEMS"+'\n')
        f.write("END_SKEW_SYSTEMS")
        line = False
    f.close()

    print("High-order files written!")

# Write linear files

DimsFile = '{}.dims.dat'.format(caseName)
GeoDatFile = '{}.geo.dat'.format(caseName)
FixBouFile = '{}.fix.bou'.format(caseName)

# *.dims.dat

f=open(DimsFile,'w')
f.write(' '.join(map(str,["NODAL_POINTS ",numNodes]))+'\n')
f.write(' '.join(map(str,["ELEMENTS     ",numIntElem]))+'\n')
f.write(' '.join(map(str,["BOUNDARIES   ",numBounElem])))
f.close()

# *.fix.bou

f=open(FixBouFile,'w')
for iboun in range(numBounElem):
    if iboun == numBounElem-1:
        f.write(' '.join(map(str,[iboun+1,linearMeshTable[numIntElem+iboun,0].astype(int)])))
    else:
        f.write(' '.join(map(str,[iboun+1,linearMeshTable[numIntElem+iboun,0].astype(int)]))+'\n')
f.close()

# *.geo.dat

f=open(GeoDatFile,'w')
line = True
while line:
    f.write("NODES_PER_ELEMENT"+'\n')
    for ielem in range(numIntElem):
        counter = 0
        for inode in linearMeshTable[ielem,1:].astype(int):
            if inode == 0:
                break
            else:
                counter += 1
        f.write(' '.join(map(str,[ielem+1,counter]))+'\n')
    f.write("END_NODES_PER_ELEMENT"+'\n')
    f.write("ELEMENTS"+'\n')
    for ielem in range(numIntElem):
        counter = 0
        for inode in linearMeshTable[ielem,1:].astype(int):
            if inode == 0:
                break
            else:
                counter += 1
        f.write(str(ielem+1))
        for inode in range(1,counter+1):
            f.write(" ")
            string = str(linearMeshTable[ielem,inode].astype(int))
            f.write(string)
        #element = linearMeshTable[ielem,1:counter+1].astype(int)
        #element = np.insert(element,0,ielem+1)
        #string = str(element).strip('[]')
        #f.write(string+'\n')
        f.write('\n')
    f.write("END_ELEMENTS"+'\n')
    f.write("COORDINATES"+'\n')
    for inode in range(numNodes):
        xyz = linearCoor[inode,0:numDim]
        xyz = np.insert(xyz,0,(inode+1))
        xyz = xyz.tolist()
        xyz[0] = int(xyz[0])
        if numDim == 2:
            f.write(str(xyz[0])+' '+str(xyz[1])+' '+str(xyz[2])+'\n')
        elif numDim == 3:
            f.write(str(xyz[0])+' '+str(xyz[1])+' '+str(xyz[2])+' '+str(xyz[3])+'\n')
    f.write("END_COORDINATES"+'\n')
    f.write("BOUNDARIES"+'\n')
    for iboun in range(numBounElem):
        counter = 0
        for inodb in linearMeshTable[numIntElem+iboun,1:].astype(int):
            if inodb == 0:
                break
            else:
                counter += 1
        boundary = linearMeshTable[numIntElem+iboun,1:counter+1].astype(int)
        boundary = np.insert(boundary,0,iboun+1)
        string = str(boundary).strip('[]')
        f.write(string+'\n')
    f.write("END_BOUNDARIES"+'\n')
    f.write("SKEW_SYSTEMS"+'\n')
    f.write("END_SKEW_SYSTEMS")
    line = False
f.close()

print("Linear mesh files written!")
print("End of program!")

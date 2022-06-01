# Post-process Alya node set file data
# Created by G.Guillamet <gerard.guillamet@bsc.es>
# Date: 6 October 2017

from operator import itemgetter

# Inputs
namemod = 'block'
nameinp = namemod + '-node.sld.set'
nameout = 'alya-solution-path-block.txt'
lines_h = 13    # Time (Line 13, position 12)

# Initialize and get number of nodes from node-set
with open(namemod + '.set.dat', 'r') as infile:
    numline = 0
    for line in infile:
        numline += 1
nnodes = numline - 2

# Reading
timeList,lineList = [],[]
with open(nameinp, 'r') as file:
    lc = 0
    for (i, line) in enumerate(file):
            lc += 1
            linecurr = line.split()
            if lc >= lines_h:
                    if linecurr[1] == 'Time':
                        timeList.append(float(linecurr[3]))
                        lineList.append(lc)

# Times and lines for postprocess
firs_time_List = timeList[0]
firs_line_List = lineList[0]
last_time_List = timeList[-1]
last_line_List = lineList[-1]

# Save data
coordxList, coordyList = [], []
with open(nameinp, 'r') as file:
    lc = 0
    for line in file:
        lc += 1
        if lc >= firs_line_List+1 and lc < firs_line_List+nnodes+1:
            linecurr = line.split()
            coordx = float(linecurr[3])
            coordy = float(linecurr[4])
            coordxList.append(coordx)
            coordyList.append(coordy)

displxList, displyList = [], []
with open(nameinp, 'r') as file:
    lc = 0
    for line in file:
        lc += 1
        if lc >= last_line_List+1 and lc < last_line_List+nnodes+1:
            linecurr = line.split()
            displx = float(linecurr[1])
            disply = float(linecurr[2])
            displxList.append(displx)
            displyList.append(disply)

# Sort Nodes according to their coordinates
data = []
for i in range(len(coordxList)):
    data.append([coordxList[i], coordyList[i], displxList[i], displyList[i]])
sorted_data = sorted(data,key=itemgetter(0))

lineouts = []
for i in range(len(coordxList)):
    lineouts.append('%1.6e %1.6e %1.6e %1.6e\n' % (sorted_data[i][0],sorted_data[i][1],sorted_data[i][2],sorted_data[i][3]))

# Write file output
with open(nameout, 'w') as file:
        file.writelines(lineouts)


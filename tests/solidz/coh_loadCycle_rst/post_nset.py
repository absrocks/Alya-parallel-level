# Post-process Alya node set file data
# Created by G.Guillamet <gerard.guillamet@bsc.es>
# Date: 6 September 2017
# Details: Posprocess only the last 3 steps

# Inputs
namemod = 'coh_loadCycle_rst_i4'
nameinp = namemod + '-node.sld.set'
nameout = 'alya-solution2.txt'
lines_h = 11    # Time (Line 11, position 10)

# Initialize and get number of nodes from node-set
with open(namemod + '.dom.set.dat', 'r') as infile:
    numline = 0
    for line in infile:
        numline += 1
nnodes = numline - 2
lineouts = []
timeList = []
lineList = []

# Reading
with open(nameinp, 'r') as file:
    lc = 0
    for (i, line) in enumerate(file):
        lc += 1
        linecurr = line.split()
        if lc >= lines_h:
            if linecurr[1] == 'Time':
                timeList.append(float(linecurr[3]))
                lineList.append(lc)

# Save data
for istep in range(len(timeList)):
    if istep >= 4:
        lineno = lineList[istep]
        lc, rforce, displa = 0, 0, 0
        with open(nameinp, 'r') as file:
            for line in file:
                lc += 1
                if lc >= lineno+1 and lc < lineno+nnodes+1:
                    linecurr = line.split()
                    displa += float(linecurr[1])
                    rforce += float(linecurr[2])
        lineouts.append('%1.6e %1.6e %1.6e\n' % (timeList[istep],displa/nnodes,rforce))

# Write file output
with open(nameout, 'w') as file:
    file.writelines(lineouts)


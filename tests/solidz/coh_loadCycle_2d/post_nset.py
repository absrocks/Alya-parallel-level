# Post-process Alya node set file data
# Created by G.Guillamet <gerard.guillamet@bsc.es>

# Inputs
nameinp = 'coh_loadCycle_2d-node.sld.set'
nameout = 'post-alyasol.txt'
nnodes  = 2
rf      = 'FRXIY'
u       = 'DISPY'

# Initialize
lineouts = []
timeList = []
rforList = []
dispList = []
lineList = []

# Reading
with open(nameinp, 'r') as file:
    lc = 0
    for (i, line) in enumerate(file):
            lc += 1
            linecurr = line.split()
            if lc >= 8:
                    if linecurr[1] == 'Time':
                            timeList.append(float(linecurr[3]))
                            lineList.append(i)

for istep in range(len(timeList)):
    lineno = lineList[istep]+1
    lc, rforce, displa = 0, 0, 0
    with open(nameinp, 'r') as file:
        for line in file:
            lc += 1
            if lc >= lineno+1 and lc < lineno+nnodes+1:
                linecurr = line.split()
                displa += float(linecurr[1])
                rforce += float(linecurr[2])
        rforList.append(rforce)
        dispList.append(displa/nnodes)
        if rf == 'FRXIX' and u == 'DISPX':
           lineouts.append('%1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e\n' % (timeList[istep],displa/nnodes,0.0,0.0,rforce,0.0,0.0))
        if rf == 'FRXIY' and u == 'DISPY':
           lineouts.append('%1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e\n' % (timeList[istep],0.0,displa/nnodes,0.0,0.0,rforce,0.0))
        if rf == 'FRXIZ' and u == 'DISPZ':
           lineouts.append('%1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e\n' % (timeList[istep],0.0,0.0,displa/nnodes,0.0,0.0,rforce))

# Write file output
with open(nameout, 'w') as file:
        file.writelines(lineouts)


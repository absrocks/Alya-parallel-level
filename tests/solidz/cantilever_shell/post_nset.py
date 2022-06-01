# Post-process Alya node set file data
# Created by G.Guillamet <gerard.guillamet@bsc.es>

# Inputs
nameinp = 'cantilever_shell-node.sld.set'
nameout = 'alya-solution.txt'
nnodes  = 6
qmax    = 4.0
l       = 1.0
lines_h = 12       # Line 12, position 11

# Initialize
lineouts = []
timeList = []
dispList = []
lineList = []

# Read step times and positions
with open(nameinp, 'r') as file:
    lc = 0
    for (i, line) in enumerate(file):
            lc += 1
            linecurr = line.split()
            if lc >= lines_h:
                    if linecurr[1] == 'Time':
                        timeList.append(float(linecurr[3]))
                        lineList.append(lc)
#
# Save data
#
for istep in range(len(timeList)):
    lineno = lineList[istep]
    lc, displa_z, displa_y, rforce = 0, 0, 0, 0
    with open(nameinp, 'r') as file:
        for line in file:
            lc += 1
            if lc >= lineno+1 and lc < lineno+nnodes+1:
                linecurr = line.split()
                displa_y += float(linecurr[1])
                displa_z += float(linecurr[2])
                rforce   += float(linecurr[3])

        lineouts.append('%1.6e %1.6e %1.6e %1.6e\n' % (timeList[istep],abs(displa_y/2),abs(displa_z/2),abs(rforce/qmax/l)))

# Write file output
with open(nameout, 'w') as file:
        file.writelines(lineouts)


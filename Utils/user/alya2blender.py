#!/usr/bin/python
import glob
import array
import sys
import math
from optparse import OptionParser

infinity = 1.0e20

#
# Read argument list: 
#
usage = "Join individual voxel frames into a Blender animation.\n     syntax: %prog [options] AlyaJobName"
parser = OptionParser(usage=usage)
parser.add_option("-l","--local",action="store_false", 
	dest="GlobalScale",default=True, 
	help="Adjust zero density to the minimum of each frame, default is adjust to minimum among all frames.")
parser.add_option("-v","--variable",action="store", type="string",   
	dest="Variable",help="Select which variable to postprocess.")
parser.add_option("-e","--endian",action="store_true", 
	dest="SwitchEndianness",default=False, 
	help="Force switching the endianness.")
parser.add_option("-n","--noendian",action="store_false", 
	dest="DetectEndianness",default=True, 
	help="Switch off detection of endianness (default is yes).")
parser.add_option("-s","--rescale",action="store_true", 
	dest="ReScale",default=False, 
	help="Rescale the data to the range [0,1].")
parser.add_option("-V","--verbose",action="store_true", 
	dest="Verbose",default=False, 
	help="Switch on verbose mode.")
parser.add_option("-x","--stats",action="store_true", 
	dest="Statistics",default=False, 
	help="Return some statistics on the data.")
parser.add_option("-m", "--min", type="float", action="store", dest="FixedMinimum", default = infinity,
                  help="Fixes the minimum value below which density is zero.")
parser.add_option("-M", "--max", type="float", action="store", dest="FixedMaximum", default = -infinity,
                  help="Fixes the Maximum valuefor the density.")

# Parse arguments
(options, args) = parser.parse_args()

jobname = args[0]
GlobalScale = options.GlobalScale
SwitchEndianness = options.SwitchEndianness
DetectEndianness = options.DetectEndianness
ReScale = options.ReScale
FixedMinimum = options.FixedMinimum
FixedMaximum = options.FixedMaximum
Verbose = options.Verbose
Statistics = options.Statistics
Variable = options.Variable

# Decide on some value
if Verbose: 
    print ("Fixed minimum initially set to "+str(FixedMinimum))
    print ("Fixed maximum initially set to "+str(FixedMaximum))

# Some preparations
dims = array.array('i')
singledims = array.array('i')
data = array.array('f')
minmax = array.array('f')

#
# Pick up which files to process
#
TotalFileList = glob.glob(jobname+"*.bvox")

#
# Try to detect endianness 
# 
if (DetectEndianness):
    print "\nDetecting endianness of raw data..."
    with open(TotalFileList[0],"rb") as f:
	try: singledims.fromfile(f, 4)
	except EOFError: print "Error detecting the endianness.\n"
	if (singledims[3]!=1): # Frames should be 1
	    print "...switch endianness activated"
	    SwitchEndianness=True
	else:
	    SwitchEndianness=False

#
# TODO :: CAREFUL BECAUSE GLOBAL SCALES WILL GO ACROSS VARIABLES, THIS MUST BE CORRECTED

#
# Detect which variables are present
#
jobNameLength = len(jobname)
PostProcVars = set()
if (len(Variable)>0):
    PostProcVars.add(Variable[0:4])
else:
    for name in TotalFileList:
        PostProcVars.add(name[jobNameLength+1:jobNameLength+6])

#
# if global is present, equalize scale among different files
#
if (GlobalScale):  
    #Read header and tail in all files, but only pick up max and min
    for variable in PostProcVars:
	fileList = glob.glob(jobname+"-"+variable+"*.bvox")
	for filename in fileList:
	    with open(filename,"rb") as f:
		singledims = array.array('i')
		try: singledims.fromfile(f, 4)
		except EOFError: print "\nGlobal scale: Error reading dims"
		if SwitchEndianness:
		    singledims.byteswap()
		f.seek(16+4*singledims[0]*singledims[1]*singledims[2]*singledims[3]) # We skip the data
	        try: minmax.fromfile(f, 2)
		except EOFError: print "\nGlobal scale: Error reading min max"
    if SwitchEndianness:
	minmax.byteswap()
    if (FixedMinimum < infinity):
	globalMin = FixedMinimum
    else:
	globalMin = min (minmax)
    if (FixedMaximum > -infinity):
	globalMax = FixedMaximum
    else:
	globalMax = max (minmax)
    GlobalScale = globalMax - globalMin
    dims = singledims # We pick up the dims from the last file read
else:
    globalMin = infinity
    globalMax = -infinity
    # We have to read one header to pick up dims
    with open(TotalFileList[0],"rb") as f:
	singledims = array.array('i')
	try: singledims.fromfile(f, 4)
	except EOFError: print "\nPreparation: Error reading dims"
    if SwitchEndianness:
	singledims.byteswap()
	dims = singledims

if Verbose:
    print "\nFound voxel resolution: "
    print "    X = "+str(dims[0])
    print "    Y = "+str(dims[1])
    print "    Z = "+str(dims[2])
    print "\nMin detected: "+str(globalMin)
    print "\nMax detected: "+str(globalMax)

#
# Now write files together
#
zero = array.array('f', [0.0])
one = array.array('f', [1.0])
shifteddata = array.array('f' , [0.0] )
for variable in PostProcVars:
    print "\nPostprocessing "+variable
    fileList = glob.glob(jobname+"-"+variable+"*.bvox")
    with open(jobname+"-"+variable+".bvox.anim","wb") as out:
	# We write first the header, although in principle this should be a loop like alya2pos.x
	dims[3] = len(fileList) #number of frames
	dims.tofile(out)
	for filename in sorted(fileList):
	    print "...frame "+filename[jobNameLength+7:-5]
	    # we must reset the temporary variables or they stack up and we cannot correct the min and max values
            singledims = array.array('i')
            data = array.array('f')
            minmax = array.array('f')
            with open(filename,"rb") as f:
                # Read dims
		try: singledims.fromfile(f, 4)
		except EOFError: print "Process: Error reading dims"
		if SwitchEndianness:
		    singledims.byteswap()
		dims = singledims
                # Read data
		try: data.fromfile(f, dims[0]*dims[1]*dims[2]*dims[3])
		except EOFError: print "Process: Error reading data"+str(dims[0]*dims[1]*dims[2]*dims[3])+" bytes"
		if SwitchEndianness:
		    data.byteswap()
                # Read max min
   	        try: minmax.fromfile(f, 2)
		except EOFError: print "Error reading min max"
		if SwitchEndianness:
                    minmax.byteswap()
                # Adjust local max min
              	if (GlobalScale):
		    theMinimum = globalMin
		    theMaximum = globalMax
		else:
                    theMinimum = minmax[0]
                    theMaximum = minmax[1]
		if (FixedMaximum > -infinity):
		    theMaximum = FixedMaximum
		if (FixedMinimum < infinity):
		    theMinimum = FixedMinimum
	        ScaleFactor = theMaximum-theMinimum
		if (Statistics):
		    n = len(data)
		    mean = sum(data) / n
		    sd = math.sqrt(sum((x-mean)**2 for x in data) / n)
		    print ("Mean is "+str(mean))
		    print ("STD is "+str(sd))
		    newdata = array.array('f')
		    bins = 100
		    histogram = {}
		    for x in range(bins):
			histogram[x] = 0
                # Finally...
		if Verbose:
		    print ("Writing")
		    print ("Min is "+str(theMinimum))
		    print ("Max is "+str(theMaximum))
		    print ("Scale is "+str(ScaleFactor))
		for item in data:
                    # Shift data according to min/max
		    shifteddata[0] =  (item+minmax[0]-theMinimum)/ScaleFactor
		    if (shifteddata[0] < theMaximum/ScaleFactor):
                        if (shifteddata[0] >=0.0):
		    	    if Statistics: 
			        newdata.append(shifteddata[0])
		                bin= int(math.floor(bins*shifteddata[0]))
			        if bin in histogram:
			            histogram[bin] += 1
                            shifteddata.tofile(out)
                        else:
                            zero.tofile(out)
                    else:
                        one[0] = theMaximum/ScaleFactor
                        one.tofile(out)
		if (Statistics):
		    n = len(newdata)
		    mean = sum(newdata) / n
		    sd = math.sqrt(sum((x-mean)**2 for x in newdata) / n)
		    print ("New Mean is "+str(mean))
		    print ("New STD is "+str(sd))
		    for key,val in histogram.items():
			print(str(key*1.0/bins)+"  "+str(key*1.0*ScaleFactor/bins+theMinimum-minmax[0])+"  "+str(histogram[key]))

print "Detected "+str( len(fileList))+" frames.\n"
print "Done\n"

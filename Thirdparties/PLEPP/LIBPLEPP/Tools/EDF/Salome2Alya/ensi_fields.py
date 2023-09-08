#!/usr/bin/python

import numpy as np
import sys
import os
import glob

#--------------------------------------------------------------------------------------------# 
def KerDat():
  fname = glob.glob("*.ker.dat")
  if(len(fname)<1):
    print "WARNING: there no exist *.nsa.wit files! \n\n"
    pass
  else:
    fname = fname[0]
    print "|_\'%s\'" % fname

    data = open(fname, "r")
    lines = data.readlines()
    data.close()

    ok = False
    for line in lines:
      if(line.find("WITNESS_POINTS")>0):     ok = True
      if(line.find("END_WITNESS_POINTS")>0): ok = False
      if(ok): print line[:-1]







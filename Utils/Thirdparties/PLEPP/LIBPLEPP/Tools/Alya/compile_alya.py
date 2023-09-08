import sys
import glob

ROOT_ALYA    = sys.argv[1] 
ROOT_COMMDOM = sys.argv[2] 

print "|_[Ple++] Directory: \'%s\' " % ROOT_COMMDOM
print "|_[Alya]  Directory: \'%s\' " % ROOT_ALYA


fname = ROOT_ALYA + "/Executables/unix" 

files = glob.glob(fname+"/conf.log") 
if(len(files)==0):
  print "\n\n   ERROR: there is not \'%s\' \n\n" % (fname+"/conf.log")
  sys.exit(1)


fin = open(fname+"/conf.log", "r") 
lines = fin.readlines() 
fin.close() 


nline = len(lines)
for i in range(nline):  
  line = lines[i]
  if( line.find('./configure') != -1 and line.find('-f') != -1 ):
    new_line = line[:-1]
    if(new_line.find("_ple") != -1):
      sys.exit() 

    split_lines = line[:].split(" ") 
    for split in split_lines: 
      if( split.find("-f")!= -1 ): 
        new_line = split.split("=") 


configure_in = new_line[-1].split("/")[-1]

fin = open(fname+"/"+new_line[1], "r") 
lines02 = fin.readlines() 
fin.close() 


nline = len(lines02)
for i in range(nline):  
  line = lines02[i]
  if( line.find('libs') != -1 and line.find('==') != -1 ):
    new_line = line[:-1]  + " -L%s -lcommdom \n" % (ROOT_COMMDOM) 
    lines02[i] = new_line 

  #if( line.find('mpif90') != -1 and line.find('==') != -1 ):
  #  new_line = line[:-1]  + " -D __COMMDOM__=$(COMMDOM) \n"
  #  lines02[i] = new_line 



configure_in = configure_in.split(".")
print fname+"/"+configure_in[0]+"_ple.txt"

fout = open(fname+"/"+configure_in[0]+"_ple.txt", "w") 
fout.writelines(lines02) 
fout.close() 


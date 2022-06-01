'''                      ALYA2CSV.PY
        CONVERTS ALYA'S PARTICLES OUTPUT TO A *.CSV FILE 
                        Edgar Olivares                       '''


# Insert and read the name of the case
value = raw_input( "Name of the case: ")
file = value + ".pts.res"

# Creating ouptut file with 9 columns
o = open( value + ".pts.csv", "w")
print '### From ' + file + ", creating a new file called  " + value + ".pts.csv ###"
o.write("time,part_id,xcoord,ycoord,zcoord,xvelo,yvelo,zvelo,part_type,subdom" + "\n")

# Open the file we want to read
f = open(file, "r")

# Writing rows as a csv file format
for l in f:
    if not l.startswith("#"):
        fields = l.split()
        for i in range(10):
            text = str(fields[i])
            if i < 9:
                o.write(text + ",")
            if i == 9:
                o.write(text + "\n")
            
f.close(); o.close()
print '### Done ###'

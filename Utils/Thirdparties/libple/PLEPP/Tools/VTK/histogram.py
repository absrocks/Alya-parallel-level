import numpy as np 


PATH = "CSVs"


from glob import glob
Filenames = glob(PATH+'/histogram01.*.csv')
n_Filenames = len(Filenames)

Fout = open("cosa.dat", "w")

Mins = []
Maxs = []
Mavr = []
for i in range(n_Filenames):
  filename = "CSVs/histogram01.%d.csv" % i 
  print filename, 
  Data = np.loadtxt(filename, delimiter=',', comments = '\"')

  idx      = Data[:,1].argmax() 
  dens_max = Data[idx,0]
  
  idx      = Data[:,1].argmin() 
  dens_min = Data[idx,0]

  dens_min = Data[:,0].min()
  dens_max = Data[:,0].max()
  averg = np.average( Data[:,0] )
  
  Mins.append(dens_min) 
  Maxs.append(dens_max) 
  Mavr.append(averg) 
  
  print i, [dens_min, averg, dens_max]

  print>> Fout, i, dens_min, averg, dens_max

Mins = np.array(Mins)
Maxs = np.array(Maxs)
Mavr = np.array(Mavr)

print np.average(Mins), 
print np.average(Mavr), 
print np.average(Maxs)



Fout.close() 

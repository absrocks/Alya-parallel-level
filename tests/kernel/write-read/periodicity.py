#!/usr/bin/env python
# coding: utf-8

# In[1]:


import vtk
import numpy as np
import matplotlib.pyplot as plt
import io


# In[2]:


path = '../fields-ok/ensi/write-read.ensi.case' 
path_to_steps = './fields/' #where to save steps
ndim = 2 #dimensionality of the problem
ncpus = 1 #on how many cpus to parallelize
binary_ensight = True #are the ensight case files binary?
#the time interval from which to save velocoity fields: t in [t_initial, t_final)
t_initial = 5 #inclusive
t_final = 6 #exclusive

if binary_ensight:
    rdr = vtk.vtkEnSightGoldBinaryReader()
else :
    rdr = vtk.vtkEnSightGoldReader()
    
rdr.SetCaseFileName(path)
rdr.ReadAllVariablesOn()
rdr.Update()
ntimesteps = rdr.GetTimeSets().GetItem(0).GetNumberOfTuples()

timesteps = []
for time_idx in range(ntimesteps):
    time_instant  = rdr.GetTimeSets().GetItem(0).GetTuple1(time_idx)
    timesteps.append(np.round(time_instant,5))

timesteps = np.array(timesteps)
print(timesteps)


# In[3]:



t20 = timesteps[(timesteps>=t_initial) & (timesteps<t_final)]
print('Saving timesteps:',t20)

# In[4]:


rdr.SetTimeValue(timesteps[0])
rdr.Update()
mesh = rdr.GetOutput().GetBlock(0)

npts = mesh.GetNumberOfPoints()
v = np.zeros((mesh.GetNumberOfPoints(), ndim*t20.shape[0]))
#p = np.zeros((mesh.GetNumberOfPoints(), t20.shape[0]))


for i in range(t20.shape[0]):
    t = t20[i]
    rdr.SetTimeValue(t)
    rdr.Update()
    mesh = rdr.GetOutput().GetBlock(0)
    
    
    veloc = mesh.GetPointData().GetArray('VELOC')
    #press = mesh.GetPointData().GetArray('PRESS')

    for j in range(veloc.GetNumberOfTuples()):
        vv = veloc.GetTuple3(j)
        #pp = press.GetTuple1(j)

        #p[j,i] = pp
        v[j,(ndim*i):(ndim*i+ndim)] = vv[0:ndim]


# In[14]:



def store_fields(stepid):
#    with open(f'press.{stepid}.in','w') as fpress:
        with open(f'{path_to_steps}veloc.{stepid}.in','w') as fveloc:
#            fpress.write(f'STEP {stepid+1}\n')
#            fveloc.write(f'STEP {stepid+1}\n')
#            string_press = io.StringIO()
            string_veloc = io.StringIO()

            if ndim==3:
                for pointid in range(npts):
    #                pstring = f'{pointid+1} {p[pointid, stepid]}\n'
                    vstring = f'{pointid+1} {v[pointid,3*stepid]} {v[pointid,3*stepid+1]} {v[pointid,3*stepid+2]}\n'

    #                string_press.write(pstring)
                    string_veloc.write(vstring)
            elif ndim==2:
                for pointid in range(npts):
    #                pstring = f'{pointid+1} {p[pointid, stepid]}\n'
                    vstring = f'{pointid+1} {v[pointid,2*stepid]} {v[pointid,2*stepid+1]}\n'

    #                string_press.write(pstring)
                    string_veloc.write(vstring)
        
        
#            fpress.write(string_press.getvalue())
            fveloc.write(string_veloc.getvalue())
#            fpress.write(f'END_STEP\n')
#            fveloc.write(f'END_STEP\n')


from multiprocessing import Pool
with Pool(processes=ncpus) as pool:
    for i, _ in enumerate(pool.imap_unordered(store_fields, range(t20.shape[0])), 1):
        print('\rdone {0:%}'.format(i/t20.shape[0]))

    


# In[16]:


#save steps structure for the dom.dat
#with open('steps.press.in','w') as f:
#    for stepid in progressbar.progressbar(range(p.shape[1])):
#        f.write(f'STEP {stepid+1}\n')
#        f.write(f'   INCLUDE {path_to_steps}/press.{stepid}.in\n')
#        f.write(f'END_STEP\n')

with open('steps.veloc.in','w') as f:
    for stepid in range(t20.shape[0]):
        f.write(f'STEP {stepid+1}\n')
        f.write(f'   INCLUDE {path_to_steps}veloc.{stepid}.in\n')
        f.write(f'END_STEP\n')
        


# In[ ]:


#save times
with open('times.in','w') as ftimes:
    for t in t20:
        ftimes.write(f'{np.round(t-t20[0],8)}\n')


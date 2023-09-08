import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import sys

#------------------------------------------------------------------------------#
import argparse
parser = argparse.ArgumentParser()
#
parser.add_argument('-L', action='store', dest='Level',
                    default=[-1], type=int, nargs='+',
                    help='level',
                    )
parser.add_argument('-M', action='store', dest='Mesh',
                    default=[0], type=int, nargs='+',
                    help='Mesh',
                    )
parser.add_argument('-P', action='store', dest='Points',
                    default=[0], type=int, nargs='+',
                    help='Points',
                    )
parser.add_argument('-F', action='store', dest='Files',
                    default=[''], type=str, nargs='+',
                    help='InputFiles',
                    )
Results = parser.parse_args()

Levels  =  np.array( Results.Level ) 
Levels  = -np.sort(-Levels)
if( np.any(Levels<0) ):
  parser.print_help()
  print "EXIT!!\n"
  sys.exit(1)

Mesh   = Results.Mesh[0]==1   
Points = Results.Points[0]==1

Fins   = Results.Files

import os.path
for f in Fins:
  if(f!=''):
    ok = os.path.exists(f)
    if(ok):
      print "'%s'" % f
    else:
      print "'%s' DOESNT EXIST!!" % f
      parser.print_help()
      print "EXIT!!\n"
      sys.exit(1)

#--------------------------------------------------------------------------||--#c
class PeanoHilbertIndex:
  def __init__(self, dim=3, level=1):
    self.nDim  = dim
    self.level = level

    self.GetIndex = self.GetIndex3D
    if(dim==2):
      self.GetIndex = self.GetIndex2D

    self.something = []

    self.Data = []
    self.i    = 0


  def Curve(self, x0, y0, xi, xj, yi, yj, n):
    if n <= 0:
        X = x0 + (xi + yi)/2
        Y = y0 + (xj + yj)/2
        self.Data.append( [Y, X, 0, self.i] )
        self.i += 1
    else:
        self.Curve(x0,               y0,               yi/2, yj/2, xi/2, xj/2, n - 1)
        self.Curve(x0 + xi/2,        y0 + xj/2,        xi/2, xj/2, yi/2, yj/2, n - 1)
        self.Curve(x0 + xi/2 + yi/2, y0 + xj/2 + yj/2, xi/2, xj/2, yi/2, yj/2, n - 1)
        self.Curve(x0 + xi/2 + yi,   y0 + xj/2 + yj,  -yi/2,-yj/2,-xi/2,-xj/2, n - 1)

    return self.Data

  def GetIndex3D(self, pt):
    index = 0
    x = pt.copy()
    for l in range(self.level):
      S    = np.where(x<0.5, 0, 1)
      S00  = (1-S[0])*(1-S[1])*(1-S[2])
      S12  =        (  S[1])*(1-S[2])
      S34  = (  S[0])*(1-S[1])
      S56  = (  S[0])*(  S[1])           # -2D
      S07  = (1-S[0])*(1-S[1])*(  S[2])  # -2D 

      idx  = (1+S[0]) * S12
      idx += (3+S[2]) * S34
      idx += (6-S[0]) * S56
      idx +=     7  * S07

      multiplier = np.power( 2, self.nDim * (self.level-1-l) )
      index += idx * multiplier

      transform    = np.zeros(3)
      transform[0] =(    ( 2.0 * x[0]    ) * S00
                   +     ( 2.0 * x[2]    ) * S12
                   +     (-2.0 * x[1]+1.0) * S34
                   + 2.0*(      -x[2]+1.0) * S56
                   + 2.0*(       x[0]    ) * S07 )

      transform[1] =(    ( 2.0 * x[2]    ) * S00
                   +     ( 2.0 * x[1]-1.0) * S12
                   + 2.0*(      -x[0]+1.0) * S34
                   +     ( 2.0 * x[1]-1.0) * S56
                   + 2.0*(      -x[2]+1.0) * S07 )

      transform[2] =(    ( 2.0 * x[1]    ) * S00
                   +     ( 2.0 * x[0]-S[0] ) * S12
                   +     ( 2.0 * x[2]-S[2] ) * S34
                   +     (-2.0 * x[0]+1.0) * S56
                   +     (-2.0 * x[1]+1.0) * S07 )

      x = transform.copy()

    return index

  def GetIndex2D(self, pt):
    index = 0
    x = pt.copy()
    something = []
    for l in range(self.level):
      S    = np.where(x<0.5, 0, 1)
      S00  = (1-S[0])*(1-S[1])*(1-S[2])
      S12  =          (  S[1])*(1-S[2])
      S03  = (  S[0])*(1-S[1])

      idx  = (1+S[0]) * S12
      idx +=       3  * S03

      multiplier = np.power( 2, self.nDim * (self.level-1-l) )
      index += idx * multiplier

      transform    = np.zeros(3)
      transform[0] =(    ( 2.0 * x[1]     ) * S00
                   +     ( 2.0 * x[0]-S[0]) * S12
                   +     (-2.0 * x[1]+1.0 ) * S03 )

      transform[1] =(    ( 2.0 * x[0]    ) * S00
                   +     ( 2.0 * x[1]-1.0) * S12
                   + 2.0*(      -x[0]+1.0) * S03 )

      x    = transform.copy()
#      aux  = np.floor(x*multiplier) 
#      aux +=(np.power(2,self.nDim*l) - 1) *  multiplier / 2 
#      something.append( aux.tolist() + [l] ) 
#
#    self.something.append( something ) 

    return index
#--------------------------------------------------------------------------||--#
class NETWORKX:
  def __init__(self):
    self.G01      = nx.Graph()
    self.NodeList = []
    self.LISTS    = {} 
    self.Boxes    = [] 


  def Set(self, COORDINATES, ELEMENTS):
    for idx, coord in enumerate(COORDINATES):
      pos = coord[:2]
      self.G01.add_node( idx+1, pos=pos )
      self.NodeList.append( idx+1 )
    self.pos = nx.get_node_attributes(self.G01,'pos')

    for element in ELEMENTS:
      n_vertices = len(element)
      edges = range(n_vertices)
      edges += [ edges[0] ]
      for i in range(n_vertices): 
        e =  edges[i], edges[i+1]  
        self.G01.add_edge( element[e[0]], element[e[1]] )

      Idx = element.astype(int) -1 
      Vtx =  COORDINATES[Idx,:] 
      vmax = np.amax(Vtx, axis=0) 
      vmin = np.amin(Vtx, axis=0) 
      box    =  vmax - vmin 
      center = (vmax + vmin)*0.5  
      self.Boxes.append( [center, box] )

  def SetParse(self, PARTSE):
    self.LISTS = {}
    for i in range(n_COORDINATES):
      part = PARTSE[i]
      if( self.LISTS.has_key(part) ):
        self.LISTS[part].append( i+1 ) # fortran stype 
      else:
        self.LISTS[part] = [i+1]

  def Draw(self):
    nx.draw_networkx_edges(self.G01, self.pos)
    COLOR = {0:'g', 1:'b', 2:'g'}
    for k,v in self.LISTS.iteritems():
      NodeList = v
      nx.draw_networkx_nodes(self.G01, pos, nodelist=self.NodeList, node_color=COLOR[k], node_size=75)

Nx = NETWORKX() 
#--------------------------------------------------------------------------||--#

#--------------------------------------------------------------------------||--#
# 
# +metis-5
#  mpmetis ELEMENTS01.dat 3
#  
# 
Fin01 = Fins[0]
COORDINATES   = np.loadtxt(Fin01)
n_COORDINATES = COORDINATES.shape[0]
print "n_COORDINATES:", n_COORDINATES

COORDINATES = COORDINATES[:,1:]
vmin = np.amin(COORDINATES, axis=0)

COORDINATES -= vmin
vmax = np.amax(COORDINATES, axis=0)
COORDINATES /= vmax

vmax = np.amax(COORDINATES, axis=0)
vmin = np.amin(COORDINATES, axis=0)
COORDINATES[:,2] = 0.0

if(Mesh): 
  Fin02 ="../MESHES/MESH01/A_ELEMENTS.alya"
  Fin02 = Fins[1]

  try:
    ELEMENTS   = np.loadtxt(Fin02, usecols=(1,2,3,4), skiprows=0)
  except IndexError:
    ELEMENTS   = np.loadtxt(Fin02, usecols=(1,2,3  ), skiprows=0)
  else: # no error occured
    ELEMENTS   = np.loadtxt(Fin02, usecols=(1,2,3,4), skiprows=0)
   #print "EXIT!!\n"
   #sys.exit(1)

  n_ELEMENTS = ELEMENTS.shape[0]
  print "n_ELEMENTS:", n_ELEMENTS

  Nx.Set(COORDINATES, ELEMENTS)


#--------------------------------------------------------------------------||--#
KEYs   = {}
CURVEs = {}
for idx, level in enumerate(Levels):
  if(level>=0):
    PH = PeanoHilbertIndex(dim=2, level=level)

    Keys = []
    for pt in COORDINATES:
      key = PH.GetIndex(pt)
      Keys.append( key )

    KEYs[level] = np.array( Keys[:] ).astype(int)

    CURVE      = np.array( PH.Curve(0.0, 0.0, 1.0, 0.0, 0.0, 1.0, level) )
    CURVE[:,2] = level

    CURVEs[level] = CURVE[Keys[:],:].copy()

min_level = np.amin(Levels)
max_level = np.amax(Levels)

PH = PeanoHilbertIndex(dim=2, level=max_level)

Keys = []
for pt in COORDINATES:
  key = PH.GetIndex(pt)
  Keys.append( key )

#--------------------------------------------------------------------------||--#
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
cmap = plt.cm.hot

fig, ax = plt.subplots()
fig.patch.set_alpha(0.0)
ax.patch.set_alpha(0.0)

n = max_level  
xticks, dx = np.linspace( vmin[0], vmax[0], num=2**n+1, retstep=True )
yticks, dy = np.linspace( vmin[1], vmax[1], num=2**n+1, retstep=True )

ax = fig.gca()
ax.grid( True, linestyle='--', linewidth=1.0, color='k')
ax.set_xticks( xticks ) 
ax.set_yticks( yticks ) 
ax.set_xticklabels( [] )
ax.set_yticklabels( [] )

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

plt.xlim(vmin[0]-0.05, vmax[0]+0.05)
plt.ylim(vmin[1]-0.05, vmax[1]+0.05)

## curve
colors = {'r': (1.0, 0.0, 0.0), 'g': (0.0, 0.5, 0.0), 'c': (0.0, 0.75, 0.75), 
          'b': (0.0, 0.0, 1.0), 'm': (0.75, 0, 0.75), 'y': (0.75, 0.75,0.00) }

for idx, level in enumerate(Levels):  
  if(level>0):
    CURVE = CURVEs[level]
    X  = CURVE[:,0]
    Y  = CURVE[:,1]
    Z  = CURVE[:,2]
    W  = CURVE[:,3]

    for x, y, z, w in zip(X,Y,Z,W):
      if(True): # occupied box  
        ax.add_artist( Rectangle(xy     = (x-dx/2,y-dy/2),
                                 color  = 'k',
                                 alpha  = 0.5, 
                                 width  = dx,
                                 height = dy) )      # Gives a square of area h*h

      if(False): # labels occupied box 
        ax.text(x, y, str( int(w) ).zfill(3),
                fontsize=6,  color='k', backgroundcolor='none', weight='bold',
                horizontalalignment='center', verticalalignment='center' )


    PH  = PeanoHilbertIndex(dim=2, level=level)
    r01 = np.array(PH.Curve( 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, PH.level) )

    segments = []
    for i,j in zip( r01[:-1,0:2], r01[1:,0:2]): segments.append( [i,j] ) 
    lc = LineCollection( segments, colors=colors.values()[idx], linestyle='solid', linewidths=2 ) 
#    ax.add_collection(lc) 

    if(False): # curve nodes labels 
      if(max_level<4): 
        for i in range(r01.shape[0]):
          x = r01[i,0]
          y = r01[i,1]
          j = r01[i,3]
          #
          ax.text(x, y, str( int(j) ).zfill(3),
                  fontsize=6,  color='k', backgroundcolor='none', weight='bold',
                  horizontalalignment='center', verticalalignment='center' )


## points 
if(Points): 
  plt.plot( COORDINATES[:,0], COORDINATES[:,1], 'o',  markersize=12.5, markeredgecolor='k', markerfacecolor='m')

  if(False): # nodes labels  
    if(max_level>0):
      for i in range(COORDINATES.shape[0]): 
        x = COORDINATES[i,0]
        y = COORDINATES[i,1]
        j = Keys[i]
        #
        ax.text(x, y, str(j).zfill(3), 
                fontsize=6,  color='darkgoldenrod', backgroundcolor='none', weight='bold', 
                horizontalalignment='center', verticalalignment='center' )

if(Mesh): 
#  if(True): Nx.Draw()

  if(True): # bounding box for every cell  
    for box in Nx.Boxes: 
      c,l = box 
      ax.add_artist( Rectangle(xy     = (c[0]-l[0]/2, c[1]-l[1]/2),
                                 color  = 'b',
                                 alpha  = 1.0,
                                 fill   = False, 
                                 lw     = 2, 
                                 width  = l[0],
                                 height = l[1]) )

     #plt.plot( c[0], c[1], '+',  markersize=12.5, markeredgecolor='b') #, markerfacecolor='m')

  if(False):
    Keys = []
    PH = PeanoHilbertIndex(dim=2, level=max_level)
    for box in Nx.Boxes:
      pt, l =  box 
      key = PH.GetIndex(pt)
      Keys.append( key )
     #plt.plot( pt[0], pt[1], 'o',  markersize=12.5, markeredgecolor='k', markerfacecolor='r')

    Keys = np.array( Keys[:] ).astype(int)

    Curve      = np.array( PH.Curve(0.0, 0.0, 1.0, 0.0, 0.0, 1.0, max_level) )  
    Curve[:,2] = max_level

  if(False): # occupied cells 
    for x,y,z,w in Curve[Keys,:]:
     #plt.plot( x,y, '+',  markersize=12.5, markeredgecolor='k') #, markerfacecolor='m')
      ax.add_artist( Rectangle(xy     = (x-dx/2,y-dy/2),
                                 color  = 'r',
                                 alpha  = 0.25,
                                 width  = dx,
                                 height = dy) )      # Gives a square of area h*h


plt.show()
#
Format = 'png' 
#print plt.gcf().canvas.get_supported_filetypes()
fig.savefig('space_filling_2dcurve01.%s' % Format, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight', format=Format) #, dpi=1200)
#
#--------------------------------------------------------------------------||--#

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
Results = parser.parse_args()

Levels  =  np.array( Results.Level ) 
Levels  = -np.sort(-Levels)
if( np.any(Levels<0) ):
  parser.print_help()
  print "EXIT!!\n"
  sys.exit(1)

Mesh   = Results.Mesh[0]==1   
Points = Results.Points[0]==1

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

  def Set(self, COORDINATES, ELEMENTS):
    for idx, coord in enumerate(COORDINATES):
      pos = coord[:2]
      self.G01.add_node( idx+1, pos=pos )
      self.NodeList.append( idx+1 )
    self.pos = nx.get_node_attributes(self.G01,'pos')

    for element in ELEMENTS:
      edge00 = element[0]
      edge01 = element[1]
      edge02 = element[2]
      self.G01.add_edge( edge00, edge01 )
      self.G01.add_edge( edge01, edge02 )
      self.G01.add_edge( edge02, edge00 )

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
COORDINATES   = np.loadtxt("../MESHES/MESH01/A_COORDINATES.alya")
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
  ELEMENTS   = np.loadtxt("../MESHES/MESH01/A_ELEMENTS.alya", usecols=(1,2,3))
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

    CURVEs[level] = CURVE[:,:].copy() #* 2**level  

max_level = np.amax(Levels)
min_level = np.amin(Levels)

PH = PeanoHilbertIndex(dim=2, level=max_level)

Keys = []
for pt in COORDINATES:
  key = PH.GetIndex(pt)
  Keys.append( key )
#--------------------------------------------------------------------------||--#
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d   import Axes3D


fig, ax = plt.subplots()
fig.patch.set_alpha(0.0)
ax.patch.set_alpha(0.0)

n = max_level  
xticks, dx = np.linspace( vmin[0], vmax[0], num=2**n+1, retstep=True )
yticks, dy = np.linspace( vmin[1], vmax[1], num=2**n+1, retstep=True )

ax = fig.gca( projection='3d' )
#ax.grid( True, linestyle='-', linewidth=1.0, color='k')
#ax.set_xticks( xticks ) 
#ax.set_yticks( yticks ) 
#
ax.set_xticklabels( [] )
ax.set_yticklabels( [] )
ax.set_zticklabels( [] )
#
# Hide axis 
ax.set_xticks([]) 
ax.set_yticks([]) 
ax.set_zticks([])
#

# Get rid of the panes
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

# Get rid of the spines
ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

dx, dy = 0.05, 0.05
#plt.xlim(vmin[0]-dx, vmax[0]+dx)
#plt.ylim(vmin[1]-dy, vmax[1]+dy)

## curve
colors = {'r': (1.0, 0.0, 0.0), 'g': (0.0, 0.5, 0.0), 'c': (0.0, 0.75, 0.75), 
          'b': (0.0, 0.0, 1.0), 'm': (0.75, 0, 0.75), 'y': (0.75, 0.75,0.00) }

dx = 1.0/2.0**((max_level-1))
print dx, 2**(max_level) 

for idx, level in enumerate(Levels):  
  if(level-1>=min_level):
    r01= CURVEs[level] * 2**max_level 
    ok =  np.logical_and( r01[:,0] >=2, r01[:,1] >=2 )  
    ok =  np.logical_and( r01[:,0] <=2**max_level-2, ok )
    ok =  np.logical_and( ok, r01[:,1] <=2**max_level-2 )
    ok = ~ok
    ax.plot( r01[:,0], r01[:,1], r01[:,2],
             '-', color=colors.values()[level],
             markersize=5., markeredgecolor='k', markerfacecolor='m')
 
    ax.plot( r01[ok,0], r01[ok,1], r01[ok,2], 
             'o', color=colors.values()[level],  
             markersize=4., markeredgecolor='k', markerfacecolor='b')

    for i,j,k in zip( r01[1:,0:3], r01[:-1,0:3], ok[:]): 
      if(k):
        x, y, z = ([i[0],j[0]]), ([i[1],j[1]]), ([i[2],j[2]])
        ax.plot( x, y, '-', lw=1.0, color='k', zdir='z', zs=0.0)
        ax.plot( x, y, 'o', zdir='z', zs=0.0,
             color=colors.values()[level],
             markersize=4., markeredgecolor='k', markerfacecolor='b')


    r01= CURVEs[level-1] * 2**max_level 
    ok =  np.logical_and( r01[:,0] >=2, r01[:,1] >=2 )
    ok =  np.logical_and( r01[:,0] <=2**max_level-2, ok )
    ok =  np.logical_and( ok, r01[:,1] <=2**max_level-2 )

    ax.plot( r01[:,0], r01[:,1], r01[:,2],
             '-', color=colors.values()[level],
             markersize=5.0, markeredgecolor='k', markerfacecolor='m')

    ax.plot( r01[ok,0], r01[ok,1], r01[ok,2],
             'o', color=colors.values()[level],
             markersize=4.0, markeredgecolor='k', markerfacecolor='r')

    for i,j,k in zip( r01[1:,0:3], r01[:-1,0:3], ok[0:]):   
      if(k):
        x, y, z = ([i[0],j[0]]), ([i[1],j[1]]), ([i[2],j[2]])
        ax.plot( x, y, '-', lw=1.0, color='k', zdir='z', zs=0.0)
        ax.plot( x, y, 'o', zdir='z', zs=0.0, 
             color=colors.values()[level],
             markersize=4.0, markeredgecolor='k', markerfacecolor='r')

ax.view_init(elev=25, azim=40)

plt.show()

#
Format = 'png' 
#print plt.gcf().canvas.get_supported_filetypes()
fig.savefig('space_filling_3dcurve02.%s' % Format, format=Format, 
            facecolor=fig.get_facecolor(), 
            edgecolor='none', bbox_inches='tight') #, dpi=1200)
#
#--------------------------------------------------------------------------||--#

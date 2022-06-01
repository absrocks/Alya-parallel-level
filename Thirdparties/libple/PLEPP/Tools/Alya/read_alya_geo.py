
#--------------------------------------------------------------------------||--#
def Read_alya_geo(fname, INIT_KEY="COORDINATES"):
    import glob 
    fname = glob.glob(fname) 
    if(not len(fname)<1):
      fname = fname[-1]
      print "|_\'%s\'" % fname
    else:
      print "Error: \'%s\'! \n" % fname 
      exit(1)

    data = open(fname, "r")
    lines = data.readlines()
    data.close()

    nline = len(lines)

    #global INIT_KEY
    INIT_KEY = INIT_KEY.replace("_", "") 
    INIT_KEY = INIT_KEY.replace("-", "") 
    INIT_KEY = INIT_KEY.replace("&", "") 
    END_KEY  = "END" + INIT_KEY 

    ok  = False 
    IDs = []
    for i in range(nline):
      line = lines[i]
      if(not line.find(INIT_KEY)<0): IDs.append(i+1)
      if(not line.find(END_KEY)<0):  IDs.append(i+0) 

    XYZ = []      
    for i in range(IDs[0], IDs[1]-1):
      line = lines[i]
      line = line.strip() 
      line = line.split()
      XYZ.append([eval(val) for val in line[1:]]) 
    
    print "  |_No elements:", len(XYZ)
    return XYZ 
#--------------------------------------------------------------------------||--#

#--------------------------------------------------------------------------||--#
def Alya_cell_types():
    # col=1: kernel/domain/elmtyp.f90        <<==
    # col=2: kernel/elsest/elsest_geogid.f90 <<=
    alya_types = []
    alya_types.append( ("TRI03", 10) )
    alya_types.append( ("QUA04", 12) )
    alya_types.append( ("TET04", 30) )
    alya_types.append( ("PYR05", 32) )
    alya_types.append( ("PEN06", 34) )
    alya_types.append( ("HEX08", 37) )

    return alya_types

def Saturne_cell_types():
  #SATURNE: ./src/fvm/fvm_defs.h -> fvm_element_t
  fvm_element_t = [] 
  #fvm_element_t.apppend( ("FVM_EDGE",             0) )
  fvm_element_t.append( ("FVM_FACE_TRIA",        1) ) 
  fvm_element_t.append( ("FVM_FACE_QUAD",        2) ) 
  #fvm_element_t.apppend( ("FVM_FACE_POLY",        3) )  
  fvm_element_t.append( ("FVM_CELL_TETRA",       4) )  
  fvm_element_t.append( ("FVM_CELL_PYRAM",       5) ) 
  fvm_element_t.append( ("FVM_CELL_PRISM",       6) )    
  fvm_element_t.append( ("FVM_CELL_HEXA",        7) )  
  #fvm_element_t.apppend( ("FVM_CELL_POLY",        8) )   
  #fvm_element_t.apppend( ("FVM_N_ELEMENT_TYPES",  9) )   
  return fvm_element_t


def Vtk_cell_types():
    vtk_types = []    #(vtk_type, vtk_id, vtkCell) <= linear cell types found in VTK 
    vtk_types.append( ("VTK_TRIANGLE",     5, vtk.vtkTriangle()   ) )
    vtk_types.append( ("VTK_QUAD",         9, vtk.vtkQuad()       ) )
    vtk_types.append( ("VTK_TETRA",       10, vtk.vtkTetra()      ) )
    vtk_types.append( ("VTK_PYRAMID",     14, vtk.vtkPyramid() ) )
    vtk_types.append( ("VTK_WEDGE",       13, vtk.vtkWedge()       ) )
    vtk_types.append( ("VTK_HEXAHEDRON",  12, vtk.vtkHexahedron() ) )

    return vtk_types

class Vtk2Alya_cell():
  def __init__(self):
    alya_types = Alya_cell_types()
    vtk_types  = Vtk_cell_types()

    self.alya2vtk = {}
    self.vtk2alya = {}
    for alya_type, vtk_type in zip(alya_types, vtk_types):
      alya_idx = alya_type[1]
      self.alya2vtk[ alya_idx ] = vtk_type

      vtk_idx  = vtk_type[1]
      self.vtk2alya[ vtk_idx ]  = alya_idx


  def get_vtk_cell(self, alya):
    vtk = self.alya2vtk.get(alya, -1)
    return vtk


  def get_alya_cell(self, vtk):
    alya = self.vtk2alya.get(vtk, -1)
    return alya


class Get_cell_type():
  def __init__(self):
    alya_types    = Alya_cell_types()
    saturne_types = Saturne_cell_types()
    #print saturne_types

    self.alya2saturne = {}
    self.n_saturne = {}
    for alya_type, saturne_type in zip(alya_types, saturne_types):
      alya_idx    = alya_type[1]
      saturne_idx = saturne_type[1]
      self.alya2saturne[ alya_idx ] = saturne_idx

      self.n_saturne[saturne_idx] = [0, saturne_type[0]] 
    
    self.cosa = ""

  def get_saturne_cell(self, alya):
    saturne = self.alya2saturne.get(alya, -1)
    self.n_saturne[saturne][0] += 1

    return saturne


  def print_types(self):
    saturne_types = Saturne_cell_types()

    tot = 0
    for key, val in self.n_saturne.items(): 
      if(val[0]>0): 
        self.cosa += "\'%s\': %d "% (val[1], val[0]) 
        tot += val[0]
    self.cosa += "\'FVM_CELLs\': %d" % tot
    return self.cosa

#--------------------------------------------------------------------------||--#


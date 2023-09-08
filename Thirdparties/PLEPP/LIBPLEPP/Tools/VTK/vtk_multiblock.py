import vtk 
import os
import sys 
import vtk_tools as tools 


def __times__( eg ):
    times = eg.GetTimeSets().GetItem(0)
    print times.GetNumberOfTuples()
  

#--------------------------------------------------------------------------||--#
#-----------------------------------------------------------------| Ensight |--#
#--------------------------------------------------------------------------||--#
class Vtk_ensight_reader():
  def __init__(self):
    self.reader = vtk.vtkGenericEnSightReader() # vtk.vtkEnSightGoldReader()

  def init(self, filename=""):
    self.filename, self.fileExtension = os.path.splitext( filename )
    print "[Vtk_ensight_reader] <- \'%s\'" % (filename)

    cdp     = vtk.vtkCompositeDataPipeline()
    self.eg = self.reader
    self.eg.SetDefaultExecutivePrototype(cdp) 

    self.eg.SetCaseFileName(filename) 
    self.eg.ReadAllVariablesOn() 
    self.eg.Update()

    self.Obj = self.eg.GetOutput()

    __times__( self.eg ) 


  def set_composite(self): 
    self.Merge        = vtk.vtkMultiBlockMergeFilter()
    self.AppendFilter = vtk.vtkAppendFilter()
    #self.AppendFilter = vtk.vtkAppendPolyData()
    tools.__set_composite__(self.Obj, self.Merge, self.AppendFilter)


  def save_vtu(self):
    print "[Vtk_ensight_reader]", 
    tools.__save_vtu__(self.filename, self.AppendFilter.GetOutput())


  def save_vtp(self):
    print "[Vtk_ensight_reader]", 
    tools.__save_vtp__(self.filename, self.AppendFilter.GetOutput())


  def save_multiblock(self):
    print "[Vtk_ensight_reader]", 
    tools.__save_vtm__(self.filename, self.eg.GetOutput())


  def BlockIdScalars(self):
    print "[Vtk_ensight_reader]", 

    #print "[BlockIdScalars]",
    #print self.Obj.GetNumberOfBlocks(), 

    ymi = tools.__youngs_material_interface__( self.Obj ) 
    #tools.__vtkBlockIdScalars__( self.Obj )

    print "<-----"

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------| MultiBlock |--#
#--------------------------------------------------------------------------||--#
class Vtk_multiblock_reader():
  def __init__(self):
    self.reader = vtk.vtkXMLMultiBlockDataReader() 


  def init(self, filename=""):
    self.filename, self.fileExtension = os.path.splitext( filename )
    print "[Vtk_multiblock_reader] <- \'%s\'" % (filename)

    self.reader.SetFileName(filename)
    self.reader.Update()
    self.Obj = self.reader.GetOutput()

  def set_composite(self):  
    self.Merge        = vtk.vtkMultiBlockMergeFilter()
    self.AppendFilter = vtk.vtkAppendFilter()
    tools.__set_composite__(self.Obj, self.Merge, self.AppendFilter)

    #unstructuredGrid = vtk.vtkUnstructuredGrid() 
    #unstructuredGrid.ShallowCopy( self.Merge.GetOutput().GetBlock(0) )


  def save_vtp(self):
    print "[Vtk_multiblock_reader]", 
 
    filename  = self.filename
    filename += "_merge"
    tools.__save_vtp__(filename, self.Merge.GetOutput().GetBlock(0) ) 


  def save_vtu(self):
    print "[Vtk_multiblock_reader]", 
    tools.__save_vtu__(self.filename, self.AppendFilter.GetOutput() ) 


#--------------------------------------------------------------------------||--#
#----------------------------------------------------------------| Examples |--#
#--------------------------------------------------------------------------||--#

#--------------------------------------------------------------------------||--#
if(False): 
  Fin  = "SURFACE.vtm"

  MB = Vtk_multiblock_reader()
  MB.init(Fin) 
  MB.set_composite()

  MB.save_vtp() 
  MB.save_vtu() 

#--------------------------------------------------------------------------||--#
if(False):
  Fin  = "star.case"

  ENSI = Vtk_ensight_reader()
  ENSI.init(Fin) 
  ENSI.set_composite() 

  ENSI.save_multiblock() 
  ENSI.save_vtu() 

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#

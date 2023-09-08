import numpy as np 
import smesh
import os
import sys

from my_writer import *
from alya_create_field import *


class Salome2dat:
  def __init__(self):
    self.points = []
    self.groups = [] 

    self.num_pts    = -1
    self.num_vols   = -1
    self.num_faces  = -1
    self.num_groups = -1
    self.num_bndrys = -1

    self.vols_ids   = [] 
    self.vols_type  = [] 
    self.faces_ids  = [] 
    self.faces_type = []  
    self.faces_grp  = [] 
    #self.pts_grp    = [] 

    self.groups_size = {}
    self.groups_names = {}
    self.entity_3d   = {} 
    self.entity_type = {}
    self.vtk_pts_ok  = False 
         
    self.VTK  = VTK_writer()
    self.Alya = Alya_field()
    
    self.without_group = -1 
    
    print "\n +Salome2dat"


  def get_data(self, mesh):
    self.__mesh_info__( mesh )

    self.__get_points__(  mesh, self.points)
    self.__get_entities__(mesh, smesh.VOLUME, self.vols_ids,  self.vols_type)

    #faces_ids_aux = [] 
    #self.__get_entities__(mesh, smesh.FACE,   faces_ids_aux, self.faces_type)
    self.__get_entities__(mesh, smesh.FACE,   self.faces_ids, self.faces_type)
    
    self.__get_group_for_each_face__(mesh, self.faces_grp)
    self.__faces_by_group_correction__(self.faces_ids, self.faces_grp) 
    
    self.__get_pts_groups__(mesh)

    self.VTK.num_pts = len(self.points)
    print "   |_Salome2dat.get_data OK!"


  def delete_face_grup(self, grup_name): 
    grup_idx = self.groups_names.get(grup_name, -1)
    if(grup_idx==-1): 
      print "ERROR:", grup_name, "doesnt exist!! \n"
      sys.exit()

    OK = []
    n_by_deletes = 0
    n_faces_grp = len( self.faces_grp ) 
    for i in range(n_faces_grp):
      grp = self.faces_grp[i]
      if(grp==grup_idx): 
        n_by_deletes += 1
      else:
        OK.append(i)

    self.num_faces -= n_by_deletes 

    self.faces_ids = [ self.faces_ids[ok] for ok in OK  ] 
    self.faces_grp = [ self.faces_grp[ok] for ok in OK  ]

    print "   |_Salome2dat.delete_face_grup \'%s\', FACES:" %(grup_name), 
    print self.num_faces+n_by_deletes, "->", self.num_faces
    self.groups_names[grup_name] *= -1 


  def get_npoints(self):
    return self.num_pts 


  def get_groups(self):
    return np.array( self.groups )


  def get_points(self):
    return np.array( self.points )


  def set_scale(self, dx=1.0, dy=1.0, dz=1.0):
    for i in range(self.num_pts): 
      x,y,z = self.points[i]
      self.points[i] = [x*dx,y*dy,z*dz]
      #print x*dx,y*dy,z*dz
      #print self.points[i] 
    print "   |_Salome2dat.set_scale: (%f,%f,%f) " % (dx,dy,dz)


  def to_vtk(self, file_name, file_type=""):
    self.file_name = os.path.splitext(file_name)[0]
    
    if(file_type==""): file_type="vols" 
    
    self.__cells_types__()

    if(not self.vtk_pts_ok):
      for point in self.points: self.VTK.set_point(point) 
      self.vtk_pts_ok = True 
    
    if(file_type=="faces"): 
      for face_ids, face_type in zip(self.faces_ids, self.faces_type): 
        vtk_cell = self.salome2vtk.get(face_type, "") 
        self.VTK.set_cell(face_ids, vtk_cell)

    if(file_type=="vols"): 
      for vol_ids, vol_type in zip(self.vols_ids, self.vols_type): 
        vtk_cell = self.salome2vtk.get(vol_type, "??") 
        #print vol_type, vtk_cell 
        self.VTK.set_cell(vol_ids, vtk_cell)

    self.VTK.set_points_prop(self.groups, "Ids") 

    self.VTK.to_vtkunstructured_grid() 
    self.VTK.writer(self.file_name+"_"+file_type+".vtk")     


  def to_alya(self, file_name):
    self.file_name = os.path.splitext(file_name)[0]

    self.__cells_types__()

    alya_types = [] 
    for vol_type in self.vols_type: alya_types.append( self.salome2alya[vol_type][1] ) 

    header = ["TYPES", "END_TYPES"]
    self.__save_raw__(alya_types,  self.file_name+"_TYPES.alya", header)
    

    a=0; b=1; c=2; d=3
    ids_tetra = [b,c,d,a]

    a,b,c,d,e,f,g,h = 0,1,2,3,4,5,6,7
    ids_hexa = [e,f,g,h,a,b,c,d]

    a,b,c,d,e,f = 0,1,2,3,4,5
    ids_penta = [d,e,f,a,b,c] 

    a,b,c,d,e = 0,1,2,3,4
    ids_pyram = [b,c,d,e,a] 

    aux = [] 
    for nodes, kind in zip(self.vols_ids,self.vols_type):  
      nodes = np.array(nodes)
      #
      if(kind==smesh.Entity_Tetra): 
        nodes = nodes[ids_tetra]
        self.VTK.get_tetra_vol(nodes-1, self.points)
      if(kind==smesh.Entity_Hexa): 
        nodes = nodes[ids_hexa]      
      if(kind==smesh.Entity_Penta): 
        nodes = nodes[ids_penta]
      if(kind==smesh.Entity_Pyramid): 
        nodes = nodes[ids_pyram]
      #      
      aux.append( nodes.tolist() )  
    
    header = ["ELEMENTS", "END_ELEMENTS"]
    self.__save_raw__(aux,   self.file_name+"_ELEMENTS.alya", header) 
    #print self.vols_ids 
    #aux01 = np.array( self.vols_ids, int )  
    #aux02 = np.zeros(aux01.shape, int)
  
    header = ["COORDINATES", "END_COORDINATES"]
    self.__save_raw__(self.points,     self.file_name+"_COORDINATES.alya", header) 

    header = ["BOUNDARIES", "END_BOUNDARIES"]
    self.__save_raw__(self.faces_ids,  self.file_name+"_BOUNDARIES.alya", header) 

    header = ["ON_BOUNDARIES", "END_ON_BOUNDARIES"] # grps
    self.__save_raw__(self.faces_grp, self.file_name+"_ON_BOUNDARIES.alya", header) 

    header = ["CHARACTERISTICS", "END_CHARACTERISTICS"] 
    #self.__save_raw__(self.groups, self.file_name+"_CHARACTERISTICS.alya", header) 
    self.__save_raw__(self.groups, self.file_name+"_CHARACTERISTICS.alya") 

    # xxx.nsa.dat 
    alya_types = [self.salome2alya[mtype][0] for mtype, mok in self.entity_3d.items() if(mok)]
    
    self.Alya.ELEMENTS         = self.num_vols
    self.Alya.BOUNDARIES       = self.num_faces # self.num_bndrys
    self.Alya.TYPES_OF_ELEMS   = " ".join(alya_types)
    self.Alya.SPACE_DIMENSIONS = 3
    self.Alya.GRUPS_NAMES      = self.groups_names
    self.Alya.GRUPS_SIZES      = self.groups_size 
    self.Alya.write(file_name)


  def save_data(self, file_name):
    self.file_name = os.path.splitext(file_name)[0]

    self.__save_raw__(self.points,    self.file_name+"_pts.dat") 
    self.__save_raw__(self.faces_ids, self.file_name+"_faces.dat") 
    self.__save_raw__(self.vols_ids,  self.file_name+"_vols.dat") 

    self.__save_raw__(self.groups,    self.file_name+"_faces_grp.dat") 
    self.__save_raw__(self.faces_type, self.file_name+"_faces_type.dat") 
    self.__save_raw__(self.vols_type,  self.file_name+"_vols_type.dat") 
    
    
  def __cells_types__(self): 
    self.salome_types = []
    self.salome_types.append( smesh.Entity_Triangle   )
    self.salome_types.append( smesh.Entity_Quadrangle )
    self.salome_types.append( smesh.Entity_Tetra      )
    self.salome_types.append( smesh.Entity_Pyramid    )
    self.salome_types.append( smesh.Entity_Penta      ) 
    self.salome_types.append( smesh.Entity_Hexa       )

    # col=1: kernel/domain/elmtyp.f90        <<==
    # col=2: kernel/elsest/elsest_geogid.f90 <<=
    self.alya_types = []
    self.alya_types.append( ("TRI03", 10) ) 
    self.alya_types.append( ("QUA04", 12) )
    self.alya_types.append( ("TET04", 30) ) 
    self.alya_types.append( ("PYR05", 32) )  
    self.alya_types.append( ("PEN06", 34) ) 
    self.alya_types.append( ("HEX08", 37) )  


    self.salome2vtk = {} 
    for salome_type, vtk_type in zip(self.salome_types, self.VTK.vtk_types): 
      self.salome2vtk[salome_type] = vtk_type 
      
    self.salome2alya = {} 
    for salome_type, alya_type in zip(self.salome_types, self.alya_types): 
      self.salome2alya[salome_type] = alya_type 
      

  def __save_raw__(self, data, file_name, header=["",""]):
    file_out = open(file_name, "w")
    if(header[0]!=""): print >> file_out, header[0]

    for i in range( len(data) ):  
      if(header[0]!=""): print >> file_out, i+1, 
      
      if(isinstance(data[i], list)): 
        for cosa in data[i]: print >> file_out, cosa, 
        print >> file_out  
      else:
        print >> file_out, data[i]

    if(header[1]!=""): print >> file_out, header[1]
    file_out.close() 

    print "   |_Salome2dat: \'%s\' " % file_name 


  def __get_pts_groups__(self, mesh):    
    cells_idx   =  mesh.GetElementsByType(smesh.FACE)   
    self.groups = -np.ones( mesh.NbNodes(), int )

    self.__get_groups__(mesh)  # >> self.cells_group.get
    for cell_id in cells_idx:
      group_idx = self.cells_group.get(cell_id, self.without_group)
      if(group_idx!=self.without_group):  
        nodes_idx = mesh.GetElemNodes(cell_id)
        for node_id in nodes_idx: self.groups[node_id-1] = group_idx #<< group for each node 

    self.num_bndrys = 0
    for idx in self.groups: 
      if(idx!=self.without_group): self.num_bndrys += 1



  def __get_groups__(self, mesh): 
    groups_names = mesh.GetGroupNames() 
    groups       = mesh.GetGroups() 
    self.cells_group = {}
    
    group_idx = 1
    for idx in range(mesh.NbGroups()): 
      group_ids = groups[idx].GetListOfID()
      print "    |_%d) \'%s\'" % (group_idx, groups_names[idx]), 
      print len( group_ids )

      self.groups_names[groups_names[idx]] = group_idx
      self.groups_size[groups_names[idx]]  = len( group_ids )

      for cell_id in group_ids: self.cells_group[cell_id] = group_idx
      group_idx += 1

    print "   |_Salome2dat.get_groups %d" % (group_idx-1)


  def __get_points__(self, salome_mesh, array):
    for idx in salome_mesh.GetNodesId(): array.append( salome_mesh.GetNodeXYZ(idx) )
    print "   |_Salome2dat.get_entity \'NODE\' %d" % self.VTK.num_pts


  def __mesh_info__(self, mesh):
    print " |_Elements:", mesh.NbElements()
    print " |_Vols: ", mesh.NbVolumes(), 
    print " (Linear %d)" % ( mesh.NbVolumesOfOrder(smesh.ORDER_LINEAR) )
    print " |_Faces:", mesh.NbFaces(),
    print " (Linear %d)" % ( mesh.NbFacesOfOrder(smesh.ORDER_LINEAR) )

    self.entity_3d[smesh.Entity_Tetra] = False
    self.entity_3d[smesh.Entity_Penta] = False
    self.entity_3d[smesh.Entity_Hexa]  = False
    self.entity_3d[smesh.Entity_Pyramid]  = False

    mesh_info = mesh.GetMeshInfo() 
    for mesh_type, mesh_size in mesh_info.items(): 
      if(mesh_size!=0): 
        print "   |_\'%s\' %d" % (mesh_type, mesh_size)
        self.entity_type[mesh_type] = mesh_size 
        if(not self.entity_3d.get(mesh_type,True)): self.entity_3d[mesh_type] = True

    #for (key, val) in self.entity_3d.items(): 
    #  if(val==True):
    #    if(key==smesh.Entity_Pyramid): 
    #      print "ERROR: \'%s\' entity unsupported!!\n\n" % key 
    #      sys.exit()

    print "  _|"
    print " |_SubMesh:", mesh.NbSubMesh()
    print " |_Groups:",  mesh.NbGroups() 

    self.num_pts    = mesh.NbNodes()
    self.num_vols   = mesh.NbVolumes()
    self.num_faces  = mesh.NbFaces()    
    self.num_groups = len(mesh.GetGroupNames())

    self.VTK.num_pts       = self.num_pts 
    self.Alya.NODAL_POINTS = self.num_pts
    
    
  def __get_entities__(self, salome_mesh, element_type, entities_nodes, entities_type=[]):
    elements_ids = salome_mesh.GetElementsByType(element_type)
    for idx in elements_ids:
      entity_nodes = salome_mesh.GetElemNodes(idx)
      entities_nodes.append( entity_nodes )    
      
      if(element_type!=smesh.NODE):
        entity_type  = salome_mesh.GetElementGeomType(idx)
        entities_type.append(  entity_type )

    print "   |_Salome2dat.get_entity \'%s\' %d" % (element_type, len(entities_nodes))


  def __get_group_for_each_face__(self, salome_mesh, faces_grp):  
    elements_id =  salome_mesh.GetElementsByType(smesh.FACE)   

    self.__get_groups__(salome_mesh)  # >> self.cells_group.get
    for idx in elements_id:
      group_id = self.cells_group.get(idx, self.without_group)
      #entity_nodes = salome_mesh.GetElemNodes(idx)
      #print idx, group_id,entity_nodes 
      #faces_grp.append( entity_nodes+[group_id, idx] )
      faces_grp.append( group_id )

  
  def __faces_by_group_correction__(self, faces, faces_grps):
    if(len(faces) != len(faces_grps)): 
      print "ERROR: size face != size face groups, %d!=%d !! \n\n" % (len(faces), len(faces_grps))
      sys.exit() 
    
    self.faces_ids = []
    self.faces_grp = []
    for face, grp in zip(faces, faces_grps):
      if(grp>0): 
        self.faces_ids.append(face) 
        self.faces_grp.append(grp) 
      else:   
        self.num_faces -= 1

    grp_size = sum(self.groups_size.values())
    if(self.num_faces != grp_size): 
      print "ERROR: size face != size face groups, %d!=%d !! \n\n" % (self.num_faces, grp_size)
      print self.groups_names.items() 
      sys.exit() 

    print "   |_Salome2dat.faces_correction \'%s\' changed: %d!" % ("FACE", self.num_faces)


  def periodic_boundaries_z(self, inlet_name, outlet_name, file_name, error=1e-3):
    """
    ->outlet->inlet->DOMAIN->outlet->inlet->
    ->master->slave->DOMAIN->master->slave->
    inlet_name  -> master:   
    outlet_name -> slave: 
    """
    #self.Alya.periodic_boundaries(inlet_name, outlet_name) 
    file_name = os.path.splitext(file_name)[0]    
    
    inlet_id  = self.groups_names.get(inlet_name,  -1) 
    outlet_id = self.groups_names.get(outlet_name, -1)     
    if(inlet_id==-1 or outlet_id==-1): 
      print "ERROR: there not exist %s/%s key!! \n\n" % (inlet_name, outlet_name)
      print self.groups_names 
      sys.exit() 

    inlets  = [] 
    outlets = []  
    for idx in range(self.num_pts):
      grp = self.groups[idx] 
      if(grp==inlet_id):  inlets.append(idx) 
      if(grp==outlet_id): outlets.append(idx) 
    
    if(not len(outlets) == len(inlets)): 
      print "ERROR: size %s/%s %d!=%d!! \n\n" % (inlet_name, outlet_name, len(outlets), len(inlets))
      sys.exit() 


    periodicity = {}
    for i in range(len(inlets)):
      master = inlets[i]
      slave  = outlets[i]
      
      #master_pt = self.points[master]
      #slave_pt  = self.points[slave]
      #if(master_pt[:-1]==slave_pt[:-1]): 

      master_pt = np.array(self.points[master]) 
      slave_pt  = np.array(self.points[slave])
      diff = master_pt[:-1] - slave_pt[:-1] 
      mod = np.sqrt( (diff**2).sum() ) < error
      if(mod):  

        periodicity[master] = slave
        inlets[i]  *= -1
        outlets[i] *= -1

    print 
    for i in range(len(inlets)):
      master = inlets[i]
      if(master>=0):
        #print "%d)"% i ,master, 
        for slave in outlets:
          if(slave>0):
            master_pt = np.array(self.points[master]) 
            slave_pt  = np.array(self.points[slave])
            
            diff = master_pt[:-1] - slave_pt[:-1] 
            mod = np.sqrt( (diff**2).sum() ) < error
            
            #if(master_pt[:-1]==slave_pt[:-1]):
            if(mod):  
              periodicity[master] = slave
              inlets[i]  *= -1
              #print slave, 
        #print 

    file_out = open(file_name+"_periodic.alya", "w")
    for i in range(len(inlets)):
      master = -inlets[i]
      slave  = periodicity.get(master, -1)
      print>> file_out, master+1, slave+1

      if(master<0): 
        print "ERROR: there not exist symmetry for point %d !!" % master, 
        print self.points[abs(master)], 
        print 
        sys.exit()       
      
    file_out.close()

    self.Alya.PERIODIC_NODES = len(inlets)

    file_out = open(file_name+"_periodic_master.dat", "w")
    for i in range(len(inlets)):
      master = -inlets[i]
      slave  = periodicity.get(master, -1)
      print>> file_out, master+1, 
      for xi in self.points[abs(master)]: print>> file_out, xi, 
      print>> file_out
    file_out.close()

    file_out = open(file_name+"_periodic_slave.dat", "w")
    for i in range(len(inlets)):
      master = -inlets[i]
      slave  = periodicity.get(master, -1)
      print>> file_out, slave+1, 
      for xi in self.points[abs(slave)]: print>> file_out, xi, 
      print>> file_out
    file_out.close()


    print "   |_Salome2dat.periodic_boundaries \'%s\' " % (file_name+"_periodic.alya")


#============================================================================# 

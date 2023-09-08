NSI_DAT = """$$Create by J. MIGUEL ZAVALA-AKE 
$--------------------------------------------------------XXX_nsi.dat 
RUN_DATA
  CODE:               1
  ALYA:               %s 
  MEMORY:             yes
END_RUN_DATA
$-------------------------------------------------------------------
PROBLEM_DATA
  TIME_COUPLING:      Global, From_critical
  TIME_INTERVAL  =    0.0, 1.0
  TIME_STEP_SIZE =    1.0e-3
  NUMBER_OF_STEPS=    %d 
$
$ MAXIMUM_NUMBER_GLOBAL: -1 -2 -3  
$ BLOCK_ITERATION: -3 
$ END_BLOCK_ITERATION
$
  NASTIN_MODULE:      ON
  END_NASTIN_MODULE
$
$  XXXXX_MODULE:      OFF
$  END_TEMPER_MODULE
$
  PARALL_SERVICE:     ON
\t PARTITION_TYPE:       FACES
\t$COMMUNICATIONS:       SYNCHRONOUS
\t$FILE_HIERARCHY:       ON
\t$FILE_OPEN_CLOSE:      Yes
\t$VIRTUAL_FILE:         On, MAXIMUM_MEMORY=0.5
\t$$TASK:                READ_PREPROCESS, BINARY                 $ TO RUN
\t$TASK:                 ONLY_PREPROCESS, BINARY, SUBDOMAIN = %d $ N_DOMs-1
  END_PARALL_SERVICE
$
 $HDFPOS_SERVICE:     ON
 $END_HDFPOS_SERVICE
END_PROBLEM_DATA
$-------------------------------------------------------------------
"""

NSI_NSIDAT = """$$Create by J. MIGUEL ZAVALA-AKE 
$------------------------------------------------------------
PHYSICAL_PROBLEM
$
  PROBLEM_DEFINITION       
    TEMPORAL_DERIVATIVES:   On 
    CONVECTIVE_TERM:        On
    VISCOUS_TERM:           FULL
    REGIME:                 INCOMPRESSIBLE
   $REGIME:                 LOWMACH, TEMPERATURE = 298.0, PRESSURE=101325.0
   $GRAVITY:                GZ=
  END_PROBLEM_DEFINITION  
$
  PROPERTIES
$  INCOMPRESSIBLE:
     DENSITY                 1000.0
     VISCOSITY               0.001   
     PRESSURE                0.0
$  LOWMACH:
$    DENSITY            = EXTERNAL 
$    LAW_VISCOSITY      = EXTERNAL
$    THERMODYN_PRESSURE = 101325.0
$    GAS_CONSTANT       = 287.0 $$ 8.3144621
  END_PROPERTIES  
$
END_PHYSICAL_PROBLEM  

$------------------------------------------------------------
NUMERICAL_TREATMENT 
  STABILIZATION:            ASGS
  ELEMENT_LENGTH:           Minimum
  TIME_INTEGRATION:         Trapezoidal, ORDER: 2, EULER=5
  SAFETY_FACTOR:            1.0 
  STEADY_STATE_TOLER:       1e-12
  NORM_OF_CONVERGENCE:      LAGGED_ALGEBRAIC_RESIDUAL
  MAXIMUM_NUMBER_OF_IT:     1
  CONVERGENCE_TOLERANCE:    1e-03
  ALGORITHM:                SCHUR
    SOLVER:                 ORTHOMIN, CONTINUITY_PRESERVING
    PRECONDITIONER:         DT
    ELEMENT_LENGTH:         MINIMUM
    TAU_STRATEGY:           CODINA
  END_ALGORITHM
  MOMENTUM
    ALGEBRAIC_SOLVER:       GMRES,  ITERA=1000, TOLER=1.0e-10, ADAPTIVE, RATIO=0.001, KRYLOV=50, CONVERGENCE
    PRECONDITIONER:         DIAGONAL
  END_MOMENTUM 
  CONTINUITY
    ALGEBRAIC_SOLVER:       CG,     ITERA=1000, TOLER=1.0e-10, ADAPTIVE, RATIO=0.001, UNSYMME, CONVERGENCE
    PRECONDITIONER:         DIAGONAL 
  END_CONTINUITY
END_NUMERICAL_TREATMENT  

$------------------------------------------------------------
OUTPUT_&_POST_PROCESS  
  POSTPROCESS  VELOCITY
  POSTPROCESS  PRESSURE
END_OUTPUT_&_POST_PROCESS   

$------------------------------------------------------------
BOUNDARY_CONDITIONS
$
  PARAMETERS
    VARIATION:         CONSTANT
    INITIAL:           INERTIAL, Value= 0.0 0.0 0.0  
   $FIX_PRESSURE: ON, ON_NODE = 1077, VALUE = 1.0 
   $CODES, NODES, PRESSURE
   $END_CODES
  END_PARAMETERS
$
  CODES, NODES 
    INCLUDE %s_init_codes.alya
  END_CODES
$
END_BOUNDARY_CONDITIONS  
$------------------------------------------------------------
"""

NSI_DOM = """$$Create by J. MIGUEL ZAVALA-AKE!
$----------------------------------------- mesh dims
DIMENSIONS
  NODAL_POINTS=       %d
  ELEMENTS=           %d
  SPACE_DIMENSIONS=   %d
  BOUNDARIES=         %d
  TYPES_OF_ELEMS=     %s
END_DIMENSIONS

$----------------------------------------- gauss points
STRATEGY
  INTEGRATION_RULE:          OPEN
  DOMAIN_INTEGRATION_POINTS: 0
END_STRATEGY

$----------------------------------------------- mesh files
GEOMETRY
  GROUPS=100
  INCLUDE %s_TYPES.alya
  INCLUDE %s_ELEMENTS.alya
  INCLUDE %s_COORDINATES.alya
  INCLUDE %s_BOUNDARIES.alya
$
$  FIELDS, NUMBE=1
$    FIELD=1, DIMENSION=1, NODE
$      INCLUDE  %s_CHARACTERISTICS.alya
$    END_FIELD
$  END_FIELDS
$
END_GEOMETRY

$--------------------------------------------------
SETS
$
 BOUNDARIES
  $INCLUDE %s_ON_BOUNDARIES.alya
 END_BOUNDARIES
$
 ELEME, ALL= 69
  $INCLUDE XXXX.alya 
 END_ELEME
$
END_SETS

$---------------------------------------------------
BOUNDARY_CONDITIONS, EXTRAPOLATE
  INCLUDE %s_ON_BOUNDARIES.alya
$
$  VALUES, FUNCTION = -1, DIMENSION = -3 
$    INCLUDE FIELD_3D.alya
$  END_VALUES
$
END_BOUNDARY_CONDITIONS
$---------------------------------------------------
"""


DIMENSIONS = """$$Create by J. MIGUEL ZAVALA-AKE!
$----------------------------------------- mesh dims
DIMENSIONS
  NODAL_POINTS=       %d
  ELEMENTS=           %d
  SPACE_DIMENSIONS=   %d
  BOUNDARIES=         %d
$ NODES=              %d
  TYPES_OF_ELEMS=     %s 
$ BOUNDARIES=         %d 
$ PERIODIC_NODES=     %d 
END_DIMENSIONS
"""

STRATEGY = """
$----------------------------------------- gauss points
STRATEGY
  DOMAIN_INTEGRATION_POINTS:  1
$ INTEGRATION_RULE:           OPEN
$ PERIODICITY:                MATRIX
END_STRATEGY
"""

GEOMETRY = """
$----------------------------------------------- mesh files
$GEOMETRY, GID, WALL_DISTANCE= 0.0, ROUGHNESS=0.0
 GEOMETRY
 GROUPS=100
 $
 $ FIELDS, NUMBER = %d
 $%s END_FIELDS
"""

FIELD ="""    FIELD= %d, DIMENSION %d, NODES
      INCLUDE %s
    END_FIELD\n"""

RUN_DATA = """$$Create by J. MIGUEL ZAVALA-AKE
$------------------------------------------------------- xxx.dat 
RUN_DATA
  ALYA:                   %s 
  OUTPUT:                 visit
$ LIVE_INFO:              Screen
$ OUTPUT_FORMAT: HDF5
END_RUN_DATA

$-------------------------------------------------------------------
PROBLEM_DATA
  MAXIMUM_NUMBER_GLOBAL:  1
  NUMBER_OF_STEPS:        %d 
$ TIME_STEP_SIZE=         0.1 
$ TIME_INTERVAL:          0, 0.7
$
  NASTAL_PROBLEM:         ON
  END_NASTAL
$
  PARALL_SERVICE:         ON
$   OUTPUT_FILE:          no
$   POSTPROCESS:          MASTER
$   PARTITION:            FACES
  END_PARALL
$
$  HDFPOS_SERVICE:         ON
$  END_HDFPOS_SERVICE
$
END_PROBLEM_DATA
$-------------------------------------------------------------------
"""

KER_DAT = """$$Create by J. MIGUEL ZAVALA-AKE 
$-------------------------------------------------------xxx.ker.dat
PHYSICAL_PROBLEM
$
$  LOWMACH:
$  MATERIAL = 1
$    DENSITY       = LOW_MACH
$    VISCOSITY     = CONSTANT, PARAMETER = 1.8e-5
$    SPECIFIC_HEAT = CONSTANT, PARAMETER = 1300.6
$    CONDUCTIVITY  = CONSTANT, PARAMETER = 0.012
$   $TURBULENT_VISCOSITY = SMAGO, PARAMETER = 0.01
$    POROSITY      = CONSTANT, PARAMETER = 0.0
$  END_MATERIAL
$
$  COUPLING
$  END_COUPLING
$
END_PHYSICAL_PROBLEM

NUMERICAL_TREATMENT
$
  MESH
   $DIVISION=1
  END_MESH
$
  ELSEST: ON
    STRATEGY: BIN
    NUMBER: 50,50,50
  END_ELSEST
$  
  $TIME_FUNCTIONS
  $  FUNCTION=XXX, DISCRETE, NUMBER= 4
  $    INCLUDE inflow-velocity.txt
  $  END_FUNCTION
  $END_TIME_FUNCTIONS
$
END_NUMERICAL_TREATMENT

OUTPUT_&_POST_PROCESS
  ON_LASTMESH
  STEPS = %d
$
 $WITNESS_POINTS, NUMBER=N?
   $X0 Y0 Z0?
 $END_WITNESS_POINTS
$
END_OUTPUT_&_POST_PROCESS
$-------------------------------------------------------------------
"""

POST_PROS = """$$Create by J. MIGUEL ZAVALA-AKE
$------------------------------------------------- xxx.post.alyadat
DATA
  FORMAT:                   ensight
  MARK_ELEMENTS:            type
  ELIMINATE_BOUNDARY_NODES: yes
  MULTIPLE_FILE:            OFF
  BOUNDARY:                 ON
  SUBDOMAINS, ALL
  END_SUBDOMAINS
 $SUBDOMAINS
 $ ID 
 $END_SUBDOMAINS
 $ONLY_STEP
 $ ID
 $END_ONLY_STEP
END_DATA
$-------------------------------------------------------------------
"""

NSA_DAT = """$$Create by J. MIGUEL ZAVALA-AKE
$------------------------------------------------- xxx.nsa.dat 
PHYSICAL_PROBLEM
  PROBLEM_DEFINTION  
    FORCED_REGIME:     COMPRESSIBLE
    VISCOSITY_TERMS:   ON
    INITIAL_CONDS:     REFERENCE_VALUES
  END_PROBLEM_DEFINITION  
  PROPERTIES
%s  END_PROPERTIES  
END_PHYSICAL_PROBLEM
  
$------------------------------------------------------------
NUMERICAL_TREATMENT    
  SAFETY_FACTOR:    0.4
  STABILIZATION:    MULTISCALE
 $STABILIZATION:    MULTISCALE, TRACKING
 $SHOCK_CAPTURING:  ON, VALUE= 0.7
END_NUMERICAL_TREATMENT  

$------------------------------------------------------------
OUTPUT_&_POST-PROCESS  
  POSTPROCESS XVELO,       STEPS=10
  POSTPROCESS YVELO,       STEPS=10
  POSTPROCESS ZVELO,       STEPS=10
  POSTPROCESS VELOCITY,    STEPS=10
  POSTPROCESS PRESSURE,    STEPS=10
  POSTPROCESS TEMPERATURE, STEPS=10
  POSTPROCESS DENSITY,     STEPS=10
 $POSTPROCESS MACH_NUMBER, STEPS=10
 $POSTPROCESS ENERGY,      STEPS=10
$
 $WITNESS_POINTS
 $ VELOX
 $ VELOY
 $ VELOZ
 $ PRESSURE
 $END_WITNESS_POINTS
$
END_OUTPUT_&_POST_PROCESS  

$------------------------------------------------------------
BOUNDARY_CONDITIONS, CONTINUE
  $PARAMETERS
  $  NONCONSTANT
  $END_PARAMETERS
  $ , TEMPORAL_FUNCTIONS=XXX
  CODES, NODES
    INCLUDE %s  
  END_CODES
END_BOUNDARY_CONDITIONS  
$------------------------------------------------------------
"""

RUN = """
Alyaclean.x 
Alya.x %s 
Alya2pos.x %s 
"""

#=============================================================================================#
import os 
import sys 
import numpy as np 

class Alya_field:
  def __init__(self):
    self.NODAL_POINTS     = -1
    self.ELEMENTS         = -2
    self.SPACE_DIMENSIONS = -3
    self.BOUNDARIES       = -4 
    self.NODES            = -5 
    self.TYPES_OF_ELEMS   = -6 
    self.BOUNDARIES       = -7
    self.PERIODIC_NODES   = -8 

    self.DOMAIN_INTEGRATION_POINTS = -1

    self.NUMBER_OF_STEPS = -1 

    self.INIT_DIMS  = {} 
    self.INIT_DIMS["VELOCITY"]     = 3
    self.INIT_DIMS["TEMPERATURE"]  = 1
    self.INIT_DIMS["PRESSURE"]     = 1
    self.INIT_DIMS["DENSITY"]      = 1
    self.INIT_DIMS["ENERGY"]       = 1

    self.INIT_FILES = {} 
    self.INIT_FILES["VELOCITY"]    = False 
    self.INIT_FILES["TEMPERATURE"] = False 
    self.INIT_FILES["PRESSURE"]    = False 
    self.INIT_FILES["DENSITY"]     = False 
    self.INIT_FILES["ENERGY"]      = False 

    self.INIT_VALS = {} 
    self.INIT_VALS["VELOCITY"]    = False 
    self.INIT_VALS["TEMPERATURE"] = False 
    self.INIT_VALS["PRESSURE"]    = False 
    self.INIT_VALS["DENSITY"]     = False 
    self.INIT_VALS["ENERGY"]      = False 

    self.GRUPS_NAMES = [] 
    self.GRUPS_SIZES = []

    print "+Alya_field"


  def set_prop_vals(self, val, prop_name):
    if(not self.INIT_VALS.get(prop_name, True)): 
      self.INIT_VALS[prop_name] = val 
    else: 
      print "ERROR: %s not found!! \n\n" % prop_name 
      sys.exit()

    if(prop_name=="VELOCITY"): 
      val = [str(algo) for algo in val]
      self.INIT_VALS[prop_name] = " ".join(val)      


  def set_prop_vec(self, prop_vec, file_name, prop_name, dim=3):
    file_name  = os.path.splitext(file_name)[0] 
    file_name += "_"+prop_name+".alya"

    self.INIT_FILES[prop_name] = file_name 
    #self.INIT_DIMS.append( dim ) 

    file_out = open(file_name, "w")
    for j in range(self.NODAL_POINTS):
      prop_vec_j = prop_vec[j]  
      print >> file_out, j+1, 
      for i in range(dim): print >> file_out, prop_vec_j[i],
      print >> file_out, ""
    print "|_Alya_field.set_prop_vec \'%s\'" % (file_name)
    

  def set_prop(self, prop, file_name, prop_name, dim=1):
    file_name  = os.path.splitext(file_name)[0] 
    file_name += "_"+prop_name+".alya"
    
    self.INIT_FILES[prop_name] = file_name
    #self.INIT_DIMS.append( dim ) 

    file_out = open(file_name, "w")
    for j in range(self.NODAL_POINTS):
      print >> file_out, j+1, 
      print >> file_out, prop[j]  
    print "|_Alya_field.set_prop \'%s\'" % (file_name)


  def write(self, file_name):
    print "|_Alya_field.ELEMENTS:", self.ELEMENTS
    print "|_Alya_field.NODAL_POINTS:", self.NODAL_POINTS 
    print "|_Alya_field.SPACE_DIMENSIONS:", self.SPACE_DIMENSIONS

    self.file_name = os.path.splitext(file_name)[0]

    props01 = (
    self.NODAL_POINTS,  
    self.ELEMENTS,        
    self.SPACE_DIMENSIONS,
    self.BOUNDARIES,   
    self.NODES,  
    self.TYPES_OF_ELEMS, 
    self.BOUNDARIES,  
    self.PERIODIC_NODES)

    file_names = [] 
    file_names.append(self.file_name+"_TYPES.alya") 
    file_names.append(self.file_name+"_ELEMENTS.alya") 
    file_names.append(self.file_name+"_COORDINATES.alya") 
    file_names.append(self.file_name+"_BOUNDARIES.alya") 
    file_names.append(self.file_name+"_ON_BOUNDARIES.alya") 

    COMMENT      = " "
    DOM_FIELDS   = ""
    NSA_FIELDS01 = "" 
    NSA_FIELDS02 = "" 
    nsa_i = 0 
    for (key, val) in self.INIT_FILES.items():
      if(val==False): 
        val = self.INIT_VALS.get(key, False)
        if(val==False): key = "$"+key
        else:           key = " "+key
        NSA_FIELDS01 += "   %s = %s \n" % (key, str(val)) 
      else:
        DOM_FIELDS   += FIELD %(nsa_i+1, self.INIT_DIMS[key], val)
        NSA_FIELDS02 += "      %s FIELD = %d \n" % (key, nsa_i+1) 
        nsa_i += 1
    if(nsa_i==0): COMMENT = "$"


    NSA_FIELDS  = "    FLUID     = MKS_AIR\n"
    NSA_FIELDS += "   $VISCOSITY = 1.75e-5 \n"
    NSA_FIELDS += "   $LAWVI     = SUTHE\n"
    NSA_FIELDS += NSA_FIELDS01
    NSA_FIELDS += "   %sINITIAL_FIELDS\n" % COMMENT
    NSA_FIELDS += NSA_FIELDS02
    NSA_FIELDS += "   %sEND_INITIAL_FIELDS\n" % COMMENT 

    file_out = open(self.file_name+"_dom.alya", "w")
    print >> file_out, DIMENSIONS % props01 
    print >> file_out, STRATEGY 
    print >> file_out, "$----------------------------------------------- mesh files"
    print >> file_out, "$GEOMETRY, GID, WALL_DISTANCE= 0.0, ROUGHNESS=0.0"
    print >> file_out, "GEOMETRY"
    print >> file_out, "  GROUPS=100"
    #
    print >> file_out, " %sFIELDS, NUMBER = %d" % (COMMENT, nsa_i) 
    print >> file_out, DOM_FIELDS
    print >> file_out, " %sEND_FIELDS" % COMMENT 
    print >> file_out, "$"       
    #
    num_file = 0 
    for name in file_names[:-1]:  
      print >> file_out, "  INCLUDE "+name
      num_file += 1    
    #
    print >> file_out, "$" 
    print >> file_out, "$ PERIODIC_NODES"
    print >> file_out, "$   INCLUDE %s" % self.file_name+"_PERIODIC.alya" 
    print >> file_out, "$ END_PERIODIC_NODES" 
    print >> file_out, "$" 
    print >> file_out, "END_GEOMETRY\n"
    print >> file_out,"$--------------------------------------------------"
    print >> file_out, "SETS"
    print >> file_out, "$ BOUNDARIES"
    print >> file_out, "$  INCLUDE " +file_names[num_file]  
    print >> file_out, "$ END_BOUNDARIES"
    print >> file_out, "END_SETS\n"
    print >> file_out, "$------------------------------"
    print >> file_out, "BOUNDARY_CONDITIONS, EXTRAPOLATE"
    print >> file_out, "    INCLUDE "+file_names[num_file] 
    print >> file_out, "END_BOUNDARY_CONDITIONS"
    print >> file_out,"$--------------------------------------------------\n"
    file_out.close()

    file_out = open(self.file_name+".nsa.dat", "w")
    print >> file_out, NSA_DAT % (NSA_FIELDS, self.file_name+"_init_codes.alya")
    file_out.close()

    """
    if(os.path.exists(self.file_name+".nsa.dat")): 
      os.system("cp %s %s_old" % (self.file_name+".nsa.dat", self.file_name+".nsa.dat")) 
      os.system("cp %s %s_old" % (self.file_name+".dom.dat", self.file_name+".dom.dat")) 

    file_out = open(self.file_name+".dat", "w")
    print >> file_out, RUN_DATA % (self.file_name, self.NUMBER_OF_STEPS) 
    file_out.close()

    file_out = open(self.file_name+".ker.dat", "w")
    print >> file_out, KER_DAT % (int(0.01*self.NUMBER_OF_STEPS)) 
    file_out.close()

    file_out = open(self.file_name+".post.alyadat", "w")
    print >> file_out, POST_PROS 
    file_out.close()
    
    file_out = open(self.file_name+".run", "w")
    print >> file_out, RUN % (self.file_name, self.file_name)
    file_out.close()

    os.system( "chmod +x %s" % (self.file_name+".run") )
    os.system( "cp %s_dom.alya %s.dom.dat"% (self.file_name, self.file_name) )
    """

    self.__to_nastin(self.file_name) 

    print "|_Alya_field.write \'%s\'" % (self.file_name+"*.alya/.dat")


  def set_initial_conditions_codes(self, physics_codes, keys_dic, vals_dic):
    file_name = self.file_name+"_init_codes.alya"

    file_out = open(file_name, "w") 
    #print>> file_out, "$$$ V.N|R.P|T.E"
    #for (key,vals) in physics_codes.items(): 
    #  key = key.replace('\"', '')   
    #  idx = vals 
    #  print>> file_out, "$$ \'%s\' Id:%d" % (key, idx) 
    #  print>> file_out, idx, 
    #  
    #  keys = keys_dic.get(key,'')  
    #  vals = vals_dic.get(key, []) 
    #  self.__set_initial_conditions_codes__(keys, vals, file_out)

    text = ""
    for key, val in self.GRUPS_NAMES.items(): 
      text += "$> \'%s\', id:%d, n_nodes:%d \n"% (key, val, self.GRUPS_SIZES.get(key, "??"))
      if(val>0): text += "\t%d 000 0.0 0.0 0.0 \n" % val 
    #print text  
    print>> file_out, text  
    print "|_Alya_field.set_initial_conditions_codes \'%s\'" % (file_name)


  def __set_initial_conditions_codes__(self, keys, vals, file_out):
    vals_matrix = -np.ones((3,5))

    bool_matrix = np.ones((3,5),bool)
    bool_matrix[0,:] = False
    bool_matrix[1,:] = False
    bool_matrix[2,:] = False

    init_codes = {}    
    init_codes["V"]  = (1,range(3))
    init_codes["v1"] = (1,0)
    init_codes["v2"] = (1,1)
    init_codes["v3"] = (1,2)
    init_codes["N"]  = (2,range(3))
    init_codes["n"]  = (2,0)
    init_codes["t1"] = (2,1)
    init_codes["t2"] = (2,2)
    init_codes["R"]  = (1,3)
    init_codes["P"]  = (2,3)
    init_codes["T"]  = (1,4)
    init_codes["E"]  = (2,4)

    for key, val in zip(keys, vals): 
      bool_matrix[ init_codes[key] ] = True
      vals_matrix[ init_codes[key] ] = val    
    bool_matrix[0,:] = bool_matrix[1,:]==bool_matrix[2,:]

    I = np.zeros((3,5),int)
    for i in range(3): I[i,:] = i 
    
    J = np.zeros((3,5),int)
    for j in range(5): J[:,j] = j 

    rows = np.where(bool_matrix==True, I, -1) 
    cols = np.where(bool_matrix==True, J, -1) 

    init_codes = -np.ones(5) 
    init_vals  = -np.ones(5)
    for row, col, val, ok in zip(rows.ravel(), cols.ravel(), vals_matrix.ravel(), bool_matrix.ravel()):  
      if ok:
        init_codes[col] = row
        init_vals[ col] = val

    print>> file_out,  "".join(str(int(idx)) for idx in init_codes),  
    print>> file_out, " ".join(str(idx) for idx in init_vals), 
    print>> file_out, "$", "".join(keys) 
    print>> file_out 
    
    return init_codes, init_vals


  def __to_nastin(self, file_name):
    file_name = os.path.splitext(file_name)[0]

    props01  = [ 
    self.NODAL_POINTS,  
    self.ELEMENTS,        
    self.SPACE_DIMENSIONS,
    self.BOUNDARIES,   
    self.TYPES_OF_ELEMS]
 
    props01 += ["../"+file_name]*7 

    NSI = file_name.upper()
    if(not os.path.exists(NSI)): os.mkdir(NSI)
    print "|_\'%s\' created!" % (NSI)
    #os.system("mv *nsi* %s" % NSI)

    file_out = open(NSI+"/"+file_name+".dom.dat", "w")
    print >> file_out, NSI_DOM % tuple(props01) 
    file_out.close()

    file_out = open(NSI+"/"+file_name+".nsi.dat", "w")
    print >> file_out, NSI_NSIDAT % ("../"+file_name) 
    file_out.close()

    file_out = open(NSI+"/"+file_name+".ker.dat", "w")
    print >> file_out, KER_DAT % (int(0.01*self.NUMBER_OF_STEPS)) 
    file_out.close()

    file_out = open(NSI+"/"+file_name+".post.alyadat", "w")
    print >> file_out, POST_PROS 
    file_out.close()

    file_out = open(NSI+"/"+file_name+".dat", "w")
    print >> file_out, NSI_DAT % (file_name, self.NUMBER_OF_STEPS, -1) 
    file_out.close()

    #file_out = open(file_name+"_nsi.run", "w")
    #print >> file_out, RUN % (file_name+"_nsi", file_name+"_nsi")
    #file_out.close()
    
    #os.system( "chmod +x %s" % (file_name+"_nsi.run") )

    #NSI = "./NSI"
    #NSI = file_name.upper()  
    #if(not os.path.exists(NSI)): os.mkdir(NSI)
    #print "|_\'%s\' created!" % (NSI)
    #os.system("mv *nsi* %s" % NSI)


#=============================================================================# 

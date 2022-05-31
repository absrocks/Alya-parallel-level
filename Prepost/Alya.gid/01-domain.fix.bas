*set var inode=0
*Set Cond ALYA_CODES_node_volume  *nodes
*Add Cond ALYA_CODES_node_surface *nodes
*Add Cond ALYA_CODES_node_line    *nodes
*Add Cond ALYA_CODES_node_point   *nodes
*loop nodes *OnlyInCond
*set var inode=inode+1
*end
*if(inode!=0) 
  ON_NODES
*loop nodes *OnlyInCond
*format " %i %i "
  *NodesNum *cond(1,int) 
*end
  END_ON_NODES
*endif
*set var iboun=0
*if(ndime==2)
*set Cond ALYA_CODES_boundary_2D *elems
*else
*set Cond ALYA_CODES_boundary_3D *elems
*endif
*loop elems *onlyincond *canrepeat
*set var iboun=iboun+1
*end elems
*if(iboun!=0) 
  ON_BOUNDARIES, UNKNOWN
*loop elems *onlyincond *canrepeat
*set var iboun=iboun+1
*format " %i %i %i %i %i %i %i %i %i %i %i %i "
*iboun *ElemsNnodeFace *globalnodes *cond(1,int) 
*end elems
  END_ON_BOUNDARIES
*endif

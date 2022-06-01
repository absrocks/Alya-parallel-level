#
#  alya-deposition
#
alya-deposition create the paraview output with all the depositions of the particules ==>  time, particule_id, particule_type, X, Y, Z, for all the particules touching the boundary

example: alya-deposition.x fensap

#
#  alya-particule
#
alya-particule create the paraview output to visualize the particules, there is an argumente to postprocess not all the particules, it is a module
==>time,part_id,xcoord,ycoord,zcoord,xvelo,yvelo,zvelo,part_type,subdom,family(part_id)

example : alya-particule.x fensap 10
if you want divide by 10 the number of particule post processed !
 

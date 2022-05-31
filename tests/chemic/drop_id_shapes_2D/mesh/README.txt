Modifying the domain and the mesh:

    1. change to mesh directory
    2. edit mesh3d.geo:
       1. gx, gy, gz determine the size of the box in x,y,z.
       2. nGapx, nGapy, nGapz determine the number of points on each edge
    (this is one more than the number of elements on that edge)
    3. execute: "source batchMesh.sh"
       1. This runs gmsh, converts it to Alya format and generates the initial
    condition.
       2. After this it copies the files to the case folder, so you don't have
    to do that manually.
       3. I think you will have to specify the full path to gmsh2alya.pl, so the
    4th line has to be: YourPathToVapor/Utils/user/gmsh2alya.pl $name -bcs=boundaries

"batchMesh.sh" executes: (info to verify)
1. gmsh, that uses mesh3d.geo and builds mesh3d.msh
2. gmsh2alya.pl, that builds mesh3d.dims.dat, mesh3d.geo.dat, mesh3d.fix.bou 
3. getCoordinates.py, that uses mesh3d.geo.dat and builds mesh3d.coord
4. initailCondition.py, that uses mesh3d.coord and builds [field_name].alya
(e.g.: VELOC.alya)
5. getElementSet.py, that usesmesh3d.dims.dat and builds mesh3d.set.elem

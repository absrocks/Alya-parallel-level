gmsh -2 mesh2d.geo
gmsh2alya.pl mesh2d -bcs=boundaries
sed -i '1d; $d' mesh2d.fix.bou 

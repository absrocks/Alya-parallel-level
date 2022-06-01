gmsh -2 -format msh2  mesh1d.geo
gmsh2alya.pl mesh1d -bcs=boundaries
sed -i -e '1d; $d' mesh1d.fix.bou 

gmsh -2 mesh2d.geo
gmsh2alya.pl mesh2d -bcs=boundaries
sed -i -e '1d; $d' mesh2d.fix.bou 
python geo2fields-opp_diff.py mesh2d.geo.dat

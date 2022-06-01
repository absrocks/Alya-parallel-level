name=mesh2d
gmsh -2 -format msh2 $name.geo
gmsh2alya.pl $name -bcs=boundaries
sed -i -e '1d; $d' $name.fix.bou 

python2 geo2fields-opp_diff.py $name.geo.dat



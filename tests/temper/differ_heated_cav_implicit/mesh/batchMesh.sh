gmsh -2 mesh2d.geo
gmsh2alya.pl mesh2d -bcs=boundaries
sed -i '1d; $d' mesh2d.fix.bou 
mv mesh2d.geo.dat mesh2d.geo.dat_unchecked
python InviertePorJacobianoTriQuad_v4.py mesh2d.geo.dat_unchecked mesh2d.geo.dat
rm mesh2d.geo.dat_unchecked
cp mesh2d.geo.dat ..
cp mesh2d.fix.bou ..

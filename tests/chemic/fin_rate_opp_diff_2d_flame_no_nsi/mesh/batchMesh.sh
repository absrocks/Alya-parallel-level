name=mesh2d
gmsh -2 $name.geo
gmsh2alya.pl $name -bcs=boundaries
python getCoordinates.py $name
python initialCondition.py $name 2S_CM2.xml "CH4:1" "N2:0.767083, O2: 0.232917" 298.0 101325 "1,1 & 3" "2,2 & 3" 
sed -i -e '1d; $d' $name.fix.bou 

cp $name.dims.dat ..
cp $name.geo.dat ..
cp $name.fix.bou ..
rm -rf ../Fields
mkdir ../Fields
cp Fields/*.alya ../Fields/
cp fields.dat ..
cp field_size.dat ..
cp speciesBoundary.dat ..

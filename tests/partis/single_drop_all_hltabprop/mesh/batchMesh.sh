module load python/2.7.13 
name='mesh3d'
gmsh -3 $name".geo"
gmsh2alya.pl $name -bcs=boundaries
python getCoordinates.py $name
python initialCondition.py $name 
python getElementSet.py $name

sed -i -e '1d; $d' $name".fix.bou" 

cp $name".dims.dat" ..
cp $name".geo.dat" ..
cp $name".fix.bou" ..
cp $name".set.elm" ..

cp *.alya ..

name='mesh3d'
#module load python/2.7.13 
gmsh -3 $name".geo"
gmsh2alya.pl $name -bcs=boundaries
python getCoordinates.py $name
python initialCondition.py $name 2S_CM2.xml "CH4:1" "N2:0.767083, O2: 0.232917"  0.67 1500.0  101325   
python getRealBoundaries.py $name 

xmin=$(python getBoundingBox.py $name x min)
xmax=$(python getBoundingBox.py $name x max)
ymin=$(python getBoundingBox.py $name y min)
ymax=$(python getBoundingBox.py $name y max)
zmin=$(python getBoundingBox.py $name z min)
zmax=$(python getBoundingBox.py $name z max)


per=$name".per"
rm $per 
touch $per
../../../../Utils/user/periodic_extractor/build/per.x $name $xmin $xmax 0 3  #arguments: name of case, lower limit of coordinate, upper limit of coordinate, index of coordinate: x=0 y=1 z=2, dimension of problem
cat $name".x.per" >> $per
rm $name".x.per"
../../../../Utils/user/periodic_extractor/build/per.x $name $ymin $ymax 1 3
cat $name".y.per" >> $per
rm $name".y.per"
../../../../Utils/user/periodic_extractor/build/per.x $name $zmin $zmax 2 3
cat $name".z.per" >> $per
rm $name".z.per"
python ../../../../Utils/user/periodic_extractor/uniquePeriodic.py $per


echo ""
echo "Number of periodic nodes:"
PER_NODES=$(wc -l < *per)
echo $PER_NODES

echo  "  PERIODIC_NODES="$PER_NODES >> $name".dims.dat"

sed -i -e '1d; $d' $name".fix.bou" 


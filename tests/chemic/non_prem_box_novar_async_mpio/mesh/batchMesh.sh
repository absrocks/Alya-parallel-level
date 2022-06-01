name='mesh3d'
#module load python/2.7.13 
gmsh -3 $name".geo"
gmsh2alya.pl $name -bcs=boundaries
python getCoordinates.py $name
python initialCondition.py $name 
python getElementSet.py $name 
#python getRealBoundaries.py $name 1 2 3

#xmin=$(python getBoundingBox.py $name x min)
#xmax=$(python getBoundingBox.py $name x max)
#ymin=$(python getBoundingBox.py $name y min)
#ymax=$(python getBoundingBox.py $name y max)
#zmin=$(python getBoundingBox.py $name z min)
#zmax=$(python getBoundingBox.py $name z max)


#per=$name".per"
#rm $per 
#touch $per
#./periodic_extractor/build/per.x $name $xmin $xmax 0 3  #arguments: name of case, lower limit of coordinate, upper limit of coordinate, index of coordinate: x=0 y=1 z=2, dimension of problem
#cat $name".x.per" >> $per
#rm $name".x.per"
#./periodic_extractor/build/per.x $name $ymin $ymax 1 3
#cat $name".y.per" >> $per
#rm $name".y.per"
#./periodic_extractor/build/per.x $name $zmin $zmax 2 3
#cat $name".z.per" >> $per
#rm $name".z.per"
#python uniquePeriodic.py $per


#echo ""
#echo "Number of periodic nodes:"
#wc -l *per

sed -i -e '1d; $d' $name".fix.bou" 

cp $name".dims.dat" ..
cp $name".geo.dat" ..
cp $name".fix.bou" ..
cp $name".set.elm" ..
#cp $name".fix.bou.real" ..
#cp $name".per" ..
cp *.alya ..

module load python/2.7.13 
name='mesh3d'
gmsh -3 -format msh2 $name".geo"
gmsh2alya.pl $name -bcs=boundaries
python getElementSet.py $name

sed -i -e '1d; $d' $name".fix.bou" 


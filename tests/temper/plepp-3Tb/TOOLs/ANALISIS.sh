
echo "RUNNIG:" $(date +'%Y-%m-%d') $(date +'%T')
echo

printf -v nameA "%s_witness1.dat" $1
printf -v nameB "%s_witness1.dat" $2 

echo $nameA 
grep "CURRENT TIME t" $1.log | tail 

echo $nameB 
grep "CURRENT TIME t" $2.log | tail                    

#Plot_wits.py -F $1.tem.wit -W TEMPE
##Plot_multignuplot.py -F $nameA -X  '( ($1>=0.0)?($1):(1/0) )' -Y  '($2)' -Rx ':' -Ry ':'
#
#gnuplot ../TOUSE/plotBNDRYLAYER01.gp 



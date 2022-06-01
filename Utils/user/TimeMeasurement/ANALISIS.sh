
#grep "CURRENT TIME STEP dt" vortex2D.log | tail
#grep "CURRENT TIME t"       vortex2D.log | tail
echo 
echo 

module load python
#python /home/bsc21/bsc21704/z2017/RUNNER/HPC02/DLR2D_02/CHT_PLEPP04/checkTimeMeasuraments01.py -F "."

#cp ../plotCoupled.gp .
#gnuplot /home/bsc21/bsc21704/z2017/RUNNER/PLEPPs/plotCoupled02.gp
#gnuplot /home/bsc21/bsc21704/z2017/RUNNER/PLEPPs/plotCoupled01.gp




#grep "TEMPER MODULE"  TAU12_03_F0512x*/*.log | grep %
#python /home/bsc21/bsc21704/z2017/RUNNER/HPC02/DLR2D_02/CHT_PLEPP04/checkTimeMeasuraments01.py -F "./TAU12_03*"
#tar -cf TAU12_03.tar TAU12_03_*/dirichl* TAU12_03_*/neumann* TAU12_03_*/NEUMANN* TAU12_03_*/DIRICHL* TAU12_03_*/*.gp quartile50th.* 

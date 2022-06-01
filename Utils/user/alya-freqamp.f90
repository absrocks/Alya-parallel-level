program freqamp
  implicit none
  real(8)        :: t1,t2,f1(20),f2(20),tinit,cutim
  real(8)        :: fsign,a,b,ta,tb,fmini,fmaxi
  integer(4)     :: ii,icolu
  integer(4)     :: one4=1
  character(100) :: filen

  call GETARG(one4,filen)
  
  open(unit=10,file=trim(filen),status='old',err=1)

  write(6,*) ' '
  write(6,20) 'alya-freqamp |--'
  write(6,*) ' ' 
  write(6,20) 'Compute frequence and amplitude'
  write(6,*) ' '

  write(6,30) 'Initial time = '
  read(5,*) tinit
  write(6,30) 'Column       = '
  read(5,*) icolu

  write(6,*) ' '

  icolu=icolu-1
  read(10,*,err=3)
  read(10,*,err=3)
  read(10,*,err=3)
  t2=-1.0e6
  do while(t2<tinit)
     read(10,*,err=3) t2,(f2(ii),ii=1,icolu)
  end do
  fsign=1.0
  do while(fsign>0.0)
     f1(icolu)=f2(icolu)
     t1=t2
     read(10,*,err=3) t2,(f2(ii),ii=1,icolu)
     fsign=f2(icolu)*f1(icolu)
  end do
  a=(f1(icolu)-f2(icolu))/(t1-t2)
  b=f1(icolu)-a*t1
  ta=-b/a

  fmini=1e6
  fmaxi=-1e6
  fsign=1.0
  do while(fsign>0.0)
     f1(icolu)=f2(icolu)  
     t1=t2
     read(10,*,err=3) t2,(f2(ii),ii=1,icolu)
     fsign=f2(icolu)*f1(icolu)
     if(f2(icolu)>fmaxi) fmaxi=f2(icolu)
     if(f2(icolu)<fmini) fmini=f2(icolu)
  end do

 
  fsign=1.0
  do while(fsign>0.0)
     f1(icolu)=f2(icolu)  
     t1=t2
     read(10,*,err=3) t2,(f2(ii),ii=1,icolu)
     fsign=f2(icolu)*f1(icolu)
     if(f2(icolu)>fmaxi) fmaxi=f2(icolu)
     if(f2(icolu)<fmini) fmini=f2(icolu)
  end do
  a=(f1(icolu)-f2(icolu))/(t1-t2)
  b=f1(icolu)-a*t1
  tb=-b/a
  write(6,40) 'Initial time = ',ta
  write(6,40) 'Final   time = ',tb
  write(6,40) 'Period       = ',tb-ta
  write(6,40) 'Frequence    = ',1.0/(tb-ta)
  write(6,40) 'Amplitude    = ',fmaxi-fmini

  write(6,*)
  write(6,20) 'Bye'
  write(6,*)
  write(6,*)

  stop

1 write(6,20) '--| Error: could not open file'
  stop
2 write(6,20) '--| Error: wrong initial time'
  stop
3 write(6,20) '--| Error: error while reading file'
  stop
20 format('--| ',a)
30 format('--| ',a,$)
40 format('--| ',a,e12.6)

end program freqamp

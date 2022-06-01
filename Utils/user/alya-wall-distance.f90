program alya_wall_distance
  implicit none
  real(8)    :: L,U,mu,rho,y1,y2,y3,y4,yplus,Re,d
  integer(4) :: i

  rho=1.2047
  mu=1.8205e-05 
  L=1.0
  U=1.0

  write(*,'(a)') ' '
  write(*,'(a)') '--| alya-wall-distance |-- '
  write(*,'(a)') ' '
  write(*,'(a)') '--| Compute the distance to the wall y for y+=1,10,30,50'
  write(*,'(a)') '--| and boundary layer thickness d'
  write(*,'(a)') '--| '
  write(*,'(a)') '--| 1. Channel flow. 2 Flat plate.'
  write(*,'(a)') '--| '
  write(*,'(a)') '--| 1. y = 5.323*H*y+*Re^{-7/8}'
  write(*,'(a)') '--| '
  write(*,'(a)') '--| 2. y = 5.199*L*y+*Re^{-9/10}'
  write(*,'(a)') '--|    d = 0.14*L/Re^{1/7}'
  write(*,'(a)') ' '

  write(*,'(a,$)') '--| Choose: 1. Channel flow. 2 Flat plate ? '
  read(*,*,err=10) i
  write(*,'(a,$)') '--| Characteristic length ..... L=   '
  read(*,*,err=10) L
  write(*,'(a,$)') '--| Characteristic velocity ... U=   '
  read(*,*,err=10) U
  write(*,'(a,$)') '--| Density ................... rho= '
  read(*,*,err=10) rho
  write(*,'(a,$)') '--| Viscosity ................. mu=  '
  read(*,*,err=10) mu

  write(*,'(a)') '--|   '  
  Re=rho*U*L/mu
  yplus=1.0_8
  if(i==1) then
     y1=5.323_8*L*yplus*Re**(-7.0_8/8.0_8)
  else
     y1=5.199_8*L*yplus*Re**(-9.0_8/10.0_8)
  end if
  yplus=10.0_8
  if(i==1) then
     y4=5.323_8*L*yplus*Re**(-7.0_8/8.0_8)
  else
     y4=5.199_8*L*yplus*Re**(-9.0_8/10.0_8)
  end if
  yplus=30.0_8
  if(i==1) then
     y2=5.323_8*L*yplus*Re**(-7.0_8/8.0_8)
  else
     y2=5.199_8*L*yplus*Re**(-9.0_8/10.0_8)
  end if
  yplus=50.0_8
  if(i==1) then
     y3=5.323_8*L*yplus*Re**(-7.0_8/8.0_8)
  else
     y3=5.199_8*L*yplus*Re**(-9.0_8/10.0_8)
  end if
  write(*,'(a,e12.6)') '--| Reynolds number .............  Re=  ', Re
  write(*,'(a,e12.6)') '--| Wall distance (y+=1) ........  d1=  ', y1
  write(*,'(a,e12.6)') '--| Wall distance (y+=10) .......  d2=  ', y4
  write(*,'(a,e12.6)') '--| Wall distance (y+=30) .......  d3=  ', y2
  write(*,'(a,e12.6)') '--| Wall distance (y+=50) .......  d4=  ', y3
  if(i==2) then
     d=0.14_8*L/Re**(1.0_8/7.0_8)
     write(*,'(a,e12.6)') '--| Boundary layer thickness ....  d =  ', d 
     !d=0.37_8*L/Re**(1.0_8/5.0_8)
     !write(*,'(a,e12.6)') '--| Boundary layer thickness ....  d =  ', d 
  end if
  write(*,'(a)') ' '
  write(*,'(a)') '--| Bye.'
  write(*,'(a)') ' '
 
  stop
10 write(*,'(a)') 'Wrong data'
  
end program alya_wall_distance

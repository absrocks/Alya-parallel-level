subroutine ibm_sediib(coord,coor1,coor2,ndime,dista,proje,ifoun)
  !-----------------------------------------------------------------------
  ! NAME
  !    linedi
  ! DESCRIPTION
  !    Minimun distance between a segmnent and a point
  !    point : point
  !    coor1,coor2 : defines the segment
  !    ndime: dimension
  !    dista: distance
  !    proje: projection of the point on the segment
  !    ifoun = 1 the projection point is inside the segment
  !    ifoun = 0 the projection point is outside the segment
  ! USED BY
  !    pofadi
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only           :  ip,rp
  implicit none
  integer(ip),   intent(in)      :: ndime
  integer(ip),   intent(out)     :: ifoun
  real(rp),      intent(in)      :: coord(ndime),coor1(ndime),coor2(ndime)
  real(rp),      intent(out)     :: dista,proje(ndime)
  integer(ip)                    :: idime
  real(rp)                       :: unive(3),vecto(3)
  real(rp)                       :: dista1,dista2,invdis
  
  dista1 = 0.0_rp

  do idime = 1,ndime
     !
     ! unit vector in the segment direction
     !
     unive(idime) = coor2(idime) - coor1(idime)
     !
     ! vector from one of the segment vertices to the point
     !
     vecto(idime) = coord(idime) - coor1(idime)
     dista1       = dista1 + (coor2(idime) - coor1(idime)) * (coor2(idime) - coor1(idime))
  end do

  if( dista1 == 0.0_rp ) then
     invdis = 1.0e6_rp
  else
     invdis = 1.0_rp / sqrt(dista1)
  end if

  do idime = 1,ndime
     unive(idime) = unive(idime) * invdis
  end do
  !
  ! distance from the first segment vertex to the point in the direction of the segment
  !
  dista2 = 0.0_rp
  do idime = 1,ndime
     dista2 = dista2 + (unive(idime) * vecto(idime))
  end do
  !
  ! normalized distance
  !
  dista = dista2 * invdis
    
  ! If the projection is inside the segment   
  if( dista >= 0.0_rp .and. dista <= 1.0_rp ) then 

     ifoun = 1     
     do idime = 1,ndime           
        proje(idime) = coor1(idime) + dista2 * unive(idime)
     end do
     dista = 0.0_rp
     do idime = 1,ndime           
        dista = dista + (coord(idime)-proje(idime)) * (coord(idime)-proje(idime))
     end do
     dista = sqrt(dista)     
  else
     ifoun  = 0     
     dista1 = 0.0_rp
     dista2 = 0.0_rp
     do idime = 1,ndime
        dista1 = dista1 + (coord(idime) - coor1(idime)) * (coord(idime) - coor1(idime))
        dista2 = dista2 + (coord(idime) - coor2(idime)) * (coord(idime) - coor2(idime))
     end do
     if( dista1 <= dista2 ) then
        dista = sqrt(dista1)
        do idime=1,ndime
           proje(idime) = coor1(idime)
        end do
     else
        dista = sqrt(dista2)
        do idime = 1,ndime
           proje(idime) = coor2(idime)
        end do
     end if
  end if

end subroutine ibm_sediib

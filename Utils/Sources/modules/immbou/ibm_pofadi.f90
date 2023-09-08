subroutine ibm_pofadi(coord,bocod,norma,ndime,dista,proje)
  !-----------------------------------------------------------------------
  ! NAME
  !    ibm_pofadi
  ! DESCRIPTION
  !    Minimun distance between a point and a face
  !    point: point
  !    bocod(dimension,vertex): triangle coordinates
  !    norma: nornmalized exterior normal to the triangle
  !    dista: minimum distance between the point and the triangle
  !    proje: nearest triangle point projection  
  ! USED BY
  !    shdiib
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only           :  ip,rp
  use def_domain, only           :  mnoib
  implicit   none
  integer(ip),   intent(in)      :: ndime
  real(rp),      intent(in)      :: coord(3),norma(3)
  real(rp),      intent(in)      :: bocod(ndime,mnoib)
  real(rp),      intent(out)     :: dista,proje(3)
  integer(ip)                    :: idime,ifoun
  real(rp)                       :: pladi,plapo(3),factor
  real(rp)                       :: dist1,dist2,dist3
  real(rp)                       :: proje1(3),proje2(3),proje3(3)
  !
  ! Distance between the point and the plane formed by the face
  ! Dot product between the exterior normal and the vector build from the first triangle vertex
  ! to the point  
  !
  pladi = 0.0_rp
  do idime = 1,ndime
     pladi = pladi + norma(idime) * ( coord(idime) - bocod(idime,1) )
  end do
  !
  ! Point projection on the plane
  !
  do idime=1,ndime
     plapo(idime) = coord(idime) - pladi * norma(idime)
  end do

  if ( pladi < 0.0_rp ) then
     factor = -1.0_rp
  else
     factor =  1.0_rp
  end if

  if( ndime == 3 ) then        
     !
     ! Determine if the projection point on plane is inside the triangle
     !
     call ibm_insiib(bocod,plapo,ndime,ifoun)
     !
     ! The projection point is inside the triangle
     !
     if( ifoun == 1 ) then
        dista = pladi
        do idime = 1,ndime
           proje(idime) = plapo(idime)
        end do

     else
        
        call ibm_sediib(coord,bocod(1,1),bocod(1,1),ndime,dist1,proje1,ifoun)
        call ibm_sediib(coord,bocod(1,2),bocod(1,2),ndime,dist2,proje2,ifoun)
        call ibm_sediib(coord,bocod(1,3),bocod(1,3),ndime,dist3,proje3,ifoun)

        if( dist1 <= dist2 .and. dist1 <= dist3 ) then
           
           dista = dist1 * factor
           do idime = 1,ndime
              proje(idime) = proje1(idime)
           end do
        else if( dist2 <= dist1 .and. dist2 <= dist3 ) then
           dista = dist2 * factor
           do idime=1,ndime
              proje(idime) = proje2(idime)
           end do
        else
           dista = dist3 * factor
           do idime = 1,ndime
              proje(idime) =  proje3(idime)
           end do
        end if
     end if

  else if( ndime == 2 ) then
     
     call ibm_sediib(coord,bocod(1,1),bocod(1,2),ndime,dist1,proje1,ifoun)        
     dista = dist1 * factor
     do idime = 1,ndime
        proje(idime) = proje1(idime)
     end do
     
  end if

end subroutine ibm_pofadi



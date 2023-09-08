subroutine ibm_qualib(elcod,chang,Q,Qmax)
  !-----------------------------------------------------------------------
  !***
  ! NAME
  !    ibm_qualib
  ! DESCRIPTION
  !    This routines determine the quality of a triangle or a tetrahedron
  !    The quality is good if the value of Q is small
  ! USED BY
  !    fielib
  !***
  !-----------------------------------------------------------------------
  use def_kintyp,   only  :  ip,rp
  use def_domain,   only  :  ndime
  use def_master,   only  :  zeror
  implicit none
  integer(ip),intent(out) :: chang
  real(rp),   intent(in)  :: elcod(ndime,*)
  real(rp),   intent(out) :: Q
  real(rp),   intent(in)  :: Qmax
  integer(ip)             :: idime
  real(rp)                :: hmax,alpha,volum,perim,surfa(4),tosur
  real(rp)                :: side(3,6)
  real(rp)                :: dista(6),cropr(3)
  
  if( ndime == 2 ) then
     
     alpha = sqrt(3.0_rp)/6.0_rp
     
     do idime = 1,ndime
        side(idime,1) = elcod(idime,2) - elcod(idime,1)
        side(idime,2) = elcod(idime,3) - elcod(idime,2)
        side(idime,3) = elcod(idime,1) - elcod(idime,3)
     end do

     side(3,1) = 0.0_rp
     side(3,2) = 0.0_rp
     side(3,3) = 0.0_rp
     
     call vecnor(side(1,1),3_ip,dista(1),2_ip)
     call vecnor(side(1,2),3_ip,dista(2),2_ip)
     call vecnor(side(1,3),3_ip,dista(3),2_ip)
     !
     ! Largest side
     !        
     hmax = maxval(dista(1:3))
     !
     ! Perimeter
     !             
     perim = 0.5_rp * ( dista(1) + dista(2) + dista(3) )
     !
     ! Area
     !     
     call vecpro(side(1,1),side(1,2),cropr,3_ip)
     call vecnor(cropr,3_ip,volum,2_ip)
     volum = 0.5_rp * volum
     
     chang = 0
     if( abs(volum) <= zeror ) then 
        Q     = Qmax
     else if( volum < 0.0_rp ) then
        chang = 1
        Q     = alpha * hmax * perim /(-volum)
     else
        Q     = alpha * hmax * perim / volum
     end if

  else if( ndime == 3 ) then

     alpha = sqrt(6.0_rp)/12.0_rp
     
     do idime=1,ndime
        side(idime,1) = elcod(idime,2) - elcod(idime,1)
        side(idime,2) = elcod(idime,3) - elcod(idime,2)
        side(idime,3) = elcod(idime,1) - elcod(idime,3)
        side(idime,4) = elcod(idime,4) - elcod(idime,1)
        side(idime,5) = elcod(idime,4) - elcod(idime,2)
        side(idime,6) = elcod(idime,4) - elcod(idime,3)
     end do

     call vecnor(side(1,1),ndime,dista(1),2_ip)
     call vecnor(side(1,2),ndime,dista(2),2_ip)
     call vecnor(side(1,3),ndime,dista(3),2_ip)
     call vecnor(side(1,4),ndime,dista(4),2_ip)
     call vecnor(side(1,5),ndime,dista(5),2_ip)
     call vecnor(side(1,6),ndime,dista(6),2_ip)
     !
     ! Largest side
     !     
     hmax = maxval(dista(1:6))
     !
     ! Surface areas
     !
     call vecpro(side(1,1),side(1,2),cropr,ndime)
     call vecnor(cropr,ndime,surfa(1),2_ip)
     surfa(1) = 0.5_rp*abs(surfa(1))

     call vecpro(side(1,1),side(1,5),cropr,ndime)
     call vecnor(cropr,ndime,surfa(2),2_ip)
     surfa(2) = 0.5_rp*abs(surfa(2))

     call vecpro(side(1,2),side(1,6),cropr,ndime)
     call vecnor(cropr,ndime,surfa(3),2_ip)
     surfa(3) = 0.5_rp*abs(surfa(3))

     call vecpro(side(1,3),side(1,4),cropr,ndime)
     call vecnor(cropr,ndime,surfa(4),2_ip)
     surfa(4) = 0.5_rp*abs(surfa(4))
          
     tosur = surfa(1) + surfa(2) + surfa(3) + surfa(4)
     !
     ! Volume
     !
     call vecpro(side(1,2),side(1,3),cropr,ndime)
     do idime=1,ndime
        volum=cropr(idime)*side(idime,4)     
     end do
     
     chang = 0
     if( abs(volum) <= zeror ) then 
        Q     = Qmax
     else if ( volum < 0.0_rp ) then        
        chang = 1
        Q     = alpha * hmax * tosur /(3.0_rp*(-volum))
     else
        Q     = alpha * hmax * tosur /(3.0_rp*(volum))
     end if
     
  end if

end subroutine ibm_qualib

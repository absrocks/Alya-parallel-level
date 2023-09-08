subroutine chm_setpsi(&
     shape_chm,spher_chm,diame_chm,lawvt_chm)
  !-----------------------------------------------------------------------
  !****f* chemic/chm_setpsi
  ! NAME 
  !    chm_setpsi
  ! DESCRIPTION
  !    Calculates the particle shape factor psi depending on the velocity mode
  !       lawvt_chm = 1   GANSER      psi = sphericity
  !       lawvt_chm = 2   WILSON      psi = (b+c)/2a    a>b>c semi-axes
  !       lawvt_chm = 3   DELLINO     psi = sphericity/circularity      
  ! USES
  ! USED BY
  !    chm_reaphy
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  implicit none
  !  
  integer(ip),  intent(in)  :: lawvt_chm
  real(rp),     intent(in)  :: spher_chm,diame_chm
  real(rp),     intent(out) :: shape_chm
  !
  real(rp) :: pi,gama
  !
  !  Initializations
  !
  pi = 4.0_rp*atan(1.0_rp) 
  !  
  !   Computes psi
  !
  if ( lawvt_chm == 1 .or. lawvt_chm == 4 ) then
     shape_chm = spher_chm

  else if ( lawvt_chm == 2 ) then
     call get_gama(diame_chm,shape_chm,gama)       ! Get a/c
     if(gama.ge.1.0_rp) then                       ! oblate
        shape_chm = 0.5_rp*(1.0_rp+1.0_rp/gama) 
     else                                          ! prolate
        shape_chm = gama
     end if
 
  else if ( lawvt_chm == 3 ) then	   
     call get_gama(diame_chm,spher_chm,gama)         ! Get a/c
     if(gama.ge.1.0_rp) then                         ! oblate
        shape_chm = shape_chm
     else                                            ! prolate
        shape_chm = shape_chm/sqrt(gama)
     end if
  end if
  !
end subroutine chm_setpsi
!
!
!
subroutine get_gama(diam,sphe,gama)
  !-----------------------------------------------------------------------
  !
  !     Gets gama = a/c 
  !
  !     NOTE: In all cases it is assumed that particles fall as 
  !     prolate ellipsoids
  !            
  !            a = b < c  prolate   (gama < 1)
  !
  !    The inversion of the area of the ellipsoid is done numerically. 
  !     Area given by:
  !
  !     A = 2*pi*(a**2 + c**2*e/tan(e) )   e = acos(gama)   
  !     d = 2*c*gama**(2/3)               (prolate)
  !     
  !    NOTE: particle diameter is multiplied by a factor. It does not affect
  !          results (a/c) and is done to facilitate convergence and prevent
  !           propagation of rounding errors (e.g. for micron size particles
  !           diam of the order 1d-6 rised to 2 or 3) 
  !
  !-----------------------------------------------------------------------  
  use def_kintyp, only      :  ip,rp
  implicit none
  real(rp) :: diam,sphe,gama
  !      
  integer(ip) :: iiter,niter
  real   (rp) :: d,pi,gmin,gmax,Ao,toler,e
  real   (rp) :: Vp,Ap
  !
  !***   Initializations
  !
  d     = diam*1d3         ! see NOTE
  niter = 1000
  toler = 1d-8
  gmin  = 1d-3
  gmax  = 1.0_rp 
  !
  !***   Volume and area
  !
  pi = 4.0_rp*atan(1.0_rp)
  Vp = 4.0_rp*pi*((0.5_rp*d)**3.0_rp)/3.0_rp
  Ap = (pi**(1.0_rp/3.0_rp))*((6.0_rp*Vp)**(2.0_rp/3.0_rp))/sphe
  !
  !***   Iterates
  !
  do iiter = 1,niter
     gama = 0.5_rp*(gmin+gmax)
     e    = acos(gama)
     Ao   = 0.5_rp*pi*d*d*(gama**(-4.0_rp/3.0_rp))*(gama*gama + (e/tan(e)))
     if(Ao.lt.Ap) then
        gmax = gama
     else
        gmin = gama
     end if
     if((iiter.gt.1).and.(abs(Ao-Ap).lt.toler)) goto 10
  end do
  call runend('Subroutine get_gama: convergence not achieved')
  !
  !***  convergence
  !
10 return
end subroutine get_gama



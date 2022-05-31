subroutine nsa_shocap(ifcdr,kshoc,ndime,  &
     shfac,graun,cartd,resid,hconv,xvelo,taupa,difeq,zesho,qufac,shote)
!-----------------------------------------------------------------------
!
! 
!
!
!
!
!-----------------------------------------------------------------------------------------------------------
  use      def_kintyp
  implicit none
  integer(ip),  intent(in) :: ifcdr,kshoc,ndime
  integer(ip)           :: kfcdr,inode,ipoin,idome,igaus,itidi,idime,jdime,itime

  real(rp), intent(in)  :: graun(ndime),cartd(ndime),resid,hconv,xvelo(ndime)
  real(rp), intent(in)  :: qufac,taupa,shfac,difeq,zesho
  real(rp), intent(out) :: shote

  real(rp)              :: grau2,xvel2,grmod,remod,ficve,shpec,difau,shfau,zesh2,shmet(3,3)
  real(rp)              :: epsdc,epssu,epssl,epsde,fiso,faniso,uvel,vvel,wvel


! !$OMP


  if (kshoc == 0) return

  grau2= 0.0_rp
  xvel2= 0.0_rp
  do idime = 1,ndime                                     
     grau2 = grau2 + graun(idime)*graun(idime)            !evaluate grad modul
     xvel2 = xvel2 + xvelo(idime)*xvelo(idime)            !evaluate veloc modul
     shmet(1:ndime,idime) = 0.0_rp                        !initialize shock metrics
  end do
  
  if (grau2 < zesho) return
  if (xvel2 < zesho) return
  
  zesh2 = zesho*zesho
  grmod = sqrt(grau2)
  remod =  abs(resid)

  if (grmod < zesho) return

  ficve = remod / grmod

  if (ficve < zesho) return

  shfau = shfac
  if (difeq > zesh2) then
     shpec = ficve * hconv / 2.0_rp / difeq
     shfau = shfac - 1.0 / shpec
     if (shfau < 0.0_rp) shfau = 0.0_rp
  end if
  
!  shfau= shfau * qufac
  epsdc= 0.5_rp * shfau * hconv * ficve 

  fiso    = epsdc
  faniso  = 0.0_rp

  if (ifcdr == 1) then
     epssu = taupa * xvel2
     epssl = epsdc - epssu
     if (epssl < 0.0_rp) epssl= 0.0_rp
     faniso= (epssl - epsdc)/xvel2
  end if
  
  uvel= xvelo(1)
  vvel= xvelo(2)
  wvel= xvelo(ndime)

  shmet(1,1) =  faniso*uvel*uvel + fiso
  shmet(2,2) =  faniso*vvel*vvel + fiso
  shmet(3,3) =  faniso*wvel*wvel + fiso

  shmet(1,2)= faniso*uvel*vvel
  shmet(1,3)= faniso*uvel*wvel
  shmet(2,3)= faniso*vvel*wvel
  shmet(2,1)= shmet(1,2)
  shmet(3,1)= shmet(1,3)
  shmet(3,2)= shmet(2,3)

  shote = 0.0_rp
  do jdime= 1,ndime
     do idime= 1,ndime
        shote= shote + (cartd(idime)*shmet(jdime,idime))*graun(jdime)
     end do
  end do

  shote = -shote

end subroutine nsa_shocap

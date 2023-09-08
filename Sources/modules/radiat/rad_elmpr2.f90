subroutine rad_elmpr2(&
     pnode,pgaus,pmate,gpsgs,&
     gpsha,elrad,elcod,gppro,&
     gpcod,lnods)
  !-----------------------------------------------------------------------
! Variables needed
! gpdif,gprea,gpden,gpvel,gppro,gpvol,gpgrd,gprhs,gpcar,gphes,elrad,elcod,gpsta,gpdiv,gplap,gpsgs,elmat,elrhs
  !****f* Radiat/rad_elmpr2
  ! NAME
  !   rad_elmpr2
  ! DESCRIPTION
  !    Compute coefficients that go into the ADR equation for P1 model
  ! OUTPUT 
  !    GPGRD(NDIME) ... grad(Gamma) coefficient
  !    GPRAD .......... Radiation of previous iterations and time steps
  !    GPVEL .......... Advection a = 0
  !    GPADV .......... Advection term a.grad(Ni) = 0
  ! USES
  ! USED BY
  !    rad_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  mnode,mgaus,ndime,ntens,kfl_naxis
  use def_radiat, only       :  &
       &                        kfl_exacs_rad
  implicit none 
  integer(ip), intent(in)    :: pnode,pgaus,pmate,lnods(pnode)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpsgs(pgaus,*)
  real(rp),    intent(in)    :: elrad(pnode,*)
  real(rp),    intent(in)    :: elcod(ndime,pnode)
  real(rp),    intent(out)   :: gppro(pgaus)
  real(rp),    intent(out)   :: gpcod(ndime,pgaus)
  integer(ip)                :: idime,inode,igaus,itime,ipoin
  real(rp)                   :: dummr,gpdkt

  !
  ! Coordinates
  !
  if(kfl_exacs_rad/=0.or.kfl_naxis==1) then
     do igaus=1,pgaus
        do idime=1,ndime
           gpcod(idime,igaus)=0.0_rp
        end do
        do inode=1,pnode
           do idime=1,ndime
              gpcod(idime,igaus)=gpcod(idime,igaus)&
                   +gpsha(inode,igaus)*elcod(idime,inode)
           end do
        end do
     end do
  end if
  !
  ! Orthogonal SGS
  !
  do igaus = 1,pgaus
     gppro(igaus) = 0.0_rp
  end do
!!$  if( kfl_ortho_rad >= 1 ) then
!!$     do igaus = 1,pgaus
!!$        do inode = 1,pnode
!!$           ipoin = lnods(inode)
!!$           gppro(igaus) = gppro(igaus) + gpsha(inode,igaus) * rapro_rad(ipoin)
!!$        end do
!!$     end do
!!$  end if

end subroutine rad_elmpr2

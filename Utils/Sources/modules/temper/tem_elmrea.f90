subroutine tem_elmrea(&
     itask,pnode,pgaus,igaui,igauf,elvel,gpden,gpcar,gprea)   
  !-----------------------------------------------------------------------
  !****f* Temper/tem_elmpro
  ! NAME 
  !    tem_elmpro
  ! DESCRIPTION
  !    Compute reaction term
  ! USES
  ! USED BY
  !    tem_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,mnode
  use def_temper, only     :  react_tem,kfl_regim_tem
  use def_kermod, only     :  gasco
  use mod_ADR,    only     :  mreac_adr

  implicit none
  integer(ip), intent(in)  :: itask,pnode,pgaus,igaui,igauf
  real(rp),    intent(in)  :: elvel(ndime,pnode)
  real(rp),    intent(in)  :: gpden(pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(out) :: gprea(pgaus,mreac_adr)
  integer(ip)              :: igaus,idime,inode,ipoin 
  real(rp)                 :: gpdiv

  if( itask == 1 ) then
     !
     ! Reaction: GPREA=r
     !
     do igaus = igaui,igauf
        gprea(igaus,1) = react_tem
     end do
     !
     ! Compressibility term: GPREA=r + rho*R*div(u)
     !
     if( kfl_regim_tem == 1 .or. kfl_regim_tem == 2 ) then
        do igaus = igaui,igauf
           gpdiv = 0.0_rp
           do inode = 1,pnode
              do idime = 1,ndime
                 gpdiv = gpdiv &
                      + elvel(idime,inode) * gpcar(idime,inode,igaus)
              end do
           end do
           gprea(igaus,1) = gprea(igaus,1) + gpden(igaus) * gasco * gpdiv
        end do
     end if

  end if

end subroutine tem_elmrea

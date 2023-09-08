subroutine chm_vdiffu(&
     iclas,pnode,pelty,pgaus,gpcar,elvel,elcod,gpden,gpdif)
  !-----------------------------------------------------------------------
  !****f* chemic/chm_vdiffu
  ! NAME 
  !    chm_vdiffu
  ! DESCRIPTION
  !    Compute vertical diffusion
  ! USES
  ! USED BY
  !    chm_elmpro
  !***
  ! NOTE: diffu_chm contains the diffusion coefficient (m^2/2). In the 
  !       conservative form we solve using densi*diffu (kg/m.s). So, it
  !       in necessary to multiply be density... 
  !
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,mgaus,mnode,hnatu,elmar
  use def_chemic, only      :  lawdi_chm,diffu_chm
  implicit none
  integer(ip),  intent(in)  :: iclas,pnode,pelty,pgaus
  real(rp),     intent(in)  :: elvel(ndime,pnode)
  real(rp),     intent(in)  :: elcod(ndime,pnode)
  real(rp),     intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),     intent(in)  :: gpden(pgaus)
  real(rp),     intent(out) :: gpdif(ndime,pgaus)
  integer(ip)               :: igaus
  
  do igaus = 1,pgaus
     gpdif(ndime,igaus) = 0.0_rp      ! Only vertical component
  end do

  if( lawdi_chm(2,iclas) == 1 ) then
     !
     ! Constant
     !
     do igaus = 1,pgaus
        gpdif(ndime,igaus) = gpden(igaus)*diffu_chm(2,iclas)  ! multiply by densi
     end do

  else if( lawdi_chm(2,iclas) == 2 ) then
     !
     ! Similarity
     !
     call runend('SIMILARITY DIFFUSION TO BE CODED')

  end if
 
end subroutine chm_vdiffu
 

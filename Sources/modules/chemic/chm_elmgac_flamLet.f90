subroutine chm_elmgac_flamLet(&
     ielem,pnode,lnods,elcod,elcon,elvel,elmas,&
     elkey,eleps,elrrt)
  !------------------------------------------------------------------------
  !****f* Chemic/chm_elmgac_cfi
  ! NAME 
  !    chm_elmgac_cfi
  ! DESCRIPTION
  !    Gather operations for the combustion models
  ! USES
  ! USED BY
  !    chm_elmcfi
  !***
  !------------------------------------------------------------------------ 
  use def_kintyp, only     :  ip,rp 
  use def_domain, only     :  ndime,mnode,coord
  use def_master, only     :  conce,untur,&
                              lescl,massk,advec
  use def_chemic, only     :  kfl_advec_chm,&
       &                      nclas_chm,&
       &                      kfl_cotur_chm,&
       &                      ADR_chm
  use def_kermod, only     :  gasco
  use mod_ADR,    only     :  BDF
  use mod_ker_proper

  implicit none 
  integer(ip), intent(in)  :: pnode,ielem
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(out) :: elcod(ndime,pnode)
  real(rp),    intent(out) :: elcon(pnode,nclas_chm,*)
  real(rp),    intent(out) :: elrrt(pnode)                          ! Transport of reaction rate by turbulent fluctuation {w_c * c}
  real(rp),    intent(out) :: elvel(ndime,pnode)
  real(rp),    intent(out) :: elmas(pnode,nclas_chm)                ! Mass source terms
  real(rp),    intent(out) :: elkey(pnode),eleps(pnode)
  integer(ip)              :: inode,ipoin,idime,itime,iclas,dummi

  !
  ! Concentration and coordinates
  !
  do inode=1,pnode
     ipoin=lnods(inode)
     do iclas=1,nclas_chm
        elcon(inode,iclas,1) = conce(ipoin,iclas,1)
     end do
     do idime=1,ndime
        elcod(idime,inode)   = coord(idime,ipoin)
     end do
  end do

  !
  ! Time integration
  !
  if( ADR_chm(1) % kfl_time_integration == 1 ) then
     do iclas = 1,nclas_chm
        do inode = 1,pnode
           ipoin = lnods(inode)
           elcon(inode,iclas,2) = conce(ipoin,iclas,3)
        end do
     end do

     if( ADR_chm(1) % kfl_time_scheme == BDF ) then
        do iclas = 1,nclas_chm
           do itime = 3,ADR_chm(iclas) % kfl_time_order + 1
              do inode = 1,pnode
                 ipoin = lnods(inode)
                 elcon(inode,iclas,itime) = conce(ipoin,iclas,itime+1)
              end do
           end do
        end do
     end if

  end if

  !
  ! Advection
  !
  if( kfl_advec_chm /= 0 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        do idime = 1,ndime
           elvel(idime,inode) = advec(idime,ipoin,1)
        end do
     end do
  else
     elvel = 0.0_rp
  end if

  !
  ! TKE and EPS for RANS model
  !
  elkey = 0.0_rp
  eleps = 0.0_rp
  elrrt = 0.0_rp 

  do inode=1,pnode
     ipoin=lnods(inode)
     elrrt(inode) = lescl(ipoin)
  end do

  if (kfl_cotur_chm == 1_ip ) then       !! K-EPSILON MODEL
     do inode=1,pnode
        ipoin=lnods(inode)
        elkey(inode) = untur(1,ipoin,1)
        eleps(inode) = untur(2,ipoin,1)
     end do
  elseif (kfl_cotur_chm == 2_ip ) then   !! K-OMEGA MODEL
     do inode=1,pnode
        ipoin=lnods(inode)
        elkey(inode) = untur(1,ipoin,1)
        eleps(inode) = 0.09_rp * untur(1,ipoin,1) * untur(2,ipoin,1)
     end do
  endif       

  !
  ! Mass source terms coefficients
  !
  do iclas = 1, nclas_chm
     do inode = 1,pnode
        ipoin = lnods(inode)
        elmas(inode,iclas) = massk(ipoin,iclas)
     enddo
  end do
    
end subroutine chm_elmgac_flamLet

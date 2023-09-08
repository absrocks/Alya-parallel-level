subroutine chm_elmgtt(&
     pnode,lnods,elcod,elcon,elpro,elvel,elvte,eltmr,elden)
  !------------------------------------------------------------------------
  !****f* partis/chm_elmgtt
  ! NAME 
  !    chm_elmgat
  ! DESCRIPTION
  !    This routine performs the gather operations
  ! USES
  ! USED BY
  !    chm_elmope
  !***
  !------------------------------------------------------------------------ 
  use def_kintyp, only     :  ip,rp 
  use def_domain, only     :  ndime,mnode,npoin,coord
  use def_master, only     :  conce
  use def_chemic, only     :  kfl_timei_chm,kfl_advec_chm,&
       &                      kfl_tiacc_chm,kfl_tisch_chm,&
       &                      nclas_chm,nspec_chm,iclas_chm,&
       &                      proje_chm,kfl_stabi_chm,kfl_advec_chm,&
       &                      veloc_chm,vterm_chm,tmrat_chm,densi_chm,&
       &                      lawvt_chm,kfl_sourc_chm,lawde_chm
  implicit none 
  integer(ip), intent(in)  :: pnode
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(out) :: elcod(ndime,pnode)
  real(rp),    intent(out) :: elcon(pnode,nspec_chm,*)
  real(rp),    intent(out) :: elpro(pnode)
  real(rp),    intent(out) :: elvel(ndime,pnode)
  real(rp),    intent(out) :: elvte(pnode)
  real(rp),    intent(out) :: eltmr(pnode)
  real(rp),    intent(out) :: elden(pnode)
  integer(ip)              :: inode,ipoin,idime,itime,iclas
  !
  ! Current concentration and coordinates
  !
  do inode=1,pnode
     ipoin=lnods(inode)
     do iclas=1,nspec_chm
        elcon(inode,iclas,1)=conce(ipoin,iclas,1)
     end do
     do idime=1,ndime
        elcod(idime,inode)=coord(idime,ipoin)
     end do  
  end do
  !
  ! Time integration
  !
  if( kfl_timei_chm == 1 ) then
     do iclas = 1,nclas_chm
        do inode = 1,pnode
           ipoin = lnods(inode)
           elcon(inode,iclas,2) = conce(ipoin,iclas,3)
        end do
     end do
     if( kfl_tisch_chm == 2 ) then
        do iclas = 1,nclas_chm
           do itime = 3,1+kfl_tiacc_chm
              do inode = 1,pnode
                 ipoin = lnods(inode)
                 elcon(inode,iclas,itime) = conce(ipoin,iclas,itime+1)
              end do
           end do
        end do
     end if
  end if 
  !
  ! Projection
  !
  if( kfl_stabi_chm >= 1 ) then
    do inode = 1,pnode
       ipoin = lnods(inode)
       elpro(inode) = proje_chm(ipoin,iclas_chm)
    end do
  end if
  !
  ! Advection
  !
  if( kfl_advec_chm < 0 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        do idime = 1,ndime
           elvel(idime,inode) = veloc_chm(idime,ipoin)
        end do
     end do
  end if
  !
  ! Terminal velocity
  !
  if( lawvt_chm >= 1 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        elvte(inode) = vterm_chm(ipoin,iclas_chm)
     end do
  end if
  !
  ! Source term
  !
  if( kfl_sourc_chm < 0 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        eltmr(inode) = tmrat_chm(iclas_chm,ipoin)
     end do
  end if
  !
  ! Density
  !
  if( lawde_chm /= 0 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        elden(inode) = densi_chm(ipoin)
     end do     
  end if

end subroutine chm_elmgtt

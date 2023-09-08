!------------------------------------------------------------------------
!> @addtogroup Partis
!! @{
!> @name    Partis memory
!! @file    pts_outvar.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!! @brief   Postprocess particle results on nodes
!! @details Postprocess particle results:
!!          - Number of particles per type
!!          - Desposited particles per type
!> @} 
!------------------------------------------------------------------------

subroutine pts_outvar(ivari)
  use def_parame
  use def_master
  use def_domain
  use def_partis
  use mod_memory
  use mod_gradie
  use mod_postpr
  use def_kermod, only :  ndivi
  use mod_projec, only :  projec_elements_to_nodes
  use mod_projec, only :  projec_boundaries_to_nodes
  implicit none
  integer(ip), intent(in) :: ivari
  integer(ip)             :: ipoin,igaus,itype,idime,iboun
  integer(ip)             :: ilagr,ielem,pelty,pnode,pgaus,inode
  real(rp)                :: elcod(ndime,mnode)
  real(rp)                :: gpvol,gpdet
  real(rp)                :: xjacm(9) 
  real(rp),    pointer    :: gelem(:,:)
  real(rp),    pointer    :: depob_tmp(:)

  if( ivari == 0 ) return

  select case ( ivari )  

  case( 1_ip )
     !
     ! Number of particles/node
     !
     if( INOTMASTER ) then
        nullify(gelem)
        call memgen(0_ip,ntyla_pts,npoin)
        call memory_alloca(mem_modul(1:2,modul),'GELEM','pts_outvar',gelem,ntyla_pts,nelem)
        !
        ! Count particles inside each element
        !
        do ilagr = 1,mlagr 
           if( lagrtyp(ilagr) % kfl_exist == -1 ) then 
              itype = lagrtyp(ilagr) % itype
              ielem = lagrtyp(ilagr) % ielem
              if( ielem > 0 .and. ielem <= nelem ) then
                 gelem(itype,ielem) = gelem(itype,ielem) + 1.0_rp
              end if
           end if
        end do
        call projec_elements_to_nodes(gelem,gevec) 
        call memory_deallo(mem_modul(1:2,modul),'GELEM','pts_outvar',gelem)
     else
        !
        ! This is necessary for outvar to get the right dimension
        !
        call memgen(0_ip,ntyla_pts,1_ip)
     end if

  case( 2_ip )
     !
     ! Deposition map
     !
     if( kfl_depos_pts == 0 ) return
     postp(1) % wopos(3,2) = 'NPOIN'
     if( INOTMASTER ) then
        call memgen(0_ip,ntyla_pts,npoin)
        call projec_boundaries_to_nodes(depob_pts,meshe(ndivi),elmar,gevec)
     else
        call memgen(0_ip,ntyla_pts,1_ip)
     end if
     
  case( 3_ip )
     !
     ! Slip wall distance
     !     
     gesca => walld_slip_pts

  case( 4_ip )
     !
     ! Friction
     !     
     gesca => friction_pts

  case( 5_ip )
     !
     ! Residence
     !     
     if( kfl_resid_pts == 0 ) return
     if( INOTMASTER ) then
        call memgen(0_ip,ntyla_pts,npoin)
        call projec_elements_to_nodes(resid_pts,gevec) 
     else
        call memgen(0_ip,ntyla_pts,1_ip)
     end if

  case( 6_ip )
     !
     ! Deposition map on boundaries
     !
     if( INOTMASTER ) then
        call memgen(0_ip,nboun,0_ip)
        do iboun = 1,nboun
           gesca(iboun) = depob_pts(1,iboun)
        end do
     end if
     
  case( 7_ip )
     !
     ! Total deposition map
     !
     if( kfl_depos_pts == 0 ) return
     
     if( INOTMASTER ) then
        nullify(depob_tmp)
        call memory_alloca(mem_modul(1:2,modul),'DEPOB_TMP','pts_outvar',depob_tmp,max(1_ip,nboun))
        do iboun = 1,nboun
           do itype = 1,ntyla_pts
              if( parttyp(itype) % kfl_exist /= 0 ) then
                 depob_tmp(iboun) = depob_tmp(iboun) + depob_pts(itype,iboun)
              end if
           end do
        end do        
        call memgen(0_ip,npoin,0_ip)
        call projec_boundaries_to_nodes(depob_tmp,meshe(ndivi),elmar,gesca)
        call memory_deallo(mem_modul(1:2,modul),'DEPOB_TMP','pts_outvar',depob_tmp)
     end if
     
  case( 8_ip )
     !
     ! Bouncing wall distance
     !     
     gesca => walld_bouncing_pts

  end select

  call outvar(&
       ivari,&
       ittim,cutim,postp(1) % wopos(1,ivari))

  if( ivari == 2 ) postp(1) % wopos(3,2) = 'NBOUN'
  
end subroutine pts_outvar

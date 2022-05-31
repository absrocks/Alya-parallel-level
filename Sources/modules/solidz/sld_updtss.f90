!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_updtss.f90
!> @author  Solidz Team
!> @date    November, 2017
!> @author  Gerard Guillamet
!> @date    August, 2018
!>          - Parallelization using OpenMP
!> @brief   Calculation of the stable time increment
!> @details Calculates the stable time increment (STI) for dynamic
!>          (transient) analysis. The stable time increment or critical
!>          time step is calculated using the speed of sound from the
!>          material and the element characteristic length.
!>
!>          Default elements
!>          dtcri = sf * Le * \sqrt(rho/E)
!>
!>          Cohesive elements
!>          dtcri = \sqrt(rho_tilde/Kcoh)
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_updtss()

  use def_kintyp,                only : ip,rp
  use def_master,                only : dtinv,INOTMASTER,kfl_timco,zeror
  use def_domain,                only : ndime,mnode,nelem
  use def_domain,                only : ltype,lelch,lmate,lnods,lorde,nnode
  use def_domain,                only : coord,hnatu,elmar
  use def_elmtyp,                only : SHELL, BAR3D, ELINT
  use mod_communications,        only : PAR_MIN
  use mod_parall,                only : par_omp_nelem_chunk
  use def_solidz,                only : kfl_timei_sld, SLD_STATIC_PROBLEM
  use def_solidz,                only : kfl_rigid_sld, celen_sld
  use def_solidz,                only : kfl_xfeme_sld, eleng_sld
  use def_solidz,                only : velas_sld
  use def_solidz,                only : dtinv_sld, dtcri_sld, safet_sld
  use def_solidz,                only : kfl_pseud_sld, dttau_sld, safet_pseud_sld
  use mod_sld_interface_element, only : ELINT_stable_time_increment

  implicit none

  integer(ip)                    :: ielem,idime,inode,ipoin
  integer(ip)                    :: pnode,pmate,pelty,porde
  real(rp)                       :: dtcri
  real(rp)                       :: dtmin
  real(rp)                       :: hmini
  real(rp)                       :: elcod(ndime,mnode)
  real(rp)                       :: hleng(3)
  real(rp)                       :: tragl(9)
  real(rp)                       :: dummr

  !
  ! Stable time increment (minimum of element time steps)
  !
  if ( kfl_timei_sld /= SLD_STATIC_PROBLEM .and. kfl_timco /= 0 ) then

     if ( INOTMASTER ) then

        dtmin = huge(1.0_rp)

        !--------------------------------------------------------------------------
        !$OMP PARALLEL  DO                                                        &
        !$OMP SCHEDULE  ( DYNAMIC , par_omp_nelem_chunk )                         &
        !$OMP DEFAULT   (NONE)                                                    &
        !$OMP PRIVATE   (ielem,pelty,pnode,porde,pmate,inode,ipoin,idime,         &
        !$OMP            elcod,tragl,hleng,dummr,hmini,dtcri)                     &
        !$OMP SHARED    (ltype,nnode,lorde,lmate,lnods,coord,elmar,hnatu,nelem,   &
        !$OMP            lelch,kfl_xfeme_sld,kfl_rigid_sld,eleng_sld,velas_sld,   &
        !$OMP            dttau_sld,celen_sld,safet_sld,                           &
        !$OMP            kfl_pseud_sld, safet_pseud_sld,                          &
#ifndef NDIMEPAR
        !$OMP            ndime,                                                   &
#endif
        !$OMP            par_omp_nelem_chunk)                                     &
        !$OMP REDUCTION (MIN:dtmin)
        !--------------------------------------------------------------------------
        !
        ! Loop over elements
        !
        elements: do ielem = 1,nelem

           pelty = ltype(ielem)

           if ( pelty > 0 .and. pelty /= SHELL .and. pelty /= BAR3D .and. lelch(ielem) /= ELINT ) then
              pnode = nnode(pelty)
              porde = lorde(pelty)
              pmate = lmate(ielem)
              !
              ! Default: elcod
              !
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 elcod(1:ndime,inode) = coord(1:ndime,ipoin)
              end do
              !
              ! Speed of sound according to the material
              !
              if ( kfl_rigid_sld == 0_ip ) then
                 call sld_velsnd(pmate)
              else
                 velas_sld(1,pmate) = 1.0_rp
              end if
              !
              ! Element characteristic lenght
              !
              call elmlen(&
                   ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                   hnatu(pelty),hleng)
              !
              ! Minimum length of the element
              !
              hmini = hleng(ndime) / real(porde,rp)
              celen_sld(ielem) = hmini
              !
              ! X-FEM
              !
              if ( kfl_xfeme_sld > 0 ) then
                 dummr = 0.0_rp
                 do idime = 1,ndime
                    dummr = dummr + hleng(idime)
                 end do
                 eleng_sld(ielem) = dummr / real(ndime,rp)
              end if
              !
              ! Critical time increment
              !
              dtcri = safet_sld*(hmini/velas_sld(1,pmate))
              if( kfl_pseud_sld == 1 ) dttau_sld(ielem) = dtcri*safet_pseud_sld/safet_sld

           else if ( lelch(ielem) == ELINT ) then
              !
              ! Interface elements
              !
              call ELINT_stable_time_increment(ielem,dtcri)

           end if
           !
           ! Min. time increment
           !
           dtmin = min( dtmin, dtcri )

        end do elements
        !$OMP END PARALLEL DO
        !--------------------------------------------------------------------------

     end if
     !
     ! Parall: Look for minimum over all subdomains (dtmin)
     !
     call PAR_MIN(dtmin,'IN MY CODE')
     !
     ! Assign 1/dt
     !
     dtcri_sld = dtmin
     if ( dtcri_sld /= 0.0_rp ) dtinv_sld = 1.0_rp/dtcri_sld
     if ( kfl_timco == 1 )      dtinv     = max(dtinv,dtinv_sld)

  end if

end subroutine sld_updtss

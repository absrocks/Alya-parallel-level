!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_outwig.f90
!> @author  houzeaux
!> @date    2020-01-29
!> @brief   Witness output
!> @details Output of witness geometries
!> @} 
!-----------------------------------------------------------------------

subroutine pts_outwig()

  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_partis
  use mod_ker_proper,     only : ker_proper
  use mod_communications, only : PAR_SUM
  use mod_witness,        only : witness_in_geometry
  use mod_witness,        only : WITNESS_RING 
  use mod_witness,        only : witness_ring_transform
  use mod_witness,        only : witness_volume_geometry
  use mod_physics,        only : physics_set_liquid_temperature
  use mod_pts_particle,   only : pts_particle_diameter
  implicit none
  integer(ip) :: iwitg,ilagr,ilagr_local,idime
  real(rp)    :: num_part,dummr,denpa,diame,dummr_vel(ndime),VolCV 

  do iwitg = 1,nwitg
     !
     ! Almost all variables need the number of paritcles in the control volume:
     !
     num_part = 0.0_rp
     do ilagr_local = 1,nlagr_local_pts
        ilagr = permu_nlagr_pts(ilagr_local)
        if( witness_in_geometry( lagrtyp(ilagr) % coord,gewit(iwitg) % kfl_geometry,gewit(iwitg) % param) ) then
           num_part = num_part + parttyp(lagrtyp(ilagr) % itype) % n_drop
        end if
     end do
     call PAR_SUM(num_part)

     if( postp(1) % npp_witng(1) > 0 ) then
        ! 
        !                      \int_t n(t) dt
        ! Number of particles: ------------
        !                        \int_t dt
        !
        if( IMASTER ) witng(1,iwitg) = witng(1,iwitg) + num_part
     end if

     if( postp(1) % npp_witng(2) > 0 ) then
        !
        !                   \sum_n d
        ! Average diameter:----------
        !                   \sum_n 1
        !
        do ilagr_local = 1,nlagr_local_pts
           ilagr = permu_nlagr_pts(ilagr_local)
           if( witness_in_geometry( lagrtyp(ilagr) % coord,gewit(iwitg) % kfl_geometry,gewit(iwitg) % param) ) then
              if( parttyp(lagrtyp(ilagr) % itype) % kfl_therm == 0 ) then
                 denpa     = parttyp(lagrtyp(ilagr) % itype) % denpa
              else
                 call physics_set_liquid_temperature( parttyp(lagrtyp(ilagr) % itype) % liq , lagrtyp(ilagr) % tempe_k)
                 denpa     = parttyp(lagrtyp(ilagr) % itype) % liq % rho 
              endif
              witng(2,iwitg) = witng(2,iwitg) + parttyp(lagrtyp(ilagr) % itype) % n_drop * pts_particle_diameter(lagrtyp(ilagr) % itype,ilagr,denpa)
           end if
        end do
        dummr = witng(2,iwitg)
        call PAR_SUM(dummr)
        if( IMASTER ) witng(2,iwitg) = dummr 

        postp(1) % witng_deldenom(2,iwitg) = num_part
     endif

     if( postp(1) % npp_witng(3) > 0 ) then
        !
        !                       \sum_n d^3
        ! Sauter mean diameter:------------
        !                       \sum_n d^2
        !
        postp(1) % witng_deldenom(3,iwitg) = 0.0_rp
        do ilagr_local = 1,nlagr_local_pts
           ilagr = permu_nlagr_pts(ilagr_local)
           if( witness_in_geometry( lagrtyp(ilagr) % coord,gewit(iwitg) % kfl_geometry,gewit(iwitg) % param) ) then
              if( parttyp(lagrtyp(ilagr) % itype) % kfl_therm == 0 ) then
                 denpa     = parttyp(lagrtyp(ilagr) % itype) % denpa
              else
                 call physics_set_liquid_temperature( parttyp(lagrtyp(ilagr) % itype) % liq , lagrtyp(ilagr) % tempe_k)
                 denpa     = parttyp(lagrtyp(ilagr) % itype) % liq % rho 
              endif
              diame                           = pts_particle_diameter(lagrtyp(ilagr) % itype,ilagr,denpa)
              witng(3,iwitg)                  = witng(3,iwitg) + parttyp(lagrtyp(ilagr) % itype) % n_drop * diame**3
              postp(1) % witng_deldenom(3,iwitg) = postp(1) % witng_deldenom(3,iwitg) + parttyp(lagrtyp(ilagr) % itype) % n_drop * diame**2
           end if
        end do
        dummr = witng(3,iwitg)
        call PAR_SUM(dummr)
        if( IMASTER ) witng(3,iwitg) = dummr 
        dummr = postp(1) % witng_deldenom(3,iwitg)
        call PAR_SUM(dummr)
        if( IMASTER ) postp(1) % witng_deldenom(3,iwitg) = dummr 

     endif

     if( postp(1) % npp_witng(4) > 0 ) then
        !
        !                      \sum_n T
        ! Average temperature:----------
        !                      \sum_n 1
        !
        do ilagr_local = 1,nlagr_local_pts
           ilagr = permu_nlagr_pts(ilagr_local)
           if( witness_in_geometry( lagrtyp(ilagr) % coord,gewit(iwitg) % kfl_geometry,gewit(iwitg) % param) ) then
              witng(4,iwitg) = witng(4,iwitg) + parttyp(lagrtyp(ilagr) % itype) % n_drop * lagrtyp(ilagr) % tempe_k
           end if
        end do
        dummr = witng(4,iwitg)
        call PAR_SUM(dummr)
        if( IMASTER ) witng(4,iwitg) = dummr 

        postp(1) % witng_deldenom(4,iwitg) = num_part

     endif

     if( postp(1) % npp_witng(5) > 0 ) then
        !
        !                   \sum_n U
        ! Average velocity:----------
        !                   \sum_n 1
        !
        do ilagr_local = 1,nlagr_local_pts
           ilagr = permu_nlagr_pts(ilagr_local)
           if( witness_in_geometry( lagrtyp(ilagr) % coord,gewit(iwitg) % kfl_geometry,gewit(iwitg) % param) ) then
              !
              ! Transform to polar coordinates for ring
              !
              dummr_vel = lagrtyp(ilagr) % veloc
              if (gewit(iwitg) % kfl_geometry == WITNESS_RING) call witness_ring_transform(lagrtyp(ilagr) % coord,gewit(iwitg) % param,vec=dummr_vel)
              do idime = 1,ndime
                 witng(4+idime,iwitg) = witng(4+idime,iwitg) + parttyp(lagrtyp(ilagr) % itype) % n_drop * dummr_vel(idime)
              enddo
           end if
        end do

        do idime = 1,ndime
           dummr_vel(idime) = witng(4+idime,iwitg)
           postp(1) % witng_deldenom(4+idime,iwitg) = num_part
        enddo
        call PAR_SUM(ndime,dummr_vel)
        if( IMASTER ) then
           do idime = 1,ndime
              witng(4+idime,iwitg) = dummr_vel(idime) 
           enddo
        endif
     endif

     if( postp(1) % npp_witng(8) > 0 ) then
        !
        !                          \sum_n d^2
        ! Average diameter square:------------
        !                           \sum_n 1
        !
        do ilagr_local = 1,nlagr_local_pts
           ilagr = permu_nlagr_pts(ilagr_local)
           if( witness_in_geometry( lagrtyp(ilagr) % coord,gewit(iwitg) % kfl_geometry,gewit(iwitg) % param) ) then
              if( parttyp(lagrtyp(ilagr) % itype) % kfl_therm == 0 ) then
                 denpa     = parttyp(lagrtyp(ilagr) % itype) % denpa
              else
                 call physics_set_liquid_temperature( parttyp(lagrtyp(ilagr) % itype) % liq , lagrtyp(ilagr) % tempe_k)
                 denpa     = parttyp(lagrtyp(ilagr) % itype) % liq % rho 
              endif
              witng(8,iwitg) = witng(8,iwitg) + parttyp(lagrtyp(ilagr) % itype) % n_drop * pts_particle_diameter(lagrtyp(ilagr) % itype,ilagr,denpa)**2
           end if
        end do
        dummr = witng(8,iwitg)
        call PAR_SUM(dummr)
        if( IMASTER ) witng(8,iwitg) = dummr 

        postp(1) % witng_deldenom(8,iwitg) = num_part
     endif

     if( postp(1) % npp_witng(9) > 0 ) then
        !
        !                             \sum_n T^2
        ! Average temperature square:------------
        !                              \sum_n 1
        !
        do ilagr_local = 1,nlagr_local_pts
           ilagr = permu_nlagr_pts(ilagr_local)
           if( witness_in_geometry( lagrtyp(ilagr) % coord,gewit(iwitg) % kfl_geometry,gewit(iwitg) % param) ) then
              witng(9,iwitg) = witng(9,iwitg) + parttyp(lagrtyp(ilagr) % itype) % n_drop * lagrtyp(ilagr) % tempe_k**2
           end if
        end do
        dummr = witng(9,iwitg)
        call PAR_SUM(dummr)
        if( IMASTER ) witng(9,iwitg) = dummr 

        postp(1) % witng_deldenom(9,iwitg) = num_part

     endif

     if( postp(1) % npp_witng(10) > 0 ) then
        !
        !                          \sum_n U^2
        ! Average velocity square:------------
        !                           \sum_n 1
        !
        do ilagr_local = 1,nlagr_local_pts
           ilagr = permu_nlagr_pts(ilagr_local)
           if( witness_in_geometry( lagrtyp(ilagr) % coord,gewit(iwitg) % kfl_geometry,gewit(iwitg) % param) ) then
              dummr_vel = lagrtyp(ilagr) % veloc
              if (gewit(iwitg) % kfl_geometry == WITNESS_RING) call witness_ring_transform(lagrtyp(ilagr) % coord,gewit(iwitg) % param,vec=dummr_vel)
              do idime = 1,ndime
                 witng(9+idime,iwitg) = witng(9+idime,iwitg) + parttyp(lagrtyp(ilagr) % itype) % n_drop * dummr_vel(idime)**2
              enddo
           end if
        end do

        do idime = 1,ndime
           dummr_vel(idime) = witng(9+idime,iwitg)
           postp(1) % witng_deldenom(9+idime,iwitg) = num_part
        enddo
        call PAR_SUM(ndime,dummr_vel)
        if( IMASTER ) then
           do idime = 1,ndime
              witng(9+idime,iwitg) = dummr_vel(idime) 
           enddo
        endif
     endif

     if( postp(1) % npp_witng(13) > 0 ) then
        !
        !               \sum_n m
        ! Average mass:-----------
        !               \sum_n 1
        !
        do ilagr_local = 1,nlagr_local_pts
           ilagr = permu_nlagr_pts(ilagr_local)
           if( witness_in_geometry( lagrtyp(ilagr) % coord,gewit(iwitg) % kfl_geometry,gewit(iwitg) % param) ) then
              witng(13,iwitg) = witng(13,iwitg) + parttyp(lagrtyp(ilagr) % itype) % n_drop * lagrtyp(ilagr) % mass_k
           end if
        end do
        dummr = witng(13,iwitg)
        call PAR_SUM(dummr)
        if( IMASTER ) witng(13,iwitg) = dummr 

        postp(1) % witng_deldenom(13,iwitg) = num_part

     endif

     if( postp(1) % npp_witng(14) > 0 ) then
        !
        !                           \sum_n m U             \sum_n m U               
        ! Average massu flux:-------------------------- = ------------ 
        !                     Area*thickness  \sum_n 1     V \sum_n 1
        !
        VolCV = witness_volume_geometry( gewit(iwitg) % kfl_geometry,gewit(iwitg) % param)
        do ilagr_local = 1,nlagr_local_pts
           ilagr = permu_nlagr_pts(ilagr_local)
           if( witness_in_geometry( lagrtyp(ilagr) % coord,gewit(iwitg) % kfl_geometry,gewit(iwitg) % param) ) then
              !
              ! Transform to polar coordinates for ring
              !
              dummr_vel = lagrtyp(ilagr) % veloc
              if (gewit(iwitg) % kfl_geometry == WITNESS_RING) call witness_ring_transform(lagrtyp(ilagr) % coord,gewit(iwitg) % param,vec=dummr_vel)
              do idime = 1,ndime
                 witng(13+idime,iwitg) = witng(13+idime,iwitg) + parttyp(lagrtyp(ilagr) % itype) % n_drop * dummr_vel(idime) * lagrtyp(ilagr) % mass_k / VolCV
              enddo
           end if
        end do

        do idime = 1,ndime
           dummr_vel(idime) = witng(13+idime,iwitg)
           postp(1) % witng_deldenom(13+idime,iwitg) = num_part
        enddo
        call PAR_SUM(ndime,dummr_vel)
        if( IMASTER ) then
           do idime = 1,ndime
              witng(13+idime,iwitg) = dummr_vel(idime) 
           enddo
        endif

     endif

  end do

end subroutine pts_outwig


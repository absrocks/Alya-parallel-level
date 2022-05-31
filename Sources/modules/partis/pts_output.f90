!------------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @name    Starts an iteration
!> @file    pts_output.f90
!> @author  Guillaume Houzeaux
!> @date    28/01/2013
!> @brief   Output of particles
!> @details Postprocess on mesh (particle density) and output
!> @}
!------------------------------------------------------------------------

subroutine pts_output()
  use def_parame
  use def_master
  use def_kermod
  use def_partis
  use mod_memory
  use mod_ker_timeline,        only : ker_timeline
  use mod_communications,      only : PAR_GATHER,PAR_GATHERV
  use mod_pts_parallelization, only : pts_parallelization_pack
  use mod_pts_parallelization, only : pts_parallelization_unpack
  use mod_messages,            only : livinf
  use mod_outfor,              only : outfor
  use mod_iofile,              only : iofile_flush_unit
  use mod_output_postprocess,  only : output_postprocess_variables
  
  implicit none

  external                 :: pts_outvar
  integer(ip)              :: ivari,ivarp,ilagr,ipars,ipart,ivard
  integer(ip)              :: kfl_ifpos,kfl_ifdep,ilimi,nvart_pts
  integer(4)               :: nvarp_pts4,nvard_pts4,nvart_pts4
  integer(ip), save        :: ittim_last=-1
  real(rp)                 :: cutim_pts
  real(rp),    allocatable :: depos_surface(:)
  integer(4)               :: num_particles
  integer(4),  pointer     :: num_particles_gat(:)
  real(rp),    pointer     :: particles_gat(:)
  real(rp),    pointer     :: particles_send(:)
  type(latyp)              :: particle
  logical(lg), pointer     :: variables_pts(:)
  integer(ip), pointer     :: perm_postprocess(:)
  integer(ip), pointer     :: perm_deposition(:)

  nullify(num_particles_gat) 
  nullify(particles_gat)
  nullify(particles_send)
  nullify(variables_pts)
  nullify(perm_postprocess)
  nullify(perm_deposition)
  
  !----------------------------------------------------------------------
  !
  ! Mesh dependent postprocess
  !
  !----------------------------------------------------------------------

  call output_postprocess_variables(pts_outvar)

  !----------------------------------------------------------------------
  !
  ! Witness
  !
  !----------------------------------------------------------------------

  if( ittyp == ITASK_ENDTIM ) then
     call pts_outwig()
  end if

  !----------------------------------------------------------------------
  !
  ! Particles output
  !
  !----------------------------------------------------------------------
  !
  ! Check what kind of postprocess
  ! KFL_ISPOS = 1: All particles info is required
  ! KFL_IFDEP = 1: Deposited particles info is required
  ! Initial solution if this is a restart is not written
  !
  kfl_ifpos = 0
  kfl_ifdep = 0
  if( mod(ittim, kfl_oufre_pts) == 0 )                   kfl_ifpos = 1
  if( kfl_injec == 1 )                                   kfl_ifpos = 1
  if( nlagr_going_out_pts > 0 .and. kfl_oudep_pts /= 0 ) kfl_ifdep = 1

  if( ittyp == ITASK_INITIA .and. kfl_rstar == 2 ) then
     kfl_ifpos = 0
     kfl_ifdep = 0
  end if
  !
  ! What type of particle to output
  !
  ilimi = -1 
  if(      kfl_injec == 1 ) then
     ilimi = -5                                                ! Only just deposited particles information is gathered
  else if( kfl_injec == 0 ) then
     if( kfl_ifpos == 0 .and. kfl_ifdep == 1 ) then
        ilimi = -2                                             ! Only deposited and vanishing particles information is gathered
     else
        ilimi = -1                                             ! Existing, deposited and vanishing particles information is gathered
     end if
  end if
  if( ittim == ittim_last ) kfl_ifpos = 0                      ! Last time was already postprocessed

  !----------------------------------------------------------------------
  !
  ! MPI_GATHER particles info. Information for postprocess and deposition
  ! are merged in order to avoid extra communication
  !
  !----------------------------------------------------------------------

  if( kfl_ifpos == 1 .or. kfl_ifdep == 1 ) then

     call ker_timeline('INI_OUTPUT')
     !
     ! Save last time we got here
     !
     if( kfl_injec == 1 ) then
        ittim_last = -1
     else
        ittim_last = ittim
     end if
     !
     ! Current time
     !
     if( kfl_injec == 0 ) then
        cutim_pts = cutim
     else
        cutim_pts = cutim-dtime
     end if
     !
     ! Gather
     !
     if( IPARALL ) then
        !
        ! Total number of variables to gather
        !
        call memory_alloca(mem_modul(1:2,modul),'PERM_POSTPROCESS','pts_inivar',perm_postprocess,nvarp_pts)
        call memory_alloca(mem_modul(1:2,modul),'PERM_DEPOSITION', 'pts_inivar',perm_deposition, nvard_pts)
        call memory_alloca(mem_modul(1:2,modul),'VARIABLES_PTS'   ,'pts_inivar',variables_pts,   mvarp_pts)

        nvart_pts = 0
        ivarp     = 0
        ivard     = 0
        do ivari = 1,mvarp_pts
           if(  ( postprocess_var_pts(ivari) .and. kfl_ifpos == 1 ) .or. &
                ( deposition_var_pts(ivari)  .and. kfl_ifdep == 1 ) ) then
              nvart_pts = nvart_pts + 1
           end if
           if( postprocess_var_pts(ivari) .and. kfl_ifpos == 1 ) then
              ivarp                   = ivarp + 1
              perm_postprocess(ivarp) = nvart_pts
              variables_pts(ivari)    = .true.
           end if
           if( deposition_var_pts(ivari)  .and. kfl_ifdep == 1 ) then
              ivard                   = ivard + 1
              perm_deposition(ivard)  = nvart_pts
              variables_pts(ivari)    = .true.
           end if
        end do
        nvarp_pts4 = int(nvarp_pts,4)
        nvard_pts4 = int(nvard_pts,4)
        nvart_pts4 = int(nvart_pts,4)

        num_particles = 0_4
        if( IMASTER ) then
           call memory_alloca(mem_modul(1:2,modul),'NUM_PARTICLES_GAT','pts_output',num_particles_gat,int(npart+1_ip,4_ip),LBOUN=0_4)
        else if( ISLAVE ) then
           num_particles = count( lagrtyp(1:mlagr) % kfl_exist <= ilimi ) * nvart_pts4
           call memory_alloca(mem_modul(1:2,modul),'PARTICLES_SEND','pts_output',particles_send,int(num_particles,ip))
           ipars = 0
           do ilagr = 1,mlagr
              if( lagrtyp(ilagr) % kfl_exist <= ilimi ) then
                 call pts_parallelization_pack(lagrtyp(ilagr),ipars,particles_send,variables_pts)
              end if
           end do
        end if
        call PAR_GATHER(num_particles,num_particles_gat)
        if( IMASTER ) then
           num_particles = sum(num_particles_gat)
           call memory_alloca(mem_modul(1:2,modul),'PARTICLES_GAT','pts_output',particles_gat,num_particles)
        end if
        call PAR_GATHERV(particles_send,particles_gat,num_particles_gat)

     end if

     !----------------------------------------------------------------------
     !
     ! Postprocess
     !
     !----------------------------------------------------------------------
     
     if( kfl_ifpos == 1 )  then

        if( IMASTER ) then

           ipars = 0
           do ipart = 1,npart
              do ilagr = 1,int(num_particles_gat(ipart),ip) / nvart_pts
                 call pts_parallelization_unpack(particle,ipars,particles_gat,variables_pts,'EXIST')
                 if( particle % kfl_exist <= ilimi ) then
                    if (kfl_ptsres_binary>0) then !save binary
                       write(lun_resul_pts) &
                            particles_gat(perm_postprocess(1)+ipars),&
                            int(particles_gat(perm_postprocess(2)+ipars),ip),&
                            int(particles_gat(perm_postprocess(3)+ipars),ip),&                        
                            int(particles_gat(perm_postprocess(4)+ipars),ip),&                        
                            (particles_gat(perm_postprocess(ivarp)+ipars),ivarp=5,nvarp_pts)
                    else !save text
                       write(lun_resul_pts,100) &
                            particles_gat(perm_postprocess(1)+ipars),&
                            int(particles_gat(perm_postprocess(2)+ipars),ip),&
                            int(particles_gat(perm_postprocess(3)+ipars),ip),&                        
                            int(particles_gat(perm_postprocess(4)+ipars),ip),&                        
                            (particles_gat(perm_postprocess(ivarp)+ipars),ivarp=5,nvarp_pts)
                    end if !(kfl_ptsres_binary>0) 
                 end if !( particle % kfl_exist <= ilimi )
                 ipars = ipars + nvart_pts
              end do
           end do
           call iofile_flush_unit(lun_resul_pts)

        else if( ISEQUEN ) then
           call memory_alloca(mem_modul(1:2,modul),'PARTICLES_SEND','pts_output',particles_send,nvarp_pts)
           do ilagr = 1,mlagr
              if( lagrtyp(ilagr) % kfl_exist <= ilimi ) then
                 ipars = 0
                 call pts_parallelization_pack(lagrtyp(ilagr),ipars,particles_send,postprocess_var_pts)

                 if (kfl_ptsres_binary>0) then !save binary
                    write(lun_resul_pts) particles_send(1), & 
                         int(particles_send(2),ip),&
                         int(particles_send(3),ip),&
                         int(particles_send(4),ip),&
                         (particles_send(ivarp),ivarp=5,nvarp_pts)
                 else ! (kfl_ptsres_binary>0) -- text
                    write(lun_resul_pts,100) &
                         particles_send(1),&
                         int(particles_send(2),ip),&
                         int(particles_send(3),ip),&
                         int(particles_send(4),ip),& 
                         (particles_send(ivarp),ivarp=5,nvarp_pts)
                 end if !(kfl_ptsres_binary>0) 
              end if !( lagrtyp(ilagr) % kfl_exist <= ilimi )
           end do
           call memory_deallo(mem_modul(1:2,modul),'PARTICLES_SEND','pts_output',particles_send)           
           call iofile_flush_unit(lun_resul_pts)

        end if

     end if
     
     !----------------------------------------------------------------------
     !
     ! Postprocess deposition
     !
     !----------------------------------------------------------------------

     if( kfl_ifdep == 1 ) then

        if( IMASTER ) then

           ipars = 0
           do ipart = 1,npart
              do ilagr = 1,int(num_particles_gat(ipart),ip) / nvart_pts
                 call pts_parallelization_unpack(particle,ipars,particles_gat,variables_pts,'EXIST')
                 if( particle % kfl_exist == PTS_PARTICLE_HITS_WALL .or. particle % kfl_exist == PTS_PARTICLE_OUTFLOW ) then         
                    write(lun_oudep_pts,101) (particles_gat(perm_deposition(ivard)+ipars),ivard=1,nvard_pts)
                 end if
                 ipars = ipars + nvart_pts
              end do 
           end do
           call iofile_flush_unit(lun_oudep_pts)

        else if( ISEQUEN ) then
           call memory_alloca(mem_modul(1:2,modul),'PARTICLES_SEND','pts_output',particles_send,nvard_pts)
           do ilagr = 1,mlagr
              if( lagrtyp(ilagr) % kfl_exist == PTS_PARTICLE_HITS_WALL .or. lagrtyp(ilagr) % kfl_exist == PTS_PARTICLE_OUTFLOW ) then        
                 ipars = 0
                 call pts_parallelization_pack(lagrtyp(ilagr),ipars,particles_send,deposition_var_pts)
                 write(lun_oudep_pts,101) (particles_send(ivard),ivard=1,nvard_pts)
              end if
           end do
           call memory_deallo(mem_modul(1:2,modul),'PARTICLES_SEND','pts_output',particles_send)
           call iofile_flush_unit(lun_oudep_pts)

        end if

     end if

     call ker_timeline('END_OUTPUT')

  end if

  !----------------------------------------------------------------------
  !
  ! Postprocess surface deposition 
  !
  !----------------------------------------------------------------------

  if( kfl_ifdep == 1 .and. kfl_depos_surface_pts == 1 .and. nlagr_hits_wall_pts > 0 ) then

     allocate( depos_surface(ntyla_pts+1) )
     call pts_deposition_surface(depos_surface)

     if( INOTSLAVE ) then
        write(lun_depsu_pts,101) cutim_pts,depos_surface(ntyla_pts+1),depos_surface(1:ntyla_pts)
        call iofile_flush_unit(lun_depsu_pts)
     end if

     deallocate( depos_surface )

  end if
  !
  ! Deallocate memory
  !
  call memory_deallo(mem_modul(1:2,modul),'NUM_PARTICLES_GAT' ,'pts_output',num_particles_gat)
  call memory_deallo(mem_modul(1:2,modul),'PARTICLES_SEND'    ,'pts_output',particles_send)
  call memory_deallo(mem_modul(1:2,modul),'PARTICLES_GAT'     ,'pts_output',particles_gat)
  call memory_deallo(mem_modul(1:2,modul),'PERM_POSTPROCESS'  ,'pts_output',perm_postprocess)
  call memory_deallo(mem_modul(1:2,modul),'PERM_DEPOSITION'   ,'pts_output',perm_deposition)
  call memory_deallo(mem_modul(1:2,modul),'VARIABLES_PTS'     ,'pts_output',variables_pts)
  !
  ! Formats
  !
100 format(1x,es16.8e3,3(1x,i16),400(1x,es16.8e3))
101 format(es16.8e3,400(',',es16.8e3))

end subroutine pts_output

subroutine pts_output_header()

  use def_master
  use def_partis
  use mod_outfor, only : outfor
  implicit none
  integer(ip) :: ii,ivarp
  character(255) :: binary_header=""

  if( nvarp_pts > 0 ) then
     if (kfl_ptsres_binary==0) then

        ii = 0
        call outfor(53_ip,lun_resul_pts,' ')

        do ivarp = 1,nvarp_pts
           if( postprocess_list_pts(ivarp) /= 0 ) then
              ii       = ii+1
              ioutp(1) = ii
              coutp(1) = postprocess_name_pts(postprocess_list_pts(ivarp))
              call outfor(35_ip,lun_resul_pts,' ')
           end if
        end do
        write(lun_resul_pts,1) (adjustr(trim(postprocess_name_pts(postprocess_list_pts(ivarp)))),ivarp=1,nvarp_pts)

     else !(kfl_ptsres_binary==0)

        write(binary_header,1) (adjustr(trim(postprocess_name_pts(postprocess_list_pts(ivarp)))),ivarp=1,nvarp_pts)
        write(lun_resul_pts) binary_header

     end if !(kfl_ptsres_binary==0)
  end if !( nvarp_pts > 0 ) 

1 format('#           ',a5,100('            ',a5))

end subroutine pts_output_header

subroutine pts_deposition_header()

  use def_master
  use def_partis
  implicit none

  if( nvard_pts > 0 ) then
     write(lun_oudep_pts,1) postprocess_name_pts(deposition_list_pts(1:nvard_pts))
  end if
  
1 format(100(a5,','))
  
end subroutine pts_deposition_header

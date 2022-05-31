!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    ker_adapel.f90
!> @author  Herbert Owen
!> @brief   ADAPt Exchange Location -modifies ywalb(iboun) = ywalb(iboun) * fact(itest) ; to guarantee that all exchange location points will be found.
!> @details - Adapt the distance used in exchange location so that all points are found
!!          - This will be necesary for complex geometries moreover it will give robustness in case elsest fails by slightly moving the point.
!!          - It only works with variable wall distance.
!!          - The points that are modified can be seen in Paraview. See below.
!> @} 
!-----------------------------------------------------------------------

subroutine ker_adapel()
  use def_domain

  use def_kintyp,                   only : ip,rp
  use def_master,                   only : kfl_paral, current_code,INOTMASTER, zeror, ID_TEMPER,ID_NASTIN, kfl_modul, lzone
  use def_kermod,                   only : kfl_delta, dexlo_ker
  use def_kermod,                   only : kfl_boexc_ker
  use mod_couplings,                only : COU_INIT_INTERPOLATE_POINTS_VALUES
  use mod_parall,                   only : par_code_zone_subd_to_color,PAR_MY_CODE_RANK
  use mod_parall,                   only : PAR_COMM_COLOR
  use mod_elsest,                   only : elsest_host_element
  use def_coupli,                   only : typ_color_coupling
  use def_coupli,                   only : FLOATING_WET_POINT
  use def_coupli,                   only : FLOATING_TARGET_ENTITY
  use def_coupli,                   only : ELEMENT_INTERPOLATION
  use mod_communications,           only : PAR_MIN
  use def_master,                   only : INOTMASTER
  use def_master,                   only : INOTSLAVE
  use def_master,                   only : lun_tempo
  use mod_parall,                   only : PAR_CODE_SIZE
  use mod_communications,           only : PAR_GATHER
  use mod_communications,           only : PAR_GATHERV
  use mod_communications,           only : PAR_BARRIER
  use mod_iofile,                   only : iofile
  use mod_iofile,                   only : iofile_flush_unit
  use mod_couplings_communications, only : COU_GENERATE_LOCAL_TRANSMISSION_MATRICES
  use mod_couplings_communications, only : COU_PARALLELIZE_TRANSMISSION_MATRICES
  use mod_coupling_memory,          only : cou_deallocate
  use mod_cou_output,               only : cou_output_timings
  use mod_par_output_partition,     only : par_output_coupling_timings
  use mod_messages,                 only : messages_live
  use mod_coupling_memory,          only : cou_initialization

  implicit none

  integer(ip)              :: ielem,inode,ipoin,idime,kount,ierr
  integer(ip)              :: pnode,pgaus,iboun,igaub,inodb
  integer(ip)              :: pelty,pblty,pnodb,pgaub,pmate
  real(rp)                 :: bocod(ndime,mnodb),elcod(ndime,mnode),gbcod(ndime)
  real(rp)                 :: baloc(ndime,ndime),eucta
  real(rp),pointer         :: xcoor(:,:)
  integer(ip)              :: icolo,jcolo,izone,kfl_all_found,size1,size2
  real(rp)                 :: dist_aux

  integer(ip),parameter    :: ntest = 11
  real(rp)                 :: fact(ntest)
  integer(ip),allocatable  :: missing_find(:)
  integer(ip),allocatable  :: level_found(:)
  integer(ip)              :: itest,i,ii,ipart
  type(typ_color_coupling) :: fake_waexl         ! to use COU_INIT_INTERPOLATE_POINTS_VALUES instead of
  ! I had to put it of size ntest otherwise I got 'SOMETHING STRANGE HAPPENS! because couplings can not be destroyed.

  !
  ! for visualization - Ideas borrowed from outstl
  !
  character(150)           :: fil_tempo

  real(rp),        pointer :: aux_vec(:,:)
  real(rp),        pointer :: aux_vec_gat(:)   ! in outstl xstl_gat(:)  but I belive it is cleaner to have it (:,:) - but it does not work


  integer(ip),parameter    :: n1_aux_vec = 9_ip
  integer(4)               :: nauxv4,nauxv4_tot    
  integer(4),      pointer :: nauxv4_gat(:)

  integer(ip)              :: kount_not_found_in_first
  real(rp)                 :: time1,time2
  
  kount_not_found_in_first = 0
  
  nullify(aux_vec)
  nullify(aux_vec_gat)
  nullify(nauxv4_gat)
  nauxv4 = 0 

  if ( kfl_delta /= 1 ) call runend('ker_adapel:This only works for variable wall distance')

  !  nullify(lexlo_ker)

  fact( 1) = 1.0_rp
  fact( 2) = 1.001_rp
  fact( 3) = 0.999_rp
  fact( 4) = 1.01_rp
  fact( 5) = 0.99_rp
  fact( 6) = 0.5_rp
  fact( 7) = 0.25_rp
  fact( 8) = 0.125_rp
  fact( 9) = 0.0625_rp
  fact(10) = 1.0_rp/32.0_rp
  fact(11) = 1.0_rp/64.0_rp

  allocate(missing_find(nboun))
  allocate(level_found(nboun))
  missing_find = 1_ip
  level_found  = ntest

  if ( kfl_delta /= 1 ) return

  test: do itest = 1,ntest

     nullify(xcoor)

     kount = 0_ip
     if( INOTMASTER ) then
        !
        ! Loop over boundaries - preliminary just to obtain kount
        !
        boun0: do iboun = 1,nboun

           if ( kfl_boexc_ker(kfl_codbo(iboun)) == 1 .and. missing_find(iboun) == 1 ) then
              !
              ! Element properties and dimensions
              !
              pblty = ltypb(iboun) 
              pnodb = nnode(pblty)
              ielem = lelbo(iboun)
              pgaub = ngaus(pblty) 
              pmate = 1

              if( nmate > 1 ) then
                 pmate = lmate(ielem)
              end if

              if ( pmate /= -1 )  kount = kount + pgaub

           end if

        end do boun0
        allocate( xcoor(ndime,kount) )    ! max habría que corregirlo tambien en ker_waexlo
        !        allocate( lexlo_ker(mgaub,nboun) )
     end if

     kount = 0_ip
     !
     ! Loop over boundaries
     !
     if( INOTMASTER ) then
        dist_aux= dexlo_ker
        boundaries: do iboun =1, nboun
           if ( kfl_boexc_ker(kfl_codbo(iboun)) == 1 .and. missing_find(iboun) == 1 ) then
              dist_aux = ywalb(iboun)     ! Variable wall distance (use ywalb instead of dexlo_ker)
              !
              ! Element properties and dimensions
              !
              pblty = ltypb(iboun) 
              pnodb = nnode(pblty)
              ielem = lelbo(iboun)
              pelty = ltype(ielem)
              pnode = nnode(pelty)
              pgaub = ngaus(pblty) 
              pgaus = ngaus(pelty)
              pmate = 1

              if( nmate > 1 ) then
                 pmate = lmate(ielem)
              end if

              if (  pmate /= -1  )  then   
                 !
                 ! Gather operations: ELCOD, BOCOD
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    do idime = 1,ndime
                       elcod(idime,inode) = coord(idime,ipoin)             
                    end do
                 end do

                 do inodb = 1,pnodb     ! obtain bocod for bouder
                    ipoin = lnodb(inodb,iboun)
                    do idime = 1,ndime
                       bocod(idime,inodb) = coord(idime,ipoin)
                    end do
                 end do

                 gauss_points: do igaub = 1,pgaub
                    !
                    ! Obtain normal (baloc(:,ndime) to the surface (following nsi_bouset)
                    !                 
                    call bouder(&
                         pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&    ! Cartesian derivative
                         bocod,baloc,eucta)                                   ! and Jacobian
                    call chenor(pnode,baloc,bocod,elcod)                      ! Check normal

                    kount = kount + 1_ip
                    !                    lexlo_ker(igaub,iboun) = kount

                    gbcod=0.0_rp
                    do inodb=1,pnodb
                       do idime=1,ndime
                          gbcod(idime) = gbcod(idime)        &
                               + elmar(pblty)%shape(inodb,igaub) * bocod(idime,inodb)
                       end do
                    end do
                    do idime=1,ndime
                       xcoor(idime,kount) = gbcod(idime) - fact(itest) * dist_aux * baloc(idime,ndime) 
                    end do
                 end do gauss_points

              end if

           end if

        end do boundaries
     end if
     !
     ! Intialize interpolation
     !
     izone = 0
     icolo = par_code_zone_subd_to_color(current_code,0_ip,0_ip)     ! icolo and jcolo I do not understan well 
     jcolo = par_code_zone_subd_to_color(current_code,0_ip,0_ip)     ! I just put them following nsi_velobl
     !
     ! Basic information needed to fill the coupling structures
     !
     call cou_initialization(fake_waexl)
     fake_waexl % number                    = 1002_ip                ! coupling number     
     fake_waexl % itype                     = ELEMENT_INTERPOLATION  ! Element interpolation
     fake_waexl % kfl_toda_costa            = 0_ip                   ! Let the distance change
     fake_waexl % color_target              = icolo                  ! target color
     fake_waexl % color_source              = jcolo                  ! source color
     fake_waexl % zone_source               = izone                  ! current_zone
     fake_waexl % zone_target               = izone                  ! current_zone
     fake_waexl % target_entity             = FLOATING_TARGET_ENTITY ! Where coupling is eventuall applied
     fake_waexl % wet   % npoin_wet         = kount
     fake_waexl % wet   % number_wet_points = 0_ip
     fake_waexl % wet   % point_type        = FLOATING_WET_POINT     ! Wt point is not a node
     fake_waexl % commd % PAR_COMM_WORLD    = PAR_COMM_COLOR(icolo,jcolo)

     if( associated(xcoor) ) then
        size1 = size(xcoor,1)
        size2 = size(xcoor,2)
        allocate(fake_waexl % wet % coord_wet(size1,size2))
        fake_waexl % wet % coord_wet = xcoor
     end if

     call PAR_BARRIER()
     call cputim(time1)
     call messages_live('ADAPTIVE WALL EXCHANGE PREPROCESS: COMPUTING COUPLING')
     call COU_INIT_INTERPOLATE_POINTS_VALUES(fake_waexl)
     call PAR_BARRIER()
     call cputim(time2)
     !
     ! check if everything was right
     !
     !ierr = 0
     !do ii = 1,kount
     !   if( fake_waexl % geome % status(ii) == 0 ) then
     !      print*,ii,PAR_MY_CODE_RANK,' is lost with coordinates',xcoor(:,ii),'kfl_paral',kfl_paral
     !      ierr = ierr + 1
     !   end if
     !end do
     !if ( ierr/=0 ) call runend('QUE PASA AQUI')
     !
     ! Output and postprocess
     !
     call cou_output_timings(fake_waexl)
     call par_output_coupling_timings(fake_waexl,itest-1_ip)   
     !
     ! check if exchange loc has been found for all gauss point on iboun
     !
     kount = 0_ip
     if( INOTMASTER ) then
        !
        ! Loop over boundaries - preliminary just to obtain kount
        !
        boun_check: do iboun = 1,nboun

           if ( kfl_boexc_ker(kfl_codbo(iboun)) == 1 .and. missing_find(iboun) == 1 ) then
              !
              ! Element properties and dimensions
              !
              pblty = ltypb(iboun) 
              pnodb = nnode(pblty)
              ielem = lelbo(iboun)
              pgaub = ngaus(pblty) 
              pmate = 1

              if( nmate > 1 ) then
                 pmate = lmate(ielem)
              end if

              if ( pmate /= -1 )  then
                 missing_find(iboun) = 0  ! Mark the boundary as found - if some of its gauss points is not found it will be changed to 1 in some lines
                 do igaub = 1,pgaub
                    kount = kount + 1
                    if( fake_waexl % geome % status(kount) == 0 ) then
                       if (itest == ntest) print*,' LAST TEST and point not found with coordinates',xcoor(:,kount)   !ACA probablemente mejor escribir la coord de pto de gauss que no encuentra
                       missing_find(iboun) = 1
                    end if
                 end do
                 if ( missing_find(iboun) == 0 ) then  ! If all gauss point in the boundary have been found
                    ywalb(iboun) = ywalb(iboun) * fact(itest)
                    level_found(iboun) = itest
                    if ( itest > 1 ) kount_not_found_in_first = kount_not_found_in_first + pgaub   ! to be used late for visualization
                 end if
              end if

           end if

        end do boun_check

     end if
     !
     ! Deallocate coupling
     !
     call cou_deallocate(fake_waexl)
     
     kfl_all_found = 1
     if( INOTMASTER ) then
        all_found: do iboun = 1,nboun
           if ( kfl_boexc_ker(kfl_codbo(iboun)) == 1 .and. missing_find(iboun) == 1 ) kfl_all_found = 0
        end do all_found
     end if
     call PAR_MIN(kfl_all_found,'IN MY CODE')
     if (kfl_all_found == 1) then
        if( INOTMASTER ) deallocate( xcoor )     ! me lo encontro intel inspector!!!
        exit test
     else if (itest == ntest) then
        write(*,*) 'all positions have been tested and some points have still not been found'
     end if

     if( INOTMASTER ) deallocate( xcoor )


  end do test
  deallocate(missing_find)
  !
  ! Visualization
  !
  kount = 0_ip
  if( INOTMASTER ) then
     allocate(aux_vec(n1_aux_vec,kount_not_found_in_first))
  end if

  if( INOTSLAVE ) then
     !
     ! Open file
     !
     fil_tempo = 'exchange_not_found.csv'
     call iofile(0_ip,lun_tempo,fil_tempo,'SYSTEM INFO')
     write(lun_tempo,*)'gp_x,gp_y,gp_z,found_x,found_y,found_z,level,inv_fact,kfl_paral'
     allocate( nauxv4_gat(0:PAR_CODE_SIZE-1) )
  end if  
  !
  ! How to use this in Paraview
  ! Open mesh and also open this csv file
  ! Run the filter Filters/ Alphabetical/ Table To Points.   
  ! Tell ParaView what columns are the X, Y and Z coordinate(gp_x,gp_y,gp_z). Be sure to not skip this step. Apply.
  ! Mark keep all data arrays
  ! ParaView probably didn't open up a 3d window (this is a bug). Split screen Horizontal (Icon, top right). 3D View.
  ! You can change the point size so that it looks nicer. 
  ! I believe the best value to plot is inv_fact
  !
  if( INOTMASTER ) then
     visualization: do iboun =1, nboun
        if ( kfl_boexc_ker(kfl_codbo(iboun)) == 1 .and. level_found(iboun) > 1 ) then
           !
           ! Element properties and dimensions
           !
           pblty = ltypb(iboun) 
           pnodb = nnode(pblty)
           ielem = lelbo(iboun)
           pelty = ltype(ielem)
           pnode = nnode(pelty)
           pgaub = ngaus(pblty) 
           pgaus = ngaus(pelty)
           pmate = 1

           if( nmate > 1 ) then
              pmate = lmate(ielem)
           end if

           if (  pmate /= -1  )  then   
              !
              ! Gather operations: ELCOD, BOCOD
              !
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 do idime = 1,ndime
                    elcod(idime,inode) = coord(idime,ipoin)             
                 end do
              end do

              do inodb = 1,pnodb     ! obtain bocod for bouder
                 ipoin = lnodb(inodb,iboun)
                 do idime = 1,ndime
                    bocod(idime,inodb) = coord(idime,ipoin)
                 end do
              end do

              gauss_points_vis: do igaub = 1,pgaub
                 !
                 ! Obtain normal (baloc(:,ndime) to the surface (following nsi_bouset)
                 !
                 kount = kount + 1
                 call bouder(&
                      pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&    ! Cartesian derivative
                      bocod,baloc,eucta)                                   ! and Jacobian
                 call chenor(pnode,baloc,bocod,elcod)                      ! Check normal


                 gbcod=0.0_rp
                 do inodb=1,pnodb
                    do idime=1,ndime
                       gbcod(idime) = gbcod(idime)        &
                            + elmar(pblty)%shape(inodb,igaub) * bocod(idime,inodb)
                    end do
                 end do

                 if (ndime==3) then
                    aux_vec(1:ndime,kount) = gbcod
                    aux_vec(4:ndime+3,kount) = gbcod(1:ndime) - ywalb(iboun) * baloc(1:ndime,ndime)
                    aux_vec(7,kount)   = dble(level_found(iboun))
                    aux_vec(8,kount)   = 1.0_rp/fact(level_found(iboun))
                    aux_vec(9,kount)   = dble(kfl_paral)
                 end if
              end do gauss_points_vis

           end if

        end if

     end do visualization
  end if
  !
  ! Gather all aux_vec
  ! 
  nauxv4 = int(kount_not_found_in_first,4) *  int(n1_aux_vec,4)
  call PAR_GATHER(nauxv4,nauxv4_gat,'IN MY CODE')

  if( INOTSLAVE ) then
     nauxv4_tot = 0
     do ipart = 0,PAR_CODE_SIZE-1
        nauxv4_tot = nauxv4_tot + nauxv4_gat(ipart)
     end do
     allocate( aux_vec_gat(nauxv4_tot) )
  end if
  call PAR_GATHERV(aux_vec,aux_vec_gat,nauxv4_gat,'IN MY CODE')
  !
  ! Master outputs aux_vec_gat
  !  
  if( INOTSLAVE ) then
     if( ndime == 3 ) then
        do ii = 1,nauxv4_tot/n1_aux_vec
           !        write(lun_tempo,'(9(e14.7,(a)))') ((aux_vec_gat( i + n1_aux_vec*(ii-1) ),','),i=1,9) ! gfortran did not like it this way
           ! but I do not understand why. In any case it accepts it as written below and it works fine
           write(lun_tempo,'(9(e14.7,(a)))') (aux_vec_gat( i + n1_aux_vec*(ii-1) ),',',i=1,9)
        end do
     end if
     deallocate( nauxv4_gat )
     deallocate( aux_vec_gat  )
  else
     deallocate( aux_vec )
  end if

  if( INOTSLAVE ) then
     call iofile_flush_unit(lun_tempo)
     close(lun_tempo)
  end if
  deallocate(level_found)
  !
  ! Output and postprocess
  !
  !call cou_output_timings(wallcoupling_waexl)
  !call par_output_coupling_timings(wallcoupling_waexl)
  
end subroutine ker_adapel

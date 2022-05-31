!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    ker_waexlo.f90
!> @author  Herbert Owen  & Matias Avila
!> @brief   Preliminary operations to apply the wall law at the 'Exchange Location' 
!> @details - See Bodart and Larsson - Wall modeled large eddy simulation in complex ... 2011
!>          - Obtain the exchange location for each wall law boundary gauss point
!>          - Initialize interpolation - PAR_INIT_INTERPOLATE_POINTS_VALUES
!> @} 
!-----------------------------------------------------------------------

subroutine ker_waexlo()
  
  use def_domain
  use def_kintyp,                   only : ip,rp
  use def_master,                   only : kfl_paral, current_code,INOTMASTER, zeror, ID_TEMPER,ID_NASTIN, kfl_modul, lzone, IMASTER
  use def_master,                   only : modul,mem_modul
  use def_kermod,                   only : kfl_delta, dexlo_ker, velel_ker, temel_ker, lexlo_ker, wallcoupling_waexl
  use def_kermod,                   only : kfl_waexl_ker,kfl_waexl_imp_ker,  shape_waexl
  use def_kermod,                   only : kfl_boexc_ker
  use def_coupli,                   only : ELEMENT_INTERPOLATION
  use def_coupli,                   only : FLOATING_WET_POINT
  use def_coupli,                   only : FLOATING_TARGET_ENTITY
  use mod_couplings,                only : COU_INIT_INTERPOLATE_POINTS_VALUES
  use mod_coupling_memory,          only : cou_initialization  
  use mod_parall,                   only : par_code_zone_subd_to_color,PAR_MY_CODE_RANK
  use mod_parall,                   only : PAR_COMM_COLOR
  use mod_elsest,                   only : elsest_host_element
  use def_kermod,                   only : delta_dom,ielse,relse,ndivi
  use mod_couplings_communications, only : COU_GENERATE_LOCAL_TRANSMISSION_MATRICES
  use mod_couplings_communications, only : COU_PARALLELIZE_TRANSMISSION_MATRICES
  use mod_messages,                 only : messages_live
  use mod_cou_output,               only : cou_output_timings
  use mod_par_output_partition,     only : par_output_coupling_timings
  use mod_communications,           only : PAR_BARRIER
  use mod_memory,                   only : memory_alloca
  use mod_memory,                   only : memory_deallo
  implicit none

  integer(ip)              :: ielem,inode,ipoin,idime,kount,ii,ierr
  integer(ip)              :: pnode,pgaus,iboun,igaub,inodb
  integer(ip)              :: pelty,pblty,pnodb,pgaub,pmate
  real(rp)                 :: bocod(ndime,mnodb),elcod(ndime,mnode),gbcod(ndime)
  real(rp)                 :: baloc(ndime,ndime),eucta
  integer(ip)              :: icolo,jcolo,izone,size1,size2
  real(rp)                 :: shapf(64),deriv(64*3),coloc(3), dista, dist_aux
  character(200)           :: aux_char
  real(rp)                 :: time1,time2
  integer(ip),   pointer   :: ielem_waexl(:)   
  real(rp),      pointer   :: cooel_ker(:,:)  

  if( kfl_waexl_ker /= 0_ip .and. kfl_delta == 1 ) call ker_adapel()  !  ywalb(iboun) = ywalb(iboun) * fact(itest) ; to guarantee that all exchange location points will be found.

  nullify(ielem_waexl)
  nullify(cooel_ker)

  if ( kfl_waexl_ker == 0_ip ) then 

     !allocate( velel_ker(ndime,1_ip) )
     !allocate( lexlo_ker(mgaub,1_ip) )
     !if (kfl_modul(ID_TEMPER) == 1)     allocate( temel_ker(1_ip,1_ip) )

  else

     kount = 0_ip
     if (IMASTER)then
        if (kfl_modul(ID_TEMPER) == 1)     allocate(temel_ker(1_ip,1_ip) )
     end if
     if( INOTMASTER ) then
        !
        ! Loop over boundaries - preliminary just to obtain kount
        !
        boun0: do iboun = 1,nboun

           !           if(  kfl_fixbo_nsi(iboun) ==  3 .or. &      ! Wall law
           !                & kfl_fixbo_nsi(iboun) == 13 .or. &    ! Wall law + open pressure
           !                & kfl_fixbo_nsi(iboun) == 18)  then    ! u.n in weak form
           if ( kfl_boexc_ker(kfl_codbo(iboun)) == 1 ) then
              !
              ! Element properties and dimensions
              !
              pblty = ltypb(iboun) 
              pnodb = nnode(pblty)
              ielem = lelbo(iboun)
              pgaub = ngaus(pblty) 
              pmate = lmate(ielem)

              if ( ( pmate /= -1 ) .and. ( ( delta_dom > zeror ) .or.  (kfl_delta == 1)) ) kount = kount + pgaub

           end if

        end do boun0
        !
        ! Permanent memory
        !
        call memory_alloca(mem_modul(1:2,modul),'VELEL_KER'  ,'ker_waexlo',velel_ker,ndime,kount)
        call memory_alloca(mem_modul(1:2,modul),'LEXLO_KER'  ,'ker_waexlo',lexlo_ker,mgaub,nboun)
        call memory_alloca(mem_modul(1:2,modul),'SHAPE_WAEXL','ker_waexlo',shape_waexl,mnode,kount)
        if (kfl_modul(ID_TEMPER) == 1) call memory_alloca(mem_modul(1:2,modul),'TEMEL_KER','ker_waexlo',temel_ker,1_ip,kount) 
        !
        ! Local memory
        !
        call memory_alloca(mem_modul(1:2,modul),'IELEM_WAEXL','ker_waexlo',ielem_waexl,kount)
        call memory_alloca(mem_modul(1:2,modul),'COOEL_KER'  ,'ker_waexlo',cooel_ker,ndime,kount)

     end if

     kount = 0_ip
     !
     ! Loop over boundaries
     !
     if( INOTMASTER ) then
        dist_aux= dexlo_ker
        boundaries: do iboun =1, nboun
           
           if ( kfl_boexc_ker(kfl_codbo(iboun)) == 1 ) then
              if (kfl_delta == 1) dist_aux = ywalb(iboun)     ! Variable wall distance (use ywalb instead of dexlo_ker)
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
              pmate = lmate(ielem)

              if ( ( pmate /= -1 ) .and. ( ( delta_dom > zeror ).or.( kfl_delta == 1 ) )  )  then 
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
                    lexlo_ker(igaub,iboun) = kount

                    gbcod=0.0_rp
                    do inodb=1,pnodb
                       do idime=1,ndime
                          gbcod(idime) = gbcod(idime)        &
                               + elmar(pblty)%shape(inodb,igaub) * bocod(idime,inodb)
                       end do
                    end do
                    do idime=1,ndime
                       cooel_ker(idime,kount) = gbcod(idime) - dist_aux * baloc(idime,ndime)  ! for the moment I will set dexlo_ker a const
                    end do         !  value. Later we can set some constant times de 1st elemente normal size or even more elaborate opts

                 end do gauss_points

              end if

           end if

        end do boundaries
     end if
     !
     ! Intialize interpolation
     !
     izone = 0_ip
     icolo = par_code_zone_subd_to_color(current_code,0_ip,0_ip)     ! icolo and jcolo I do not understan well 
     jcolo = par_code_zone_subd_to_color(current_code,0_ip,0_ip)     ! I just put them following nsi_velobl 
     !
     ! Basic information needed to fill the coupling structures
     ! 
     call cou_initialization(wallcoupling_waexl)
     wallcoupling_waexl % number                    = 1001_ip                ! Coupling number
     wallcoupling_waexl % itype                     = ELEMENT_INTERPOLATION  ! Element interpolation
     wallcoupling_waexl % kfl_toda_costa            = 0_ip                   ! For the moment I prefer Alya to stop so that I know what is hapening. Guillaume had put 1 force
     wallcoupling_waexl % color_target              = icolo                  ! target color   
     wallcoupling_waexl % color_source              = jcolo                  ! source color
     wallcoupling_waexl % zone_source               = izone                  ! current_zone
     wallcoupling_waexl % zone_target               = izone                  ! current_zone
     wallcoupling_waexl % target_entity             = FLOATING_TARGET_ENTITY ! Where coupling is eventuall applied
     wallcoupling_waexl % wet   % npoin_wet         = kount            
     wallcoupling_waexl % wet   % number_wet_points = 0_ip
     wallcoupling_waexl % wet   % point_type        = FLOATING_WET_POINT     ! Wet point is not a node

     wallcoupling_waexl % commd % PAR_COMM_WORLD    = PAR_COMM_COLOR(icolo,jcolo)
     !     print *, 'izone, icolo, jcolo, PAR', izone, icolo, jcolo,  PAR_COMM_COLOR(icolo,jcolo)
     if( associated(cooel_ker) ) then
        size1 = size(cooel_ker,1)
        size2 = size(cooel_ker,2)
        allocate(wallcoupling_waexl % wet % coord_wet(size1,size2))
        wallcoupling_waexl % wet % coord_wet = cooel_ker
     end if
     
     call PAR_BARRIER()
     call cputim(time1)
     call messages_live('WALL EXCHANGE PREPROCESS: COMPUTING COUPLING')
     call COU_INIT_INTERPOLATE_POINTS_VALUES(wallcoupling_waexl)
     call PAR_BARRIER()
     call cputim(time2)
     !print*,'B=',time2-time1
     !
     ! check if everything was right
     !
     ierr = 0
     do ii = 1,kount
        if( wallcoupling_waexl % geome % status(ii) == 0 ) then
           print*,ii,PAR_MY_CODE_RANK,' is lost with coordinates',cooel_ker(:,ii),'kfl_paral',kfl_paral
           ierr = ierr + 1
        end if
     end do
     if ( ierr/=0 ) call runend('QUE PASA AQUI')
     !
     ! Communication matrix creation
     !   
     call PAR_BARRIER()
     call cputim(time1)
     call messages_live('WALL EXCHANGE PREPROCESS: TRANSMISSION MATRICES')
     call COU_GENERATE_LOCAL_TRANSMISSION_MATRICES(wallcoupling_waexl)
     call PAR_BARRIER()
     call cputim(time2)
     !print*,'C=',time2-time1
     !
     ! Output and postprocess
     !
     call cou_output_timings(wallcoupling_waexl)
     call par_output_coupling_timings(wallcoupling_waexl)
     !
     !call COU_PARALLELIZE_TRANSMISSION_MATRICES(wallcoupling_waexl)
     !
     ! obtain shapes and ielem for each (iboun, igaub)
     !
     if( INOTMASTER ) then
        !
        if (kfl_waexl_imp_ker/=0) then   ! elsest does not work in parallel - for the implicit case, since you always search in the first element,
           ! there is no problem
           do ii =1, kount
              !         print *, ii, kfl_paral, 'entering'
              call elsest_host_element(&
                   ielse,relse,1_ip,meshe(ndivi),cooel_ker(:,ii),ielem,&
                   shapf,deriv,coloc,dista)
              !         print *, ii, kfl_paral, 'exiting'
              if( ielem > 0 ) then
                 pelty = abs(ltype(ielem))
                 pnode = nnode(pelty)
                 ielem_waexl(ii) = ielem
                 shape_waexl(1:pnode,ii) = shapf(1:pnode)
              else
                 write(aux_char,*)'ker_waexlo:element not found for cooel_ker(:,ii)',cooel_ker(:,ii),ii,kfl_paral
                 call runend(aux_char)
              end if
           end do
           ! if implicit exchange location, confirm point is inside first element

           do iboun =1, nboun
              !              if(  kfl_fixbo_nsi(iboun) ==  3 .or. &      ! Wall law
              !                   & kfl_fixbo_nsi(iboun) == 13 .or. &    ! Wall law + open pressure
              !                   & kfl_fixbo_nsi(iboun) == 18)  then    ! u.n in weak form
              pblty = ltypb(iboun) 
              pnodb = nnode(pblty)
              ielem = lelbo(iboun)
              do igaub = 1,pgaub
                 ii =  lexlo_ker(igaub,iboun)
                 if (ielem_waexl(ii)/=ielem) then
                    print *, 'ker_waexlo:bad_element, iproce, ielem,connec_elem:', kfl_paral,ielem_waexl(ii),ielem 
                    call runend('ker_waexlo: could not perform exch location implicit because point outside first element')     
                 end if
              end do
              !              end if
           end do
        end if
     end if
     !     print *, 'saliiii rutina', kfl_paral
     call memory_deallo(mem_modul(1:2,modul),'IELEM_WAEXL','ker_waexlo',ielem_waexl)
     call memory_deallo(mem_modul(1:2,modul),'COOEL_KER'  ,'ker_waexlo',cooel_ker)

  end if

end subroutine ker_waexlo

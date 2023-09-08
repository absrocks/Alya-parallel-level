!-----------------------------------------------------------------------
!> @addtogroup NastinInput
!> @{
!> @file    nsi_updbcs.f90
!> @date    29/01/2018
!> @author  Guillaume Houzeaux
!> @brief   Update boundary conditions
!> @details This routine updates the velocity boundary conditions:
!>          1. Before a time step begins
!>          2. Before a global iteration begins
!>          3. Before an inner iteration begins
!>          4. Force boundary conditions in UNKNO
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_updbcs(itask)

  use def_parame,                  only :  pi
  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_ker_space_time_function
  use mod_communications,          only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_couplings,               only : COU_PUT_VALUE_ON_TARGET
  use mod_nsi_schur_operations,    only : nsi_rotsch
  use mod_projec,                  only : projec_mass_conservation
  use mod_memory,                  only : memory_alloca
  use mod_memory,                  only : memory_alloca_min
  use mod_memory,                  only : memory_deallo
  use mod_output,                  only : output_mesh_gid_format
  use mod_couplings,               only : couplings_impose_dirichlet
  use mod_std
  use mod_ker_subdomain,           only : ker_subdomain_motion_exists

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,idime,iboun,dummi,itotv,izdom,jpoin
  integer(ip)             :: inodb,inode,ifunc,ibopo,jtask,jdime,kpoin
  integer(ip)             :: iperi,izone,iflow,kk,ierro,isubd,ielem
  integer(ip)             :: i_left,i_right,ind_left,ind_right
  real(rp)                :: venew,veold,cuti2,vefun(3)
  real(rp)                :: rot_matrix(ndime,ndime),udotn,given_exnor(3)
  real(rp)                :: bvess_new(3),xnorm,relax,bvnew,bvold
  real(rp)                :: xx,alpha,ulocnorm,center_gravity(3)
  real(rp)                :: flow_rate_stf_value
  real(rp),    external   :: funcre
  logical(lg)             :: update_flow_rate
  logical(lg), pointer    :: lboun_nsi(:)
  logical(lg), pointer    :: lpoin_nsi(:)
  real(rp),    pointer    :: bvess_loc(:,:)

  nullify(lboun_nsi)
  nullify(lpoin_nsi)
  nullify(bvess_loc)

  jtask = abs(itask)
  
  select case( jtask )

  case( 0_ip )

     !-------------------------------------------------------------------
     !
     ! At the beginning of the run
     !
     !-------------------------------------------------------------------

     izone = lzone(ID_NASTIN)

     if( kfl_confi_nsi == 0 ) then 
        !
        ! Prescribe pressure: Automatic
        !    
        if( INOTMASTER ) then 
           inode = huge(ip)
           do ipoin = 1,npoin
              if( lpoty(ipoin) /= 0 ) then
                 if( lninv_loc(ipoin) < inode ) inode = lninv_loc(ipoin)
              end if
           end do
        end if
        call parari('MIN',0_ip,1_ip,inode)
        if( INOTMASTER ) then 
           nodpr_nsi = 0 
           do ipoin = 1,npoin
              if( lninv_loc(ipoin) == inode ) then
                 nodpr_nsi = ipoin
              end if
           end do
        end if

     else if( kfl_confi_nsi > 0 ) then 
        !
        ! Prescribe pressure: On node
        ! If MM is used, the numbering refers to the MM(0) level, that is the original one
        !     
        if( INOTMASTER ) then 
           inode     = nodpr_nsi
           nodpr_nsi = 0
           kpoin     = 0
           do while( kpoin < npoin )
              kpoin = kpoin + 1
              ipoin = kpoin
              if( lninv_loc(ipoin) == inode ) then
                 nodpr_nsi = ipoin
                 kpoin = npoin
              end if
           end do
        end if

     else if( kfl_confi_nsi == -1 ) then 
        !
        ! Automatic bc for pressure: Prescribe when velocity is free
        !     
        if( INOTMASTER  ) then
           do ipoin = 1,npoin
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 if( kfl_fixpr_nsi(1,ipoin) == 0 ) then
                    dummi = 0
                    do idime = 1,ndime
                       if( kfl_fixno_nsi(idime,ipoin) == 0 ) dummi = dummi + 1
                    end do
                    if( dummi == ndime ) kfl_fixpr_nsi(1,ipoin)=1  
                    !! You should also check above that the node is not periodic
                 end if
              end if
           end do
        end if
     end if
     !
     ! Prescribed node cannot be a periodic slave node
     !
     if( kfl_confi_nsi >= 0 .and. INOTMASTER ) then
        if( kfl_matdi_nsi /= 0 ) then
           do iperi = 1,nperi
              if( lperi(2,iperi) == nodpr_nsi ) then
                 nodpr_nsi = lperi(1,iperi)
              end if
           end do
        else
           do iperi = 1,nperi        
              if( lperi(2,iperi) == nodpr_nsi .or. lperi(1,iperi) == nodpr_nsi ) then
                 call runend('NSI_UPDBCS: DO NOT PRESCRIBE PRESSURE ON A MASTER OR SLAVE NODE')
              end if
           end do
        end if
     end if
     !
     ! User modifications to boundary conditions
     !
     if( INOTMASTER  ) then
        call nsi_modbcs()
     end if
     !
     ! Modify KFL_FIXPR_NSI if flow is confined
     !
     if( INOTMASTER ) then
        if( kfl_confi_nsi >= 0 .and. nodpr_nsi > 0 ) then
           ibopo = lpoty(nodpr_nsi)
           !           if( ibopo == 0 ) call runend('NSI_UPDBCS: PRESSURE IS PRESCRIBED ON INTERIOR NODE')
           kfl_fixpr_nsi(1,nodpr_nsi) = 1
        end if
     end if
     !
     ! Pressure Dirichlet
     !
     if( nodpr_nsi > 0 .and. INOTMASTER ) then
        bpess_nsi(1,nodpr_nsi,1)   = valpr_nsi
        kfl_fixpp_nsi(1,nodpr_nsi) = 1
     end if
     ! 
     ! Bemol: modify if using Parall
     !
     if( bemol_nsi > 0.0_rp .and. INOTMASTER ) then
        do iboun = 1,nboun
           if( kfl_fixbo_nsi(iboun) == 13 ) then
              do inodb = 1,nnode(ltypb(iboun))
                 ipoin = lnodb(inodb,iboun)
                 kfl_fixpr_nsi(1,ipoin) = 1
              end do
           end if
        end do
     end if
     !
     ! No-penetration in weak form
     !
     if( INOTMASTER ) then
        call memgen(1_ip,npoin,0_ip)
        do iboun = 1,nboun
           if( kfl_fixbo_nsi(iboun) == 18 ) then
              do inodb = 1,nnode(abs(ltypb(iboun)))
                 ipoin = lnodb(inodb,iboun)
                 gisca(ipoin) = 1
              end do
           end if
        end do
        call parari('SLX',NPOIN_TYPE,npoin,gisca)
        do ipoin = 1,npoin
           if( gisca(ipoin) /= 0 ) then
              kfl_fixpr_nsi(1,ipoin) = 0
           end if
        end do
        call memgen(3_ip,npoin,0_ip)
     end if

     if( INOTMASTER.and.exfpr_nsi == 1 ) then
        write (*,*) 'CONDITION EXFPR, EXTENDING FIXPR NODES '
        call memgen(1_ip,npoin,0_ip)
        do ipoin = 1,npoin
           if( kfl_fixpr_nsi(1,ipoin) == 1 ) then
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)
                 if( kfl_fixpr_nsi(1,jpoin) == 0 ) gisca(jpoin) = 1
              end do
           end if
        end do
        call parari('SLX',NPOIN_TYPE,npoin,gisca)
        do ipoin = 1,npoin
           if( gisca(ipoin) /= 0 .and. lpoty(ipoin)/=0) then
              kfl_fixpr_nsi(1,ipoin) = 1
           end if
        end do
        call memgen(3_ip,npoin,0_ip)
     end if
     !
     ! Check outflows with too few nodes
     !
     call nsi_detect_outflows()

  case( 1_ip )

     !-------------------------------------------------------------------
     !  
     ! Before a time step
     !     
     !-------------------------------------------------------------------

     if( INOTMASTER .and. kfl_conbc_nsi == 0 ) then
        !
        ! Space/Time functions
        !
        if( number_space_time_function > 0 ) then
           do ipoin = 1,npoin
              if( kfl_funno_nsi(ipoin) < 0 ) then 
                 ifunc = -kfl_funno_nsi(ipoin)            
                 call ker_space_time_function(&
                      ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,vefun(1:ndime))
                 do idime = 1,ndime
                    vefun(idime) = vefun(idime) * bvess_nsi(idime,ipoin,2)
                 end do
                 if( kfl_timei_nsi /= 0 .and. kfl_tiacc_nsi == 2 .and. kfl_tisch_nsi == 1 ) then
                    do idime = 1,ndime
                       veold = veloc(idime,ipoin,ncomp_nsi)
                       bvess_nsi(idime,ipoin,1) = 0.50_rp*(vefun(idime)+veold)
                    end do
                 else
                    do idime = 1,ndime
                       bvess_nsi(idime,ipoin,1) = vefun(idime)
                    end do
                 end if
              end if
           end do
        end if
        !
        ! Transient fields
        !
        if( kexist_tran_fiel > 0 ) then
           do ipoin = 1,npoin
              !
              ! This is a bit "male trout" (trucho): codes are read in nsi.dat and when bc's are associated to a field, the field number
              ! is assigned to kfl_funno_nsi(ipoin) but adding 1000 to it. as interval_funno = 1000, this is the trigger to know if a field
              ! is associated with the bc. this should be made more elegant in the future.
              !
              if( kfl_funno_nsi(ipoin) > interval_funno ) then 
                 ifunc = kfl_funno_nsi(ipoin) - interval_funno
                 kk = k_tran_fiel(ifunc)
                 xx = x_tran_fiel(ifunc)
                 do idime = 1,ndime
                    vefun(idime) = xfiel(ifunc) % a(idime,ipoin,kk) * xx + &
                         xfiel(ifunc) % a(idime,ipoin,kk+1) * (1.0_rp-xx)
                 end do
                 !
                 ! These lines are identical to the ones from space time function perhaps it would be nicer to create a small subroutine  
                 !
                 if( kfl_timei_nsi /= 0 .and. kfl_tiacc_nsi == 2 .and. kfl_tisch_nsi == 1 ) then
                    do idime = 1,ndime
                       veold = veloc(idime,ipoin,ncomp_nsi)
                       bvess_nsi(idime,ipoin,1) = 0.50_rp*(vefun(idime)+veold)
                    end do
                 else
                    do idime = 1,ndime
                       bvess_nsi(idime,ipoin,1) = vefun(idime)
                    end do
                 end if
              end if
           end do
        end if
        !
        ! Nodes with fixity 7 (open or close uses) bvess_nsi which has already been updated if it comes from a transient field
        !
        if( kfl_exist_fixi7_nsi==1 .and. itask==1 ) then    

           call memgen(1_ip,npoin,0_ip)   !allocate gisca
           call open_close(1,kfl_fixno_nsi,bvess_nsi,ndime)

           do ipoin = 1,npoin
              if( abs(kfl_fixno_nsi(1,ipoin)) == 7 ) then
                 if ( kfl_fixrs_nsi(lpoty(ipoin)) /= 0 ) then
                    print*,'ipoinwith skew',lninv_loc(ipoin)
                    call runend('nsi_updbcs:ipoinwith skew')
                 end if
                 !
                 ! Open of close point - the decision has already been taken in open_close and is stored in gisca
                 !
                 if( gisca(ipoin) > 0 ) then
                    kfl_fixno_nsi(1:ndime-1,ipoin) = 7
                    veloc(1:ndime-1,ipoin,2) = bvess_nsi(1:ndime-1,ipoin,1)  ! correct initial guess for outer iterations  - x & y components
                    if( abs(kfl_fixno_nsi(ndime,ipoin)) == 7 ) then          ! z component
                       kfl_fixno_nsi(ndime,ipoin) = 7 !if z component is fixed(771) it will not be modified
                       veloc(ndime,ipoin,2) = bvess_nsi(ndime,ipoin,1)  ! correct initial guess for outer iterations - this will be used by tem & tur_updbcs
                    end if
                    kfl_fixpr_nsi(1,ipoin) = 0    ! for fixpr I don not care to use a 7 as it is what happens with velocity that decides

                    veloc(1:ndime,ipoin,2) = bvess_nsi(1:ndime,ipoin,1)    ! not sure if it is ok for CN
                 else
                    kfl_fixno_nsi(1:ndime-1,ipoin) = -7    
                    if( abs(kfl_fixno_nsi(ndime,ipoin)) == 7 ) kfl_fixno_nsi(ndime,ipoin) = -7
                    kfl_fixpr_nsi(1,ipoin) = 1
                 end if
              end if
           end do
           call memgen(3_ip,npoin,0_ip)   !deallocate gisca

        end if
        !
        ! BVESS: Dirichlet time dependent condition, time functions
        !
        do ipoin = 1,npoin
           ifunc = kfl_funno_nsi(ipoin)

           if( kfl_fixno_nsi(1,ipoin) == 9 ) then
              !
              ! U(t) = U_inf - w x r
              !
              bvess_nsi(:,ipoin,1)=bvess_nsi(:,ipoin,2)
              if(ndime==2) then
                 bvess_nsi(1,ipoin,1) = bvess_nsi(1,ipoin,1) &
                      + fvela_nsi(3) * coord(2,ipoin)
                 bvess_nsi(2,ipoin,1) = bvess_nsi(2,ipoin,1) &
                      - fvela_nsi(3) * coord(1,ipoin)
              else
                 bvess_nsi(1,ipoin,1) = bvess_nsi(1,ipoin,1) &
                      + fvela_nsi(3) * coord(2,ipoin) &
                      - fvela_nsi(2) * coord(3,ipoin)
                 bvess_nsi(2,ipoin,1) = bvess_nsi(2,ipoin,1) &
                      - fvela_nsi(3) * coord(1,ipoin) &
                      + fvela_nsi(1) * coord(3,ipoin)
                 bvess_nsi(3,ipoin,1) = bvess_nsi(3,ipoin,1) &
                      - fvela_nsi(1) * coord(2,ipoin) &
                      + fvela_nsi(2) * coord(1,ipoin)        
              end if

           else if( abs(kfl_fixno_nsi(1,ipoin)) == 8 ) then
              !
              !             R(theta)
              ! No inertial ========> Lab
              !
              !             R(-theta)
              !         Lab ========> No inertial
              !
              ! u(t) = R(-theta).U_inf - w x r
              ! U(t) = R(theta).( u + w x r ) (used in INERTIAL postprocess, to be used with postprocess of Mesh rotation)
              !
              call rotmat(ndime,-cutim*fvnoa_nsi,fvdia_nsi,rot_matrix)
              do idime = 1,ndime
                 bvess_new(idime) = 0.0_rp
                 do jdime = 1,ndime
                    bvess_new(idime) = bvess_new(idime) + rot_matrix(idime,jdime)*bvess_nsi(jdime,ipoin,2)
                 end do
              end do
              if( ndime == 2 ) then
                 bvess_new(1) = bvess_new(1) &
                      + fvela_nsi(3) * coord(2,ipoin)
                 bvess_new(2) = bvess_new(2) &
                      - fvela_nsi(3) * coord(1,ipoin)
              else
                 bvess_new(1) = bvess_new(1) &
                      + fvela_nsi(3) * coord(2,ipoin) &
                      - fvela_nsi(2) * coord(3,ipoin)
                 bvess_new(2) = bvess_new(2) &
                      - fvela_nsi(3) * coord(1,ipoin) &
                      + fvela_nsi(1) * coord(3,ipoin)
                 bvess_new(3) = bvess_new(3) &
                      - fvela_nsi(1) * coord(2,ipoin) &
                      + fvela_nsi(2) * coord(1,ipoin)        
              end if

              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 udotn = 0.0_rp
                 xnorm = 0.0_rp
                 do idime = 1,ndime
                    xnorm = xnorm + bvess_new(idime) * bvess_new(idime)
                 end do
                 if( xnorm > zeror ) xnorm = 1.0_rp / sqrt(xnorm)
                 do idime = 1,ndime
                    udotn = udotn + exnor(idime,1,ibopo) * bvess_new(idime) 
                 end do
                 udotn = udotn * xnorm
                 !
                 udotn = -1.0_rp                  !  All Dirichlet !!!!
                 !
                 if( udotn <= 0.1_rp ) then       !  Inflow
                    if( kfl_fixno_nsi(1,ipoin) == -8 ) then
                       relax = 0.1_rp
                       do idime = 1,ndime
                          kfl_fixno_nsi(idime,ipoin) =  8
                          kfl_fixpr_nsi(1,ipoin)     =  0
                          bvess_nsi(idime,ipoin,1)   =  relax * bvess_new(idime) + (1.0_rp-relax) * bvess_nsi(idime,ipoin,1)
                       end do
                    else
                       do idime = 1,ndime
                          kfl_fixno_nsi(idime,ipoin) =  8
                          kfl_fixpr_nsi(1,ipoin)     =  0
                          bvess_nsi(idime,ipoin,1)   =  bvess_new(idime)
                       end do
                    end if
                 else                             !  Outflow
                    do idime = 1,ndime
                       kfl_fixno_nsi(idime,ipoin) = -8
                       kfl_fixpr_nsi(1,ipoin)     =  1
                    end do
                 end if
              end if

           else if( ifunc > 0 .and. ifunc <= interval_funno ) then
              !
              ! Function: U(t)=f(t)*U(0)
              !
              do idime = 1,ndime
                 if( kfl_fixno_nsi(idime,ipoin) == 1 ) then     
                    venew = bvess_nsi(idime,ipoin,2) * funcre( &
                         time_function(ifunc) % parameters,    &
                         time_function(ifunc) % npara,         &
                         time_function(ifunc) % kfl_type,      &
                         cutim)
                    if( kfl_timei_nsi /= 0 .and. kfl_tiacc_nsi == 2 .and. kfl_tisch_nsi == 1 ) then
                       veold = veloc(idime,ipoin,ncomp_nsi)
                       bvess_nsi(idime,ipoin,1) = 0.50_rp*(venew+veold)
                    else
                       bvess_nsi(idime,ipoin,1) = venew                 
                    end if
                 end if
              end do

           else if( kfl_fixno_nsi(1,ipoin) == 5 ) then
              !
              ! SEM: synthetic eddy method
              !
              !    bvess_nsi(1,ipoin,1) = 3.7_rp
              !    do idime = 2,ndime
              !       bvess_nsi(idime,ipoin,1) = 0.0_rp                 
              !    end do             

           else if( kfl_fixno_nsi(1,ipoin) == 6 ) then
              !
              ! SPACE-TIME boundary condition read from file
              !
              !if (ittim > 0 ) then 

              !   ipass = int( mod(cutim,nbtdt_nsi*(nbtim_nsi-1_ip)) / nbtdt_nsi ,ip)   ! Location of the plane
              !   alpha = mod(cutim,nbtdt_nsi) / nbtdt_nsi                              ! Interpolation between planes

              !   do idime = 1,ndime
              !      kfl_fixno_nsi(idime,ipoin) = 6
              !      bvess_nsi(idime,ipoin,1) =  alpha * bnval_nsi(ipass*nbval_nsi+iboun_nsi(ipoin),idime) &
              !           + (1.0_rp - alpha) * bnval_nsi((ipass + 1_ip)*nbval_nsi+iboun_nsi(ipoin),idime)
              !   end do
              !else
              !   do idime = 1,ndime
              !      kfl_fixno_nsi(idime,ipoin) = 6
              !      bvess_nsi(idime,ipoin,1) = bnval_nsi(iboun_nsi(ipoin),idime)
              !   end do
              !end if


              if (ittim > 0 ) then
                 alpha = mod(cutim, nbtdt_nsi) / nbtdt_nsi                ! weight of right value  (mod(t,dt)/dt) 
                 i_left = mod( int( cutim / nbtdt_nsi, ip), nbtim_nsi)     ! Index of left value (0..nbtim_nsi-1)
                 if (i_left /= (nbtim_nsi-1_ip) ) then                             ! Index of right value (0..nbtim_nsi-1)
                    i_right = i_left + 1_ip
                 else
                    i_right = 0_ip
                 end if

                 !
                 ! indeces in table e.g.:
                 ! 
                 ind_left      = i_left  * nbval_nsi + iboun_nsi(ipoin) 
                 ind_right     = i_right * nbval_nsi + iboun_nsi(ipoin)  

                 do idime = 1,ndime
                    kfl_fixno_nsi(idime,ipoin) = 6
                    bvess_nsi(idime,ipoin,1) =  (1.0_rp-alpha) * bnval_nsi(ind_left ,idime) &
                         +  alpha          * bnval_nsi(ind_right,idime)
                 end do
              else
                 do idime = 1,ndime
                    kfl_fixno_nsi(idime,ipoin) = 6
                    bvess_nsi(idime,ipoin,1) = bnval_nsi(iboun_nsi(ipoin)  ,idime)
                 end do
              end if
           end if
        end do
        !
        ! BVNAT: Neumann space/time or time functions
        !
        if( itask < 0 ) then
           cuti2 = cutim + dtime
        else
           cuti2 = cutim
        end if
        do iboun = 1,nboun
           if( kfl_fixbo_nsi(iboun) /= 0 ) then
              center_gravity = 0.0_rp
              do inodb= 1,nnode(ltypb(iboun))
                 center_gravity(1:ndime) = center_gravity(1:ndime) + coord(1:ndime,lnodb(inodb,iboun))/nnode(ltypb(iboun))
              end do
              ifunc = kfl_funbo_nsi(iboun)
              if( ifunc > 0 ) then
                 venew = bvnat_nsi(1,iboun,2) *funcre(   &
                      time_function(ifunc) % parameters, & 
                      time_function(ifunc) % npara,      &
                      time_function(ifunc) % kfl_type,   &
                      cuti2)
                 if( kfl_timei_nsi /= 0 .and. kfl_tiacc_nsi == 2 .and. kfl_tisch_nsi == 1 ) then
                    bvnat_nsi(1,iboun,1) = 0.50_rp * ( venew + bvnat_nsi(1,iboun,1) )
                 else
                    bvnat_nsi(1,iboun,1) = venew                 
                 end if
              else if( ifunc < 0 ) then
                 ifunc = -kfl_funbo_nsi(iboun)   
                 call ker_space_time_function(&
                      ifunc,center_gravity(1),center_gravity(2),center_gravity(ndime),cutim,bvnew)
                 bvnew = bvnew * bvnat_nsi(1,iboun,2)
                 if( kfl_timei_nsi /= 0) then
                    if ( (kfl_tiacc_nsi == 2 .and. kfl_tisch_nsi == 1) .or. (kfl_tiacc_nsi == 4 .and. kfl_tisch_nsi == 4) ) then
                       bvold = bvnat_nsi(1,iboun,1)
                       bvnat_nsi(1,iboun,1) = 0.50_rp*(bvnew+bvold)
                    else
                       bvnat_nsi(1,iboun,1) = bvnew
                    end if
                 end if
              else if (kfl_cadan_nsi == 1) then

                 if (kfl_fixbo_nsi(iboun) == 2) bvnat_nsi(1,iboun,1) = press_cadan_nsi

              end if
           end if
        end do
        !
        ! Mesh velocity with code 3
        !
        if( ker_subdomain_motion_exists() ) then
           do ipoin = 1,npoin
              do idime = 1,ndime
                 if( kfl_fixno_nsi(idime,ipoin) == 3 ) then
                   bvess_nsi(idime,ipoin,1) = velom(idime,ipoin)
                 end if
              end do
           end do
        end if
        !
        ! Impose boundary conditions on VELOC
        !
        if( INOTMASTER ) then
           call nsi_rotsch(1_ip,veloc(:,:,2)) ! Global to local
           do ipoin = 1,npoin
              do idime = 1,ndime
                 if( kfl_fixno_nsi(idime,ipoin) > 0 ) then
                    veloc(idime,ipoin,2) = bvess_nsi(idime,ipoin,1)
                 end if
              end do
           end do
           call nsi_rotsch(2_ip,veloc(:,:,2)) ! Local to global

        end if
        !
        ! Exact solution
        !
        if( kfl_exacs_nsi /= 0 ) call nsi_exaerr(1_ip)
     end if
     !
     ! Dirichlet coupling
     !
     call couplings_impose_dirichlet(solve(1),veloc(:,:,2))
     call couplings_impose_dirichlet(solve(2),press(:,2))
     !
     ! Flow rate constraint
     !
     if( maxval(kfl_flow_rate_codes_nsi(1:mflow_nsi)) > 0 ) then

        if( INOTMASTER ) then
           call memory_alloca(mem_modul(1:2,modul),'BVESS_LOC','nsi_updbcs',bvess_loc,ndime,npoin)
           call memory_alloca(mem_modul(1:2,modul),'LBOUN_NSI','nsi_updbcs',lboun_nsi,nboun)
           call memory_alloca(mem_modul(1:2,modul),'LPOIN_NSI','nsi_updbcs',lpoin_nsi,npoin)
        else
           call memory_alloca_min(bvess_loc)
           call memory_alloca_min(lboun_nsi)
           call memory_alloca_min(lpoin_nsi)           
        end if

        do iflow = 1,size(kfl_flow_rate_codes_nsi,KIND=ip)

           update_flow_rate = .false.
           if( ittim == 0 .and. kfl_flow_rate_codes_nsi(iflow) > 0 ) update_flow_rate = .true.
           if( kfl_flow_rate_stfn_nsi(iflow) > 0 ) update_flow_rate = .true. 

           if( update_flow_rate ) then

              flow_rate_stf_value = 1.0_rp

              if( INOTMASTER ) then

                 lpoin_nsi      = .false.
                 bvess_loc(:,:) = bvess_nsi(:,:,2)
                 
                 if( nboun > 0 ) then 
                    lboun_nsi = .false.
                    do iboun = 1,nboun
                       if( kfl_codbo(iboun) == kfl_flow_rate_codes_nsi(iflow) ) lboun_nsi(iboun) = .true.
                    end do
                 end if
                 
                 do ipoin = 1,npoin
                    if(  kfl_codno(1,ipoin) == kfl_flow_rate_codes_nsi(iflow) .and. &
                         kfl_codno(2,ipoin) == mcodb+1 ) then
                       lpoin_nsi(ipoin) = .true.
                    end if
                 end do
                 
                 if( kfl_flow_rate_stfn_nsi(iflow) > 0 ) then
                    ifunc = kfl_flow_rate_stfn_nsi(iflow)
                    ipoin = 1
                    call ker_space_time_function(&
                         ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,vefun(1:ndime))
                    flow_rate_stf_value = vefun(1)
                 end if
                 !if (kfl_flow_rate_normal_nsi(iflow) > 0) then
                 !   given_exnor= flow_rate_normals_nsi(iflow,1:ndime)
                 !   call nsi_rotsch(2_ip,bvess_loc,given_exnor) ! Local to global, when normal defined
                 !else
                    call nsi_rotsch(2_ip,bvess_loc) ! Local to global                    
                 !end if

              end if

              !if (kfl_flow_rate_normal_nsi(iflow) > 0) then
              !   given_exnor= flow_rate_normals_nsi(iflow,1:ndime)
              !   call projec_mass_conservation(&
              !        bvess_loc,lboun_nsi,lpoin_nsi,'LOCAL MASS',&
              !        flow_rate_values_nsi(iflow)*flow_rate_stf_value,&
              !        'IN MY CODE WITHOUT MASTER',given_exnor,ERROR=ierro)
              !else
                 call projec_mass_conservation(&
                      bvess_loc,lboun_nsi,lpoin_nsi,'LOCAL MASS',&
                      flow_rate_values_nsi(iflow)*flow_rate_stf_value,&
                      'IN MY CODE WITHOUT MASTER',ERROR=ierro)
              !end if

              if( INOTMASTER ) then
                 if( ierro == 0 ) then    
                    do ipoin = 1,npoin
                       if( lpoin_nsi(ipoin) )  veloc(:,ipoin,2) = bvess_loc(:,ipoin)
                       !
                       !                 Leave this commented:
                       !                       ulocnorm= 0.0_rp
                       !                       do idime=1,ndime
                       !                          ulocnorm= ulocnorm + bvess_loc(idime,ipoin)*bvess_loc(idime,ipoin)
                       !                       end do
                       !                       ulocnorm= sqrt(ulocnorm)
                       !                       if (ulocnorm == 0.0_rp) kfl_fixno_nsi(1:ndime,ipoin) = 1_ip 
                       !
                    end do
                 else
                    call output_mesh_gid_format(meshe(ndivi),'FLOW_RATE_PROBLEM_'//trim(intost(iflow)),lun_outpu,lpoin_nsi)                    
                    write(*,*) 'PROBLEMMMMMMMMMMMMMMMMMMMMMMMMMM-nsi_updbcs'
                    write(lun_outpu,*) 'ERROR IN FLOW RATE'
                    write(lun_outpu,*) iflow,flow_rate_values_nsi(iflow),flow_rate_stf_value
                    if( associated(lboun_nsi) ) write(lun_outpu,*) 'LBOUN= ',count(lboun_nsi,KIND=ip)
                    if( associated(lpoin_nsi) ) write(lun_outpu,*) 'LPOIN= ',count(lpoin_nsi,KIND=ip)
                 end if

                 call nsi_rotsch(1_ip,bvess_loc) ! Global to local

                 if( ierro == 0 ) then
                    do ipoin = 1,npoin
                       if( lpoin_nsi(ipoin) ) bvess_nsi(:,ipoin,1) = bvess_loc(:,ipoin)
                    end do
                 end if

              end if

           end if
        end do

        call memory_deallo(mem_modul(1:2,modul),'LBOUN_NSI','nsi_updbcs',lboun_nsi)
        call memory_deallo(mem_modul(1:2,modul),'LPOIN_NSI','nsi_updbcs',lpoin_nsi)             
        call memory_deallo(mem_modul(1:2,modul),'BVESS_LOC','nsi_updbcs',bvess_loc)

     end if

  case( 2_ip )

     !-------------------------------------------------------------------
     !
     ! Before a global iteration
     !  
     !-------------------------------------------------------------------


  case( 3_ip )

     !-------------------------------------------------------------------
     !
     ! Before an inner iteration
     ! 
     !-------------------------------------------------------------------

     if( INOTMASTER ) then
     end if

  case(4_ip)

     !-------------------------------------------------------------------
     !
     ! Enforce boundary conditions
     !
     !-------------------------------------------------------------------

     if( NSI_MONOLITHIC ) then
        call nsi_rotunk(1_ip,unkno) ! Gobal to local
        do ipoin = 1,npoin
           itotv = (ipoin-1)*(ndime+1)
           do idime = 1,ndime
              itotv = itotv + 1
              if(  kfl_fixno_nsi(idime,ipoin) == 1 .or. &
                   kfl_fixno_nsi(idime,ipoin) == 8 .or. &
                   kfl_fixno_nsi(idime,ipoin) == 9 .or. &
                   kfl_fixno_nsi(idime,ipoin) == 6   ) then
                 unkno(itotv) = bvess_nsi(idime,ipoin,1)
              end if
           end do
        end do
        call nsi_rotunk(2_ip,unkno) ! Local to global
     else
        call nsi_rotsch(1_ip,unkno) ! Global to local
        itotv = 0
        do ipoin = 1,npoin
           do idime = 1,ndime
              itotv = itotv + 1
              if(  kfl_fixno_nsi(idime,ipoin) == 1 .or. &
                   kfl_fixno_nsi(idime,ipoin) == 8 .or. &
                   kfl_fixno_nsi(idime,ipoin) == 9 .or. &
                   kfl_fixno_nsi(idime,ipoin) == 6   ) then
                 unkno(itotv) = bvess_nsi(idime,ipoin,1)
              end if
           end do
        end do
        call nsi_rotsch(2_ip,unkno) ! Local to global
     end if

  end select

end subroutine nsi_updbcs




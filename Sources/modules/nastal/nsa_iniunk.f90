!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_iniunk.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Set up the initial condition for the problem variables
!> @details Set up the initial condition for the problem variables
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_iniunk
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_iniunk
  ! NAME 
  !    nsa_iniunk
  ! DESCRIPTION
  !    This routine 
  ! USED BY
  !    nsa_begste
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_parame
  use      def_nastal
  use      def_kermod
  use      mod_ker_proper
  use mod_commdom_alya, only: INONE
  implicit none
  integer(ip)        :: idime,icomp,ipoin,kpoin,itime
  real(rp)           :: velmi,xrano(3)
  real(rp)           :: dummr,rgasc,xhecv
  integer(ip)        :: dummi

  kfl_stead_nsa = 0
  call nsa_updbcs(zero)

  avtim_nsa = 0.0_rp  ! Accumulated time for time-averaging variables
  avvel_nsa = 0.0_rp

  ! Compute mean velocities only if total pressure conditions are present
  icomp=min(TIME_N,ncomp_nsa)
  if (nuinlet_nsa > 0) call nsa_evaluate_uinlet(1_ip,icomp)

  if( kfl_rstar /= 0 ) then

     conve_nsa(1) = 1.0_rp
     conve_nsa(2) = 0.0_rp
     diffu_nsa    = 0.0_rp
     react_nsa    = 0.0_rp

     umosg_nsa = 0.0_rp
     densg_nsa = 0.0_rp
     enesg_nsa = 0.0_rp

     !
     ! Read restart file
     !
     call nsa_restar(1_ip)

     avtim_nsa = cutim ! Accumulated time for time-averaging variables

     if( INOTMASTER ) then

        do ipoin = 1,npoin

           do idime=1,ndime
              veloc(idime,ipoin,icomp)= umome(idime,ipoin,icomp) / densi(ipoin,icomp) 
           end do

           velmi = 0.0_rp

           do idime=1,ndime
              velmi = velmi + &
                   veloc(idime,ipoin,icomp) * veloc(idime,ipoin,icomp) 
           end do
           velmi = sqrt(velmi)

           tempe(ipoin,icomp) = (energ(ipoin,icomp)/densi(ipoin,icomp) - 0.5_rp*velmi*velmi) / cvcoe_nsa

           call nsa_stalaw(2_ip,0_ip,densi(ipoin,icomp),press(ipoin,icomp),&
                tempe(ipoin,icomp),coord(ivert_nsa,ipoin),dummr,dummr,cpcoe_nsa)        ! compute p from rho and T given
        end do

        !     do itime=2+kfl_tiacc_nsa,4,-1
        !           do ipoin=1,npoin
        !              veloc(1:ndime,ipoin,itime) = veloc(1:ndime,ipoin,itime-1)
        !              umome(1:ndime,ipoin,itime) = umome(1:ndime,ipoin,itime-1)
        !              press(        ipoin,itime) = press(        ipoin,itime-1)
        !              densi(ipoin,itime) = densi(ipoin,itime-1)
        !              tempe(ipoin,itime) = tempe(ipoin,itime-1)
        !              energ(ipoin,itime) = energ(ipoin,itime-1)
        !           end do
        !     end do
     end if

  else 

     if( INOTMASTER ) then

        ! Convection-Diffusion-Reaction Equation
        conve_nsa(1) = 1.0_rp
        conve_nsa(2) = 0.0_rp
        diffu_nsa    = 0.0_rp
        react_nsa    = 0.0_rp

        ! Subscale vectors
        umosg_nsa = 0.0_rp
        densg_nsa = 0.0_rp
        enesg_nsa = 0.0_rp

        if (kfl_inifi_nsa(1) == 0 .or. kfl_inifi_nsa(2) == 0 .or. kfl_inifi_nsa(3) == 0 ) then
           !
           ! Initialize variables from reference values
           !        
           if (kfl_infun_nsa > 0) then
              !
              ! Non-constant (space dependent) fields
              !  
              if (kfl_infun_nsa < 4) then
                 !
                 ! Stratified atmospheres, typical of meteo problems
                 !    
                 if (kfl_inico_nsa==0) then
                    do ipoin = 1,npoin                       
                       veloc(1:ndime,ipoin,icomp) = 0.0_rp
                    end do
                 else if (kfl_inico_nsa==1) then
                    do idime=1,ndime
                       do ipoin = 1,npoin                       
                          veloc(idime,ipoin,icomp) = veloc_nsa(idime)        
                       end do
                    end do
                 end if
                 !
                 ! Meteo initial fields are either computed or read:
                 !
                 call nsa_inimet(1_ip)

              else
                 call runend('NSA_INIUNK: INITIAL FIELDS FUNCTION NOT PROGRAMMED')
              end if

           else if (kfl_infun_nsa == 0) then

              !
              ! Constant (space independent) values given in the input file
              !                           

              if (kfl_inico_nsa==0) then
                 do ipoin = 1,npoin                       
                    veloc(1:ndime,ipoin,icomp) = 0.0_rp
                 end do
              else if (kfl_inico_nsa==1) then
                 do idime=1,ndime
                    do ipoin = 1,npoin                       
                       veloc(idime,ipoin,icomp) = veloc_nsa(idime)        
                    end do
                 end do
              end if

              !
              ! Properties from the kernel: mu, molecular weight and Cp
              !
              if (kfl_prope /= 0 ) then
                 !call ker_proper('VISCO','NPOIN',dummi,dummi,visco(:,icomp),dummi,dummi,xrano,xrano)
                 !call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp_nsa,dummi,dummi,xrano,xrano)
                 call ker_proper('VISCO','NPOIN',dummi,dummi,visco(:,icomp))
                 call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp_nsa)
                 if (kfl_coupl(ID_NASTAL,ID_CHEMIC) == 0 ) then
                    wmean = mowei_nsa
                 endif
              else
                 wmean     = mowei_nsa
                 shecp_nsa = cpcoe_nsa
              endif

              if (kfl_inibu_nsa == 0) then
                 !
                 ! The normal no-initial-fields run
                 !
                 do ipoin = 1,npoin                       

                    xhecv = shecp_nsa(ipoin) - runiv_nsa / wmean(ipoin,1)

                    tempe(ipoin,icomp) = tempe_nsa
                    press(ipoin,icomp) = press_nsa
                    densi(ipoin,icomp) = densi_nsa


                    if (kfl_iposi_nsa == 1) call nsa_iniunk_iposi(ipoin,icomp)

                    !
                    ! Properties: viscosity mu  
                    !
                    if (kfl_prope == 0 ) then
                       visco(ipoin,icomp) = visco_nsa
                    end if

                    velmi = 0.0_rp
                    do idime=1,ndime
                       velmi = velmi + veloc(idime,ipoin,icomp) * veloc(idime,ipoin,icomp)
                       umome(idime,ipoin,icomp) = densi_nsa * veloc(idime,ipoin,icomp)
                    end do
                    energ(ipoin,icomp) = densi_nsa * (xhecv * tempe_nsa + 0.5_rp*velmi)

                 end do
              else if (kfl_inibu_nsa == 1) then
                 !
                 ! laxliu initial built-in field 
                 ! 
                 ! squared centered at (llcen_nsa(1),llcen_nsa(2)), side=1
                 !
                 
                 do ipoin = 1,npoin                       

                    if (coord(1,ipoin) .le. llcen_nsa(1)) then
                       if (coord(2,ipoin) .le. llcen_nsa(2)) then
                          !
                          ! Quadrant 1 - LB
                          !
                          densi(ipoin,icomp)     = llrho_nsa(    1)
                          press(ipoin,icomp)     = llpre_nsa(    1)
                          veloc(1:2,ipoin,icomp) = llvel_nsa(1:2,1)

!                          kfl_fixno_nsa(2,ipoin) = 0
!                          kfl_fixno_nsa(3,ipoin) = 0
!                          kfl_fixno_nsa(4,ipoin) = 0


                       else
                          !
                          ! Quadrant 4 - LT
                          !
                          densi(ipoin,icomp)     = llrho_nsa(    4)
                          press(ipoin,icomp)     = llpre_nsa(    4)
                          veloc(1:2,ipoin,icomp) = llvel_nsa(1:2,4)

!                          kfl_fixno_nsa(1,ipoin) = 0
!                          kfl_fixno_nsa(3,ipoin) = 0
!                          kfl_fixno_nsa(4,ipoin) = 0

                       end if
                    else 
                       if (coord(2,ipoin) .le. llcen_nsa(2)) then
                          !
                          ! Quadrant 2 - RB
                          !
                          densi(ipoin,icomp)     = llrho_nsa(    2)
                          press(ipoin,icomp)     = llpre_nsa(    2)
                          veloc(1:2,ipoin,icomp) = llvel_nsa(1:2,2)
                       else
                          !
                          ! Quadrant 3 - RT
                          !
                          densi(ipoin,icomp)     = llrho_nsa(    3)
                          press(ipoin,icomp)     = llpre_nsa(    3)
                          veloc(1:2,ipoin,icomp) = llvel_nsa(1:2,3)
                       end if
                    end if

                    xhecv = shecp_nsa(ipoin) - runiv_nsa / wmean(ipoin,1)
                    rgasc = runiv_nsa / wmean(ipoin,1)
                    tempe(ipoin,icomp) = press(ipoin,icomp) / densi(ipoin,icomp) / rgasc

                    ! 
                    ! Laxliu problem fixes bvess to the initial values
                    !
                    
                    bvess_nsa(1:ndime,ipoin,1)= veloc(1:ndime,ipoin,icomp)
                    bvess_nsa(1+ndime,ipoin,1)= densi(        ipoin,icomp)
                    bvess_nsa(2+ndime,ipoin,1)= tempe(        ipoin,icomp)


                    !
                    ! Properties: viscosity mu  
                    !
                    if (kfl_prope == 0 ) then
                       visco(ipoin,icomp) = visco_nsa
                    end if

                    velmi = 0.0_rp
                    do idime=1,ndime
                       velmi = velmi + veloc(idime,ipoin,icomp) * veloc(idime,ipoin,icomp)
                       umome(idime,ipoin,icomp) = densi(ipoin,icomp) * veloc(idime,ipoin,icomp)
                    end do
                    energ(ipoin,icomp) = densi(ipoin,icomp) * (xhecv * tempe(ipoin,icomp) + 0.5_rp*velmi)

                 end do

              end if


              !
              ! Correct boundary conditions and recompute variables 
              !        
              call nsa_setvar(zero,icomp)

           end if

        else 

           !
           ! Initialize variables from given fields
           !        
           if (kfl_inico_nsa==0) then
              do ipoin = 1,npoin                       
                 veloc(1:ndime,ipoin,icomp) = 0.0_rp
              end do
           else if (kfl_inico_nsa==1) then
              do idime=1,ndime
                 do ipoin = 1,npoin                       
                    veloc(idime,ipoin,icomp) = veloc_nsa(idime)        
                 end do
              end do
           end if

           !
           ! bvess_nsa = initial_fields 
           ! (ONLY for those nodes with no Dirichlet boundary condition imposed)
           ! !!!!! ESTO ESTABA MAL, PORQUE SI HABIA BC AL FINAL NO FIJABA NADA
           ! !!!!! HE COMENTADO LO DE MIRAR FIXNO
           !
           if (kfl_inifi_nsa(1) > 0) then ! mom. eq.: velocity given
              do ipoin = 1,npoin 
                 do idime=1,ndime
                    !                    if(kfl_fixno_nsa(idime,ipoin)<= 0 .or. kfl_fixno_nsa(idime,ipoin) > 2) then
                    bvess_nsa(idime,ipoin,1)= xfiel(-nfiel_nsa(1))%a(idime,ipoin,1)
                    !                    end if
                 end do
              end do
           end if

           if (kfl_inifi_nsa(2) > 0) then ! cont. eq. : density or pressure given 
              do ipoin = 1,npoin       
                 !                 if(kfl_fixno_nsa(ndime+1,ipoin)<= 0) then
                 bvess_nsa(ndime+1,ipoin,1)= xfiel(-nfiel_nsa(2))%a(1,ipoin,1)
                 !                 end if
              end do
           end if

           if (kfl_inifi_nsa(3) > 0) then ! ener. eq. : temperature or energy given
              do ipoin = 1,npoin   
                 !                 if(kfl_fixno_nsa(ndime+2,ipoin)<= 0) then
                 bvess_nsa(ndime+2,ipoin,1)= xfiel(-nfiel_nsa(3))%a(1,ipoin,1)
                 !                 end if
              end do
           end if

           !
           ! Properties from the kernel: molecular weight and Cp
           !
           if (kfl_prope /= 0 ) then
              !call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp_nsa,dummi,dummi,xrano,xrano)
              call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp_nsa)
              if (kfl_coupl(ID_NASTAL,ID_CHEMIC) == 0 ) then
                 wmean = mowei_nsa
              endif
           else
              wmean     = mowei_nsa
              shecp_nsa = cpcoe_nsa
           endif

           do ipoin = 1,npoin

              xhecv = shecp_nsa(ipoin) - runiv_nsa / wmean(ipoin,1)

              if (kfl_brunt_nsa == 2) then
                 stapa_nsa(1) = real(brunt_nsa(ipoin),rp)
              else
                 stapa_nsa(1) = brure_nsa
              end if
              if (kfl_inifi_nsa(1) == 1) then
                 do idime=1,ndime
                    veloc(idime,ipoin,icomp)= bvess_nsa(idime,ipoin,1)
                 end do
              end if
              velmi = 0.0_rp
              do idime=1,ndime
                 velmi = velmi + &
                      veloc(idime,ipoin,icomp) * veloc(idime,ipoin,icomp) 
              end do
              velmi = sqrt(velmi)

              if (kfl_inifi_nsa(3) == 1) then
                 tempe(ipoin,icomp) = bvess_nsa(ndime+2,ipoin,1)
                 if (kfl_inifi_nsa(2) .le. 1) then
                    densi(ipoin,icomp) = bvess_nsa(ndime+1,ipoin,1)              

                    call nsa_stalaw(2_ip,1_ip,densi(ipoin,icomp),press(ipoin,icomp), &
                         tempe(ipoin,icomp),coord(ivert_nsa,ipoin),dummr,wmean(ipoin,1),shecp_nsa(ipoin))        ! compute p from rho and T given
                 else if (kfl_inifi_nsa(2) == 2) then                            
                    press(ipoin,icomp) = bvess_nsa(ndime+1,ipoin,1)

                    if (kfl_iposi_nsa == 1) call nsa_iniunk_iposi(ipoin,icomp)                                   ! SO FAR, IT CAN BE DONE FOR P

                    call nsa_stalaw(1_ip,1_ip,densi(ipoin,icomp),press(ipoin,icomp), &
                         tempe(ipoin,icomp),coord(ivert_nsa,ipoin),dummr,wmean(ipoin,1),shecp_nsa(ipoin))        ! compute rho from p and T given
                 end if

                 energ(ipoin,icomp) =  densi(ipoin,icomp) &
                      * (xhecv * tempe(ipoin,icomp) + 0.5_rp*velmi*velmi)

              else if (kfl_inifi_nsa(3) == 2) then
                 energ(ipoin,icomp) = bvess_nsa(ndime+2,ipoin,1)
                 if (kfl_inifi_nsa(2) .le. 1) then
                    densi(ipoin,icomp) = bvess_nsa(ndime+1,ipoin,1)              
                    tempe(ipoin,icomp) = (energ(ipoin,icomp)/densi(ipoin,icomp) - 0.5_rp*velmi*velmi) / xhecv

                    call nsa_stalaw(2_ip,1_ip,densi(ipoin,icomp),press(ipoin,icomp),&
                         tempe(ipoin,icomp),coord(ivert_nsa,ipoin),dummr,wmean(ipoin,1),shecp_nsa(ipoin))        ! compute p from rho and T given
                 else if (kfl_inifi_nsa(2) == 2) then              

                    rgasc = runiv_nsa / wmean(ipoin,1)               

                    press(ipoin,icomp) = bvess_nsa(ndime+1,ipoin,1)
                    tempe(ipoin,icomp) = press(ipoin,icomp) / densi(ipoin,icomp) / rgasc

                    call nsa_stalaw(1_ip,1_ip,densi(ipoin,icomp),press(ipoin,icomp),&
                         tempe(ipoin,icomp),coord(ivert_nsa,ipoin),dummr,wmean(ipoin,1),shecp_nsa(ipoin))        ! compute rho from p and T given
                 end if

              end if
              !              if (kfl_inifi_nsa(2) == 2 .or. kfl_relat_nsa==1) then
              !                 press(ipoin,icomp) = bvess_nsa(ndime+2,ipoin,1)
              !                 call nsa_stalaw(3,densi(ipoin,icomp),press(ipoin,icomp),&
              !                      tempe(ipoin,icomp),coord(ivert_nsa,ipoin),dummr)
              !              end if
              !              if (kfl_inifi_nsa(3) == 1) then
              !                 tempe(ipoin,icomp) = bvess_nsa(ndime+2,ipoin,1) 
              !                 call nsa_stalaw(2,densi(ipoin,icomp),press(ipoin,icomp),&
              !                      tempe(ipoin,icomp),coord(ivert_nsa,ipoin),dummr)
              !              end if

              umome(1:ndime,ipoin,icomp) =  densi(ipoin,icomp) &
                   * veloc(1:ndime,ipoin,icomp)

           end do

        end if

        if (kfl_hysta_nsa > 0) then
           if (kfl_inkee_nsa(2) > 0 .and. kfl_inkee_nsa(3) > 0 .and. kfl_inkee_nsa(4) > 0) then
              do ipoin = 1,npoin                                     
                 rekee_nsa(ndime+1, ipoin) = xfiel(-nfiel_nsa(4))%a(1,ipoin,1)  ! hydrostatic density
                 rekee_nsa(ndime+2, ipoin) = xfiel(-nfiel_nsa(5))%a(1,ipoin,1)  ! hydrostatic temperature
                 rekee_nsa(ndime+3, ipoin) = xfiel(-nfiel_nsa(6))%a(1,ipoin,1)  ! hydrostatic pressure
              end do
           else
              !           call aaa() ! Compute hydrostatic density, hydrostatic temperature and hydrostatic pressure on nodes and store in rekee
           end if
        else if (kfl_inkee_nsa(2) > 0 .or. kfl_inkee_nsa(3) > 0 .or. kfl_inkee_nsa(4) > 0) then
           if (kfl_inkee_nsa(2) > 0) then
              do ipoin = 1,npoin
                 rekee_nsa(ndime+1, ipoin) = xfiel(-nfiel_nsa(4))%a(1,ipoin,1)  ! hydrostatic density
              end do
           end if
           if (kfl_inkee_nsa(3) > 0) then
              do ipoin = 1,npoin
                 rekee_nsa(ndime+2, ipoin) = xfiel(-nfiel_nsa(5))%a(1,ipoin,1)  ! hydrostatic temperature
              end do
           end if
           if (kfl_inkee_nsa(4) > 0) then
              do ipoin = 1,npoin
                 rekee_nsa(ndime+3, ipoin) = xfiel(-nfiel_nsa(6))%a(1,ipoin,1)  ! hydrostatic pressure
              end do
           end if
        end if


        do itime=2+kfl_tiacc_nsa,4,-1
           do ipoin = 1,npoin
              veloc(1:ndime,ipoin,itime) = veloc(1:ndime,ipoin,itime-1) 
              umome(1:ndime,ipoin,itime) = umome(1:ndime,ipoin,itime-1) 
              press(        ipoin,itime) = press(        ipoin,itime-1)  
              densi(ipoin,itime) = densi(ipoin,itime-1) 
              tempe(ipoin,itime) = tempe(ipoin,itime-1) 
              energ(ipoin,itime) = energ(ipoin,itime-1) 
           end do
        end do

     end if

  end if

  !  
  ! check initial minmax for the velocity module, use vmach_nsa vector and store values in vemxm_nsa
  !
  if( INOTMASTER ) then
     do ipoin = 1,npoin
        vmach_nsa(ipoin)=0.0_rp
        do idime= 1,ndime
           vmach_nsa(ipoin)= vmach_nsa(ipoin) + veloc(idime,ipoin,1)*veloc(idime,ipoin,1)
        end do
        vmach_nsa(ipoin)=sqrt(vmach_nsa(ipoin))
     end do
  end if

  kfl_zevel_nsa = 0
  !  call minmax(one,npoin,zero,vmach_nsa,vemxm_nsa(1),vemxm_nsa(2))
  !  if (vemxm_nsa(2) == 0.0_rp) then
  !     kfl_zevel_nsa(1) = 1
  !     kfl_zevel_nsa(2) = kfl_lopre_nsa
  !     kfl_lopre_nsa    = 0
  !  end if

  call nsa_updunk(10_ip)      ! Initialize all the fields for state variables (:,:,1),(:,:,2),(:,:,3)
  !
  ! Output initial solution
  !
  !!call nsa_output(-1_ip)

  !if (kfl_pro2d_nsa == 1) then
  !
  ! Initialize profile postprocessing
  !
  ! call nsa_outpro(zero)
  ! call nsa_outpro(one)

  ! end if

  call nsa_coupli(ITASK_INIUNK)

end subroutine nsa_iniunk

!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_iniunk_iposi.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Set up the initial condition for the positional fields
!> @details Set up the initial condition for the positional fields
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_iniunk_iposi(ipoin,icomp)
  use      def_master
  use      def_domain
  use      def_parame
  use      def_nastal

  ! so far programmed only for PRESSURE in a sphere

  integer(ip)  :: ipoin
  integer(ip)  :: icomp
  real(rp)     :: xradi,x(3),rgasc

  if (iposi_nsa(2)%kflag == 1) then
     if (iposi_nsa(2)%geometry == 1) then
        x = 0.0_rp
        x(1) = (coord(1,ipoin) - iposi_nsa(2)%center(1)) 
        x(2) = (coord(2,ipoin) - iposi_nsa(2)%center(2)) 
        x(3) = (coord(3,ipoin) - iposi_nsa(2)%center(3)) 

        xradi= x(1)*x(1)+x(2)*x(2)
        if (ndime == 3) xradi = xradi + x(3)*x(3)

        xradi = sqrt(xradi)

        if (xradi <= iposi_nsa(2)%radius) then
           rgasc = runiv_nsa / wmean(ipoin,1)
           press(ipoin,icomp)=  iposi_nsa(2)%value 
           densi(ipoin,icomp)=  press(ipoin,icomp) / rgasc / tempe(ipoin,icomp) 
        end if
     end if
  end if
 

end subroutine nsa_iniunk_iposi

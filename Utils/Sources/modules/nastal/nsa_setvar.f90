!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_setvar.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Compute derived fields
!> @details Compute derived fields
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_setvar(ktask,icomp)
  use      def_master
  use      def_domain
  use      def_parame
  use      def_kermod
  use      mod_ker_proper
  use      def_nastal

  implicit none
  integer(ip)       :: ktask,icomp,ipoin,kpoin,idime,kfixi,idofn,imode,dummi
  integer(ip)       :: ievat,ifstg,kinfl,kfl_fixau(ndofn_nsa)
  real(rp)          :: vmodu,vauxi,rdumy,sound,xmach,&
       soden,pauxi,oltem,olpre,olden,olvel(ndime),xcons(ndofn_nsa),xphys(ndofn_nsa),&
       xadve(ndime),xvmsh(ndime),facol,facne,dummy(ndime,ndime),xrano(3)
  real(rp)          :: xconc !SM mixing ratio of water vapor only
  real (rp)         :: auxvi(ndofn_nsa)
  real(rp)          :: xhecv,adgam,ubulk
  real(rp)          :: xdtot

  !Initialize xconc:
  xconc = 0.0_rp
  imode = 0

  !Initialization of dummy variables for viscosity
  auxvi  = 0.0_rp
  rdumy = 0.0_rp

  if (ktask == zero .and. INOTMASTER) then                                            ! called by iniunk

     !
     ! SETTING VARIABLES FROM INIUNK: FROM REFERENCE VALUES
     !
     ! Properties from the kernel: molecular weight,Cp & visco
     !
     if (kfl_prope /= 0 ) then
        !call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp_nsa,dummi,dummi,xrano,xrano)
        !call ker_proper('VISCO','NPOIN',dummi,dummi,visco(:,icomp),dummi,dummi,xrano,xrano)
        call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp_nsa)
        call ker_proper('VISCO','NPOIN',dummi,dummi,visco(:,icomp))
        if (kfl_coupl(ID_NASTAL,ID_CHEMIC) == 0 ) then
           wmean = mowei_nsa
        endif
        imode = 1 
     else
        wmean     = mowei_nsa
        shecp_nsa = cpcoe_nsa
     endif

     do ipoin = 1,npoin

        xhecv = shecp_nsa(ipoin) - runiv_nsa / wmean(ipoin,1)

        ! Local working copy
        do idime=1,ndime
           xphys(idime) = veloc(idime,ipoin,icomp)           
           xadve(idime) = xphys(idime)
           xvmsh(idime) = 0.0_rp
        end do
        ! When coupled with alefor, substract the mesh velocity velom to the advection velocity xadve
        if( kfl_coupl(ID_NASTAL,ID_ALEFOR) /= 0 ) then  
           do idime=1,ndime
              xvmsh(idime) = velom(idime,ipoin)
              xadve(idime) = xadve(idime) - xvmsh(idime)
           end do
        end if

        xphys(ndime+1) = press(        ipoin,icomp)

        do idime=1,ndime
           xcons(idime) = umome(idime,ipoin,icomp)
        end do
        xcons(ndime+1) = densi(        ipoin,icomp)
        xcons(ndime+2) = energ(        ipoin,icomp)
        xphys(ndime+2) = tempe(        ipoin,icomp)

        kfl_fixau(1:ndime)=kfl_fixno_nsa(1:ndime,ipoin)
        kfl_fixau(ndime+1)=kfl_fixno_nsa(ndime+1,ipoin)
        kfl_fixau(ndime+2)=kfl_fixno_nsa(ndime+2,ipoin)

        if(kfl_fixau(1)==5)then
           !
           ! First, check inflow-outflow prescription
           ! Check whether it is inflow or outflow, sub or supersonic
           !                 

           call nsa_chkinf(&
                kinfl,xmach,ipoin,xadve(1:ndime),xphys(1:ndime),xphys(ndime+1),xcons(ndime+1))

           !
           ! Uncomment these lines if fixno for condition "5" should be set once at the beginning.
           ! It is preferred not to do it, but just in case... 
           !


!!$           if (kinfl==1) then
!!$              ! inflow
!!$              kfl_fixno_nsa(1:ndime,ipoin)=1
!!$              kfl_fixno_nsa(ndime+1,ipoin)=0
!!$              kfl_fixno_nsa(ndime+2,ipoin)=1
!!$              if (xmach >= 1.0_rp) then 
!!$                 ! supersonic, so fix also the density
!!$                 kfl_fixno_nsa(ndime+1,ipoin)=1
!!$              end if
!!$           else
!!$              ! outflow
!!$              kfl_fixno_nsa(1:ndime,ipoin)=0
!!$              kfl_fixno_nsa(ndime+2,ipoin)=0
!!$              if (xmach < 1.0_rp) then 
!!$                 ! subsonic, so fix only the density
!!$                 kfl_fixno_nsa(ndime+1,ipoin)=1
!!$              end if
!!$           end if
!!$
           if (kinfl==1) then
              ! inflow
              kfl_fixau(1:ndime)=1
              kfl_fixau(ndime+1)=0
              kfl_fixau(ndime+2)=1
              if (xmach >= 1.0_rp) then 
                 ! supersonic, so fix also the density
                 kfl_fixau(ndime+1)=1
              end if
           else
              ! outflow
              kfl_fixau(1:ndime)=0
              kfl_fixau(ndime+2)=0
              if (xmach < 1.0_rp) then 
                 ! subsonic, so fix only the density
                 kfl_fixau(ndime+1)=1
              end if
           end if



        end if


        !
        ! Velocity prescriptions
        !
        ifstg= 0

        if(kfl_fixau(1)==2)then
           call nsa_rotunk( one,ipoin,ipoin,xadve(1:ndime))
           do idime=1,ndime
              if(kfl_fixau(idime) == 2)  then
                 xadve(idime) = 0.0_rp                       ! normal to zero
              else if(kfl_fixau(idime) == 1)  then
                 xadve(idime) = bvess_nsa(idime,ipoin,1)     ! combined normal with fixity... weird case
              end if
           end do
           call nsa_rotunk(mone,ipoin,ipoin,xadve(1:ndime))

           do idime= 1,ndime
              xphys(idime) = xadve(idime) + xvmsh(idime)
           end do

        else 

           do idime=1,ndime
              if(kfl_fixau(idime)==1)then
                 xadve(idime) = bvess_nsa(idime,ipoin,1)
                 xphys(idime) = xadve(idime) + xvmsh(idime)
                 ifstg=ifstg+1
              end if
           end do

        end if

        if (kfl_fixau(ndime+1) == 1 ) then
           !
           ! Continuity equation prescription: density
           xcons(ndime+1) = bvess_nsa(ndime+1,ipoin,1)
        else if(kfl_fixau(ndime+1) == 2 ) then
           !
           ! Stop the code, this is a deprecated possibility
           call runend ('NSA_SETVAR: FIX PRESSURE THROUGH THE ENERGY EQUATION, NOT CONTINUITY.')
        end if


        !
        ! Energy equation prescription: energy, temperature or pressure
        if (kfl_fixau(ndime+2) > 0) then
           if (kfl_fixau(ndime+2) == 1) xphys(ndime+2) = bvess_nsa(ndime+2,ipoin,1)   ! temperature prescribed
           if (kfl_fixau(ndime+2) == 2) then
              xphys(ndime+1) = bvess_nsa(ndime+2,ipoin,1)                             ! pressure prescribed
              call nsa_stalaw(3,imode,xcons(ndime+1),xphys(ndime+1),xphys(ndime+2),&  
                   rdumy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))                       ! temperature computed
           end if
           if (kfl_fixau(ndime+2) == 3) xphys(ndime+2) = tstag_nsa          ! stagnation temp prescription
        end if

        !
        ! Compute velocity module
        !
        vmodu = xphys(1) * xphys(1) + xphys(2) * xphys(2)
        if (ndime == 3) vmodu = vmodu + xphys(3) * xphys(3)


        !
        ! Prescribe Total Pressure and compute static pressure (and fix it)           
        ! Done in the continuity equation
        !
        if(kfl_fixau(ndime+1) == 3 ) then
           adgam = shecp_nsa(ipoin) / xhecv  ! gamma
!!!    la verdadera uinlet se calcula en una subru ad-hoc
!!!           uinlet_nsa = sqrt(vmodu)  ! OJO: esto es para probar, pues debería ser una media de la entrada y no la local
           vauxi= (1.0_rp + (adgam - 1.0_rp) * uinlet_nsa * uinlet_nsa / (2.0_rp * adgam * rgasc_nsa * xphys(ndime+2)))
           vauxi= vauxi ** (adgam / (adgam - 1.0_rp))
           vauxi= bvess_nsa(ndime+1,ipoin,1) / vauxi                     ! pressure computed
           xphys(ndime+1) = vauxi                     ! pressure computed
           call nsa_stalaw(1,imode,xcons(ndime+1),xphys(ndime+1),xphys(ndime+2),&  
                rdumy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))                       ! density computed
        end if

        if (kfl_brunt_nsa == 2) then
           stapa_nsa(1) = brunt_nsa(ipoin)
        else
           stapa_nsa(1) = brure_nsa
        end if

        if (kfl_fixau(ndime+1) >= 0 ) then

           call nsa_stalaw(2,imode,xcons(ndime+1),xphys(ndime+1),xphys(ndime+2),&
                rdumy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))

        end if

        if ( kfl_isent_nsa == 1) then                  ! isentropic flows
           call nsa_stalaw(3,imode,xcons(ndime+1),xphys(ndime+1),xphys(ndime+2),&
                rdumy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))
        end if

        do idime=1,ndime
           xcons(idime) = xcons(ndime+1) * xphys(idime)
        end do
        xcons(ndime+2) = xcons(ndime+1) *(xhecv *  xphys(ndime+2) + 0.5_rp * vmodu)

        if (kfl_prope == 0 ) then
           call nsa_lawvis(ipoin,icomp,rdumy,rdumy,rdumy) ! Recalculation of viscosity
        end if

        press(ipoin,icomp) = xphys(ndime+1) 
        tempe(ipoin,icomp) = xphys(ndime+2)
        energ(ipoin,icomp) = xcons(ndime+2)
        densi(ipoin,icomp) = xcons(ndime+1)

        do idime=1,ndime
           umome(idime,ipoin,icomp) = xcons(idime)
           veloc(idime,ipoin,icomp) = xphys(idime) 
        end do

     end do

  else if (ktask == three) then   ! called by upcons (CURRENT CONS. UNKNOWNS ARE STORED IN UNKNO!!)

     if (kfl_inlet_nsa /= 0_ip) then
        call nsa_inflow(kfl_inlet_nsa,ubulk)
     end if

     if (INOTMASTER) then

        !
        ! Properties from the kernel: molecular weight,Cp & visco
        !
        if (kfl_prope /= 0 ) then
           !call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp_nsa,dummi,dummi,xrano,xrano)
           !call ker_proper('VISCO','NPOIN',dummi,dummi,visco(:,icomp),dummi,dummi,xrano,xrano)
           call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp_nsa)
           call ker_proper('VISCO','NPOIN',dummi,dummi,visco(:,icomp))
           if (kfl_coupl(ID_NASTAL,ID_CHEMIC) == 0 ) then
              wmean = mowei_nsa
           endif
           imode = 1 
        else
           wmean     = mowei_nsa
           shecp_nsa = cpcoe_nsa
        endif

        !
        ! Meteo with initial conditions computed enters here.
        ! called by upcons (CURRENT CONS. UNKNOWNS ARE STORED IN UNKNO!!)
        !
        do ipoin = 1,npoin

           ! Local work copy
           do idime=1,ndime
              xphys(idime) = veloc(idime,ipoin,icomp) !u,v
           end do

           kfl_fixau(1:ndime)=kfl_fixno_nsa(1:ndime,ipoin)
           kfl_fixau(ndime+1)=kfl_fixno_nsa(ndime+1,ipoin)
           kfl_fixau(ndime+2)=kfl_fixno_nsa(ndime+2,ipoin)

           !
           ! Store xphys for pressure (total pressure regardless of the set of variables that we are solving for):
           !
           xphys(ndime+1) = press(ipoin,icomp)

           ievat = (ipoin-1)*ndofn_nsa + ndime
           if (kfl_foreg_nsa > 0) call runend('NSA_SETVAR: NO FORCED INCOMPRESSIBLE FLOW CAN BE SOLVED USING NASTAL')

           if (kfl_unkse_nsa < 10) then
              xphys(ndime+2) = tempe(ipoin,icomp)
              xcons(ndime+2) = unkno(ievat+2)
           else if (kfl_unkse_nsa < 20) then
              !
              ! In meteo we are here (kfl_unkse_nsa=10)
              !
              ! At this point unkno(:) of the new time step (current)
              ! has already been computed by the solver as:
              !
              !     unkno(:) = unkno(:) + dt*rhs(:)/mass(:)
              !
              !Store the physical solution variables from unkno:
              xphys(ndime+2) = unkno(ievat+2)     !Tempe   or tempe'
           end if
           xcons(ndime+1) = unkno(ievat+1)        !Density or densi'

           if(kfl_ncons_nsa == 0) then
              !
              ! Momentum in conservative formulation:
              !
              do idime= 1,ndime
                 ievat = (ipoin-1)*ndofn_nsa + idime
                 xcons(idime) = unkno(ievat)
              end do
           else
              !
              ! Momentum in NON-conservative formulation:
              !
              xdtot = xcons(ndime+1)
              do idime= 1,ndime
                 ievat = (ipoin-1)*ndofn_nsa + idime
                 xcons(idime) = unkno(ievat)*xdtot !U and V momentum (dime=1, idime=ndime)
              end do
           end if

           vmodu = xphys(1) * xphys(1) + xphys(2) * xphys(2)
           if (ndime == 3) vmodu = vmodu + xphys(3) * xphys(3)

           if (kfl_isent_nsa == 1) then
              xhecv     = shecp_nsa(ipoin)   - runiv_nsa / wmean(ipoin,1)
              call nsa_stalaw(3,imode,xcons(ndime+1),xphys(ndime+1),xphys(ndime+2),rdumy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))                 
              xcons(ndime+2) = xcons(ndime+1) * (xhecv * xphys(ndime+2) + 0.5_rp * vmodu)           
           end if

           ! 1. Compute physical variables from conservative ones

           oltem         = tempe(ipoin,1)
           olden         = densi(ipoin,1)
           olpre         = press(ipoin,1)
           olvel(1:ndime)= veloc(1:ndime,ipoin,1)

           if(kfl_fixau(1)==9) then
              call nsa_rotunk( one,ipoin,ipoin,olvel(1:ndime))
           end if
           if (xcons(ndime+1) < zensa) then
              xcons(ndime+1) = olden
           end if

           do idime=1,ndime
              xphys(idime) = xcons(idime) / xcons(ndime+1)
              xadve(idime) = xphys(idime)
              xvmsh(idime) = 0.0_rp
           end do

           vmodu = xphys(1) * xphys(1) + xphys(2) * xphys(2)
           if (ndime == 3) vmodu = vmodu + xphys(3) * xphys(3)
           xhecv = shecp_nsa(ipoin) - runiv_nsa / wmean(ipoin,1)

           ! When coupled with alefor, substract the mesh velocity velom to the advection velocity xadve
           if( kfl_coupl(ID_NASTAL,ID_ALEFOR) /= 0 ) then  
              do idime=1,ndime
                 xvmsh(idime) = velom(idime,ipoin)
                 xadve(idime) = xadve(idime) - xvmsh(idime)
              end do
           end if

           if (kfl_unkse_nsa < 10) then
              if ( kfl_isent_nsa == 1) then                  ! isentropic flows
                 call nsa_stalaw(3,imode,xcons(ndime+1),xphys(ndime+1),xphys(ndime+2),rdumy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))                 
              else
                 xhecv            = shecp_nsa(ipoin)   - runiv_nsa / wmean(ipoin,1)
                 xphys(ndime+2)   = (xcons(ndime+2)   /  xcons(ndime+1)  - 0.5_rp * vmodu) / xhecv
              end if
           else if (kfl_unkse_nsa < 20) then
              xhecv          = shecp_nsa(ipoin)   - runiv_nsa / wmean(ipoin,1)
              xcons(ndime+2) = xcons(ndime+1) *(xhecv * xphys(ndime+2) + 0.5_rp * vmodu)
           end if

           if (xphys(ndime+2) < zensa) then
              xphys(ndime+2) = oltem
              xhecv          = shecp_nsa(ipoin)   - runiv_nsa / wmean(ipoin,1)
              xcons(ndime+2) = xcons(ndime+1) *(xhecv * xphys(ndime+2) + 0.5_rp * vmodu)
           end if

           call nsa_stalaw(2,imode,xcons(ndime+1),xphys(ndime+1),xphys(ndime+2),rdumy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))

           !
           ! Impose boundary conditions
           !
           kfixi = 0
           do idofn= 1,ndofn_nsa
              kfixi = kfixi + kfl_fixau(idofn)
           end do
           if (kfixi > 0) then                                                  ! Prescriptions ?

              ! 2. Correct physical variables according to bc imposed

              if(kfl_fixau(1)==5)then
                 !
                 ! First, check inflow-outflow prescription
                 ! Check whether it is inflow or outflow, sub or supersonic
                 !                 
                 call nsa_chkinf(&
                      kinfl,xmach,ipoin,xadve(1:ndime),xphys(1:ndime),xphys(ndime+1),xcons(ndime+1))

                 if (kinfl==1) then
                    ! inflow
                    kfl_fixau(1:ndime)=1
                    kfl_fixau(ndime+2)=1
                    if (xmach >= 1.0_rp) then 
                       ! supersonic, so fix also the density
                       kfl_fixau(ndime+1)=1
                    end if
                 else
                    ! outflow
                    if (xmach < 1.0_rp) then 
                       ! subsonic, so fix only the density
                       kfl_fixau(ndime+1)=1
                    end if
                 end if
              end if

              ! Non-reflecting (or open boundaries) boundary conditions
              if(kfl_fixau(1)==9) then
                 facol= 0.5_rp
                 facne= 0.5_rp
                 ! rotate to local normal-tangents basis                 
                 do idime=1,ndime
                    xphys(idime)= -xphys(idime)
                    olvel(idime)= -olvel(idime)
                 end do
                 call nsa_rotunk( one,ipoin,ipoin,xphys(1:ndime))

                 adgam = shecp_nsa(ipoin) / xhecv  ! gamma

                 sound= sqrt(adgam * (facol*olpre+facne*xphys(ndime+1)) / (facol*olden+facne*xcons(ndime+1)))
                 soden= facol*         olden*sqrt(adgam * olpre          / olden         )&
                      + facne*xcons(ndime+1)*sqrt(adgam * xphys(ndime+1) / xcons(ndime+1))
                 ! pressure (auxiliar value)
                 pauxi          = 0.5_rp * (olpre + xphys(ndime+1) + soden*(olvel(1) - xphys(1)))
                 ! density
                 xcons(ndime+1) = xcons(ndime+1) + (pauxi - xphys(ndime+1))/sound/sound
                 ! normal velocity
                 xphys(      1) = olvel(1) + (olpre - pauxi)/soden
                 ! pressure 
                 xphys(ndime+1) = pauxi
                 ! rotate back to cartesian basis
                 call nsa_rotunk(mone,ipoin,ipoin,xphys(1:ndime))
                 do idime=1,ndime
                    xphys(idime)= -xphys(idime)
                 end do
                 ! 
                 ! compute the rest of the variables

                 !   temperature
                 !
                 call nsa_stalaw(3,imode,xcons(ndime+1),xphys(ndime+1),xphys(ndime+2),rdumy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))

                 !   energy and momentum
                 !
                 xcons(1) = xphys(1) * xcons(ndime+1)
                 xcons(2) = xphys(2) * xcons(ndime+1)
                 vmodu = xphys(1) * xphys(1) + xphys(2) * xphys(2)
                 if (ndime == 3) then
                    vmodu = vmodu + xphys(ndime) * xphys(ndime)
                    xcons(ndime) = xphys(ndime) * xcons(ndime+1)
                 end if
                 if (kfl_unkse_nsa < 10) then
                    xhecv        = shecp_nsa(ipoin)   - runiv_nsa / wmean(ipoin,1)
                    xphys(ndime+2) = (xcons(ndime+2) / xcons(ndime+1) - 0.5_rp * vmodu) / xhecv
                 else if (kfl_unkse_nsa < 20) then
                    xhecv        = shecp_nsa(ipoin)   - runiv_nsa / wmean(ipoin,1)
                    xcons(ndime+2) = xcons(ndime+1) *(xhecv * xphys(ndime+2) + 0.5_rp * vmodu)
                 end if

                 !
                 ! Continuity equation prescription:                   
                 !
              else if(kfl_fixau(ndime+1) == 1 ) then                    ! density
                 xcons(ndime+1) = bvess_nsa(ndime+1,ipoin,1)                 
                 oltem= tempe(ipoin,1)
                 do idime=1,ndime
                    xphys(idime) = xcons(idime) / xcons(ndime+1)
                 end do
                 vmodu = xphys(1) * xphys(1) + xphys(2) * xphys(2)
                 if (ndime == 3) vmodu = vmodu + xphys(3) * xphys(3)

                 if ( kfl_isent_nsa == 1) then                  ! isentropic flows
                    call nsa_stalaw(3,imode,xcons(ndime+1),xphys(ndime+1),xphys(ndime+2),rdumy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))
                 else
                    if (kfl_unkse_nsa < 10) then
                       xhecv           = shecp_nsa(ipoin)   - runiv_nsa / wmean(ipoin,1)
                       xphys(ndime+2)  = (xcons(ndime+2) / xcons(ndime+1) - 0.5_rp * vmodu) / xhecv
                    else if (kfl_unkse_nsa < 20) then
                       xhecv          = shecp_nsa(ipoin)   - runiv_nsa / wmean(ipoin,1)
                       xcons(ndime+2) = xcons(ndime+1) *(xhecv * xphys(ndime+2) + 0.5_rp * vmodu)
                    end if
                 end if

                 if (xphys(ndime+2) < zensa) then
                    xphys(ndime+2) = oltem
                    xhecv          = shecp_nsa(ipoin)   - runiv_nsa / wmean(ipoin,1)
                    xcons(ndime+2) = xcons(ndime+1)*(xhecv * xphys(ndime+2) + 0.5_rp * vmodu)
                 end if

              else if(kfl_fixau(ndime+1) == 2 ) then      ! pressure

                 xphys(ndime+1) = bvess_nsa(ndime+1,ipoin,1)           
                 xphys(ndime+2) = tempe(ipoin,1)

                 call nsa_stalaw(1,imode,xcons(ndime+1),xphys(ndime+1),xphys(ndime+2),rdumy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))

                 oltem= tempe(ipoin,1)
                 do idime=1,ndime
                    xphys(idime) = xcons(idime) / xcons(ndime+1)
                 end do
                 vmodu = xphys(1) * xphys(1) + xphys(2) * xphys(2)
                 if (ndime == 3) vmodu = vmodu + xphys(3) * xphys(3)
                 
                 xhecv          = shecp_nsa(ipoin)   - runiv_nsa / wmean(ipoin,1)                 
                 xcons(ndime+2) = xcons(ndime+1)*(xhecv * xphys(ndime+2) + 0.5_rp * vmodu)
              else if(kfl_fixau(ndime+1) == 3 ) then
                 !
                 ! Prescribe Total Pressure and compute static pressure (and fix it)           
                 ! Done in the continuity equation
                 !
                 adgam = shecp_nsa(ipoin) / xhecv  ! gamma
                 
                 vmodu = xphys(1) * xphys(1) + xphys(2) * xphys(2)
                 if (ndime == 3) vmodu = vmodu + xphys(3) * xphys(3)
!!!    la verdadera uinlet se calcula en una subru ad-hoc
!!!                 uinlet_nsa = sqrt(vmodu)  ! OJO: esto es para probar, pues debería ser una media de la entrada y no la local
                 
                 vauxi= (1.0_rp + (adgam - 1.0_rp) * uinlet_nsa * uinlet_nsa / (2.0_rp * adgam * rgasc_nsa * xphys(ndime+2)))
                 vauxi= vauxi ** (adgam / (adgam - 1.0_rp))
                 vauxi= bvess_nsa(ndime+1,ipoin,1) / vauxi                     ! pressure computed
                 xphys(ndime+1) = vauxi                         ! pressure computed
                 call nsa_stalaw(1,imode,xcons(ndime+1),xphys(ndime+1),xphys(ndime+2),&  
                      rdumy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))                       ! density computed


              end if
              !
              !End b.c. at fixno in position "ndime+1": densi or press
              !

              !
              ! Velocity prescriptions 
              !
              ifstg= 0
              kfixi= 0
              if(kfl_fixau(1)==2)then
                 kfixi= 1
                 call nsa_rotunk( 1,ipoin,ipoin,xadve(1:ndime))
                 do idime=1,ndime
                    if(kfl_fixau(idime) == 2)  then
                       xadve(idime) = 0.0_rp                       ! normal to zero
                    else if(kfl_fixau(idime) == 1)  then
                       xadve(idime) = bvess_nsa(idime,ipoin,1)     ! combined normal with fixity... weird case
                    end if
                 end do
                 call nsa_rotunk(-1,ipoin,ipoin,xadve(1:ndime))
                 do idime= 1,ndime
                    xphys(idime) = xadve(idime) + xvmsh(idime)
                 end do

              else 
                 do idime=1,ndime
                    if(kfl_fixau(idime)==1)then
                       kfixi= 1
                       if (kfl_inlet_nsa /= 0_ip .and. bvess_nsa(idime,ipoin,1) /= 0.0_rp) bvess_nsa(idime,ipoin,1) = ubulk
                       xadve(idime) = bvess_nsa(idime,ipoin,1)                      
                       xphys(idime) = xadve(idime) + xvmsh(idime)
                       ifstg=ifstg+1
                    end if
                 end do
              end if

              if (kfixi == 1) then

                 do idime=1,ndime
                    xcons(idime) = xcons(ndime+1) *  xphys(idime)
                 end do

                 vmodu = xphys(1) * xphys(1) + xphys(2) * xphys(2)
                 if (ndime == 3) vmodu = vmodu + xphys(3) * xphys(3)
                 if ( kfl_isent_nsa == 1) then                  ! isentropic flows
                    call nsa_stalaw(3,imode,xcons(ndime+1),xphys(ndime+1),xphys(ndime+2),rdumy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))
                 else
                    xhecv          = shecp_nsa(ipoin)   - runiv_nsa / wmean(ipoin,1)
                    xphys(ndime+2) = (xcons(ndime+2) / xcons(ndime+1) - 0.5_rp * vmodu) / xhecv
                 end if
                 oltem= xphys(ndime+2)
                 if (xphys(ndime+2) < zensa) then
                    xphys(ndime+2) = oltem
                 end if

              end if

              !
              ! Energy equation prescription:
              !           
              if((kfl_fixau(ndime+2) == 1) .or. (kfl_fixau(ndime+2) == 3)) then
                 vmodu = xphys(1) * xphys(1) + xphys(2) * xphys(2)
                 if (ndime == 3) vmodu = vmodu + xphys(3) * xphys(3)
                 xphys(ndime+2) = bvess_nsa(ndime+2,ipoin,1)
                 if ( kfl_fixau(ndime+2) == 3) xphys(ndime+2) = tstag_nsa  ! stagnation temp prescription
                 xhecv          = shecp_nsa(ipoin)   - runiv_nsa / wmean(ipoin,1)
                 xcons(ndime+2) = xcons(ndime+1) * (xhecv * xphys(ndime+2) + 0.5_rp * vmodu)
              end if

              ! 3. Finally, recompute density (or pressure) from fixed pressure (or density) 
              ! If T has not changed, this is redundant
              if (kfl_fixau(ndime+1) == 1) then
                 ! correct pressure
                 call nsa_stalaw(&
                      2,imode,xcons(ndime+1),xphys(ndime+1),xphys(ndime+2),rdumy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))
              else if ((kfl_fixau(ndime+1) == 2).or.(kfl_fixau(ndime+1) == 3)) then
                 ! correct density
                 call nsa_stalaw(&
                      1,imode,xcons(ndime+1),xphys(ndime+1),xphys(ndime+2),rdumy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))
              end if

           end if

           !
           ! Apply sponge (Rayleigh) if needed:
           !
           sponge: if(kfl_sponge_nsa > 0) then
              if( kfl_benme_nsa == 3 .or. kfl_benme_nsa == 4 .or. kfl_benme_nsa == 5 .or. &
                   kfl_benme_nsa == 5 .or. kfl_benme_nsa == 6 .or. kfl_benme_nsa == 7 .or. &
                   (kfl_benme_nsa >= 200 .and. kfl_sponge_nsa > 0)) then

                 !
                 ! Density:
                 !
                 xcons(ndime+1) =  xcons(ndime+1)*bspon_nsa(1,ipoin) + rekee_nsa(ndime+1,ipoin)*bspon_nsa(2,ipoin)

                 !
                 ! Temperature (potential temperature):
                 !
                 xphys(ndime+2) =  xphys(ndime+2)*bspon_nsa(1,ipoin) + rekee_nsa(ndime+2,ipoin)*bspon_nsa(2,ipoin)

                 !
                 ! Pressure:
                 !
                 call nsa_stalaw(2,imode,xcons(ndime+1),xphys(ndime+1),xphys(ndime+2),rdumy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))

                 do idime=1,ndime 
                    !
                    ! Momentum
                    !
                    xcons(idime) = xcons(idime)*bspon_nsa(1,ipoin) + rekee_nsa(idime,ipoin)*bspon_nsa(2,ipoin)

                    !              
                    !Velocity
                    !
                    xphys(idime) = xcons(idime)/xcons(ndime + 1)
                    vmodu = vmodu + xphys(idime) * xphys(idime)

                 end do

                 !
                 ! Energy:
                 !
                 !xcons(ndime+2) = xcons(ndime+1) * (cvcoe_nsa * xphys(ndime+2) + 0.5_rp * vmodu)
              end if
           end if sponge

           ievat = (ipoin-1)*ndofn_nsa + ndime
           press(ipoin,icomp) = xphys(ndime+1)

           if (kfl_foreg_nsa == 0) then
              tempe(ipoin,icomp) = xphys(ndime+2)

              if(kfl_unkse_nsa < 10) then
                 !Energy
                 unkno(ievat+2) = xcons(ndime+2)
              else if(kfl_unkse_nsa < 20) then !tempe
                 unkno(ievat+2) = xphys(ndime+2)
              end if

           end if
           !Density
           unkno(ievat+1) = xcons(ndime+1)

           if(kfl_ncons_nsa == 0) then
              do idime= 1,ndime
                 ievat = (ipoin-1)*ndofn_nsa + idime
                 unkno(ievat)             = xcons(idime)
                 veloc(idime,ipoin,icomp) = xphys(idime) 
              end do
           else
              do idime= 1,ndime
                 ievat = (ipoin-1)*ndofn_nsa + idime
                 unkno(ievat)             = xphys(idime)
                 veloc(idime,ipoin,icomp) = xphys(idime) 
              end do
           end if

           if (kfl_prope == 0 ) then
              call nsa_lawvis(ipoin,icomp,rdumy,rdumy,rdumy) ! Recalculation of viscosity
           end if

        end do

     end if ! From INOTMASTER

  else if (ktask == four .and. INOTMASTER) then   

     !     ACAAAAAAAAA LO DEL XADVE


     !
     ! called by upcons to compute CONS from PRIMITIVES (phys)
     ! CONS are stored in UNKNO to be recovered by setvar(three)     
     !

     !
     ! Properties from the kernel: molecular weight,Cp & visco
     !
     if (kfl_prope /= 0 ) then
        !call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp_nsa,dummi,dummi,xrano,xrano)
        !call ker_proper('VISCO','NPOIN',dummi,dummi,visco(:,icomp),dummi,dummi,xrano,xrano)
        call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp_nsa)
        call ker_proper('VISCO','NPOIN',dummi,dummi,visco(:,icomp))
        if (kfl_coupl(ID_NASTAL,ID_CHEMIC) == 0 ) then
           wmean = mowei_nsa
        endif
        imode = 1 
     else
        wmean     = mowei_nsa
        shecp_nsa = cpcoe_nsa
     endif

     do ipoin = 1,npoin

        xhecv = shecp_nsa(ipoin) - runiv_nsa / wmean(ipoin,1)

        ievat = (ipoin-1)*ndofn_nsa + ndime
        if (kfl_foreg_nsa == 0) then
           xphys(ndime+2) = unkno(ievat+2)   ! tempe
        end if
        xphys(ndime+1) = unkno(ievat+1)      ! press

        kfl_fixau(1:ndime)=kfl_fixno_nsa(1:ndime,ipoin)
        kfl_fixau(ndime+1)=kfl_fixno_nsa(ndime+1,ipoin)
        kfl_fixau(ndime+2)=kfl_fixno_nsa(ndime+2,ipoin)

        do idime= 1,ndime
           ievat = (ipoin-1)*ndofn_nsa + idime
           xphys(idime) = unkno(ievat)       ! veloc
        end do

        if (kfl_brunt_nsa == 2) then
           stapa_nsa(1) = brunt_nsa(ipoin)
        else
           stapa_nsa(1) = brure_nsa
        end if

        if (kfl_foreg_nsa == 0) then

           ! 1. Compute conservative variables from physical (primitive) ones

           vmodu = xphys(1) * xphys(1) + xphys(2) * xphys(2)
           if (ndime == 3) vmodu = vmodu + xphys(3) * xphys(3)

           olpre= press(ipoin,1) 
           if (xphys(ndime+1) < zensa) then
              xphys(ndime+1) = olpre
           end if

           oltem= tempe(ipoin,1) 
           if (xphys(ndime+2) < zensa) then
              xphys(ndime+2) = oltem
           end if

           call nsa_stalaw(1,imode,xcons(ndime+1),xphys(ndime+1),xphys(ndime+2),rdumy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))
           xcons(ndime+2) =  xcons(ndime+1) *(xhecv * xphys(ndime+2) + 0.5_rp * vmodu)
           do idime=1,ndime
              xcons(idime) =  xcons(ndime+1) * xphys(idime)
           end do

        end if

        ievat = (ipoin-1)*ndofn_nsa + ndime
        if (kfl_foreg_nsa == 0) then
           unkno(ievat+2) = xcons(ndime+2)
        end if
        unkno(ievat+1) = xcons(ndime+1)

        do idime= 1,ndime
           ievat = (ipoin-1)*ndofn_nsa + idime
           unkno(ievat) = xcons(idime)
        end do

     end do

  else if (ktask == five .and. INOTMASTER) then   

     !
     ! called by roback to compute PRIMITIVES from CONSERVATIVES (phys)
     ! after imposing boundary conditions, so no fixno is needed
     ! because roback has leave the conservative variables in unkno
     !
     !
     ! Properties from the kernel: molecular weight,Cp & visco
     !
     if (kfl_prope /= 0 ) then
        !call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp_nsa,dummi,dummi,dummy(:,1),dummy(:,1))
        !call ker_proper('VISCO','NPOIN',dummi,dummi,visco(:,icomp),dummi,dummi,xrano,xrano)
        call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp_nsa)
        call ker_proper('VISCO','NPOIN',dummi,dummi,visco(:,icomp))
        if (kfl_coupl(ID_NASTAL,ID_CHEMIC) == 0 ) then
           wmean = mowei_nsa
        endif
        imode = 1 
     else
        wmean     = mowei_nsa
        shecp_nsa = cpcoe_nsa
     endif

     do ipoin = 1,npoin

        oltem         = tempe(ipoin,1)
        olden         = densi(ipoin,1)
        olpre         = press(ipoin,1)
        olvel(1:ndime)= veloc(1:ndime,ipoin,1)

        kfl_fixau(1:ndime)=kfl_fixno_nsa(1:ndime,ipoin)
        kfl_fixau(ndime+1)=kfl_fixno_nsa(ndime+1,ipoin)
        kfl_fixau(ndime+2)=kfl_fixno_nsa(ndime+2,ipoin)

        ievat = (ipoin-1)*ndofn_nsa + ndime
        xcons(ndime+2) = unkno(ievat+2)   ! energy
        xcons(ndime+1) = unkno(ievat+1)      ! density
        if (xcons(ndime+1) < zensa) then
           xcons(ndime+1) = olden
        end if

        do idime= 1,ndime
           ievat = (ipoin-1)*ndofn_nsa + idime
           xcons(idime) = unkno(ievat)                      ! momentum
           xphys(idime) = xcons(idime) / xcons(ndime+1)     ! velocity
        end do
        vmodu = xphys(1) * xphys(1) + xphys(2) * xphys(2)
        if (ndime == 3) vmodu = vmodu + xphys(3) * xphys(3)

        xhecv            = shecp_nsa(ipoin)   - runiv_nsa / wmean(ipoin,1)
        xphys(ndime+2)   = (xcons(ndime+2)   /  xcons(ndime+1)  - 0.5_rp * vmodu) / xhecv  ! temperature

        if (xphys(ndime+2) < zensa) then
           xphys(ndime+2) = oltem
           xhecv          = shecp_nsa(ipoin)   - runiv_nsa / wmean(ipoin,1)
           xcons(ndime+2) = xcons(ndime+1) *(xhecv * xphys(ndime+2) + 0.5_rp * vmodu)
        end if

        call nsa_stalaw(&
             2,imode,xcons(ndime+1),xphys(ndime+1),xphys(ndime+2),&
             rdumy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))                                  ! pressure

        press(ipoin,icomp)= xphys(ndime+1)
        tempe(ipoin,icomp)= xphys(ndime+2)
        do idime=1,ndime
           veloc(idime,ipoin,icomp)= xphys(idime)
        end do

     end do

  end if


end subroutine nsa_setvar


subroutine nsa_inimet(itask)
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_inimet
  ! NAME 
  !    nsa_inimet
  ! DESCRIPTION
  !    This routine sets up the initial condition for the meteo benchmarks
  ! USED BY
  !    nsa_iniunk
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      mod_postpr
  use      def_nastal
  implicit none
  !Local variables:
  integer(ip) :: itask,idime,icomp,idofn,ipoin,itime,ibubb,nbubb,kshbu
  real(rp)    :: velmi,xfact,xfac2,rgacp,xinve,xradi,xdifi,xrano(3),xtemp,rauxi,exner,&
       dheig,ex1,ex2,ex3,raux2

  real(rp)    :: ptot, dtot, ttot    !total          press, density, and temp
  real(rp)    :: ppert, dpert, tpert !perturbation   press, density, and temp
  real(rp)    :: p_ref, d_ref, t_ref, pi_ref !reference (HS) press, density, and temp

  real(rp) :: dummy
  integer(ip) :: outvar

  if (nzone > 1) call runend("NSA_INIMET: THIS SUB IS NOT PREPARED TO RUN WITH ZONES.")

  if (itask == 0) then 

     !
     ! no initial field given: update boundary conditions  and compute stratified fields
     !
     do ipoin=1,npoin                 

        do idime=1,ndime
           bvess_nsa(idime,ipoin,1) = veloc(idime,ipoin,icomp) 
        end do

        tempe(ipoin,icomp) = tempe_nsa                  
        if (kfl_fixno_nsa(ndofn_nsa,ipoin) > 0) then
           tempe(ipoin,icomp) = bvess_nsa(ndofn_nsa,ipoin,1)
        end if

        if (kfl_brunt_nsa == 2) then
           stapa_nsa(1) = brunt_nsa(ipoin)
        else
           stapa_nsa(1) = brure_nsa
        end if

        !initial pressure
        if(kfl_coupl_nsa == 1) then
           call nsa_stalaw(20,0,densi(ipoin,icomp),press(ipoin,icomp),tempe_nsa,coord(ivert_nsa,ipoin),&
                conce(ipoin,1,icomp),0.0_rp,0.0_rp,0.0_rp)
        else
           call nsa_stalaw(20,0,densi(ipoin,icomp),press(ipoin,icomp),tempe_nsa,coord(ivert_nsa,ipoin),&
                conce(1,1,1),0.0_rp,0.0_rp,0.0_rp)
        end if
        !derive initial density from pressure and initial potential temperature (constant)
        if(kfl_coupl_nsa == 1) then
           call nsa_stalaw(-1,0,densi(ipoin,icomp),press(ipoin,icomp),tempe_nsa,coord(ivert_nsa,ipoin),&
                conce(ipoin,1,icomp),0.0_rp,0.0_rp,0.0_rp)
        else
           call nsa_stalaw(-1,0,densi(ipoin,icomp),press(ipoin,icomp),tempe_nsa,coord(ivert_nsa,ipoin),&
                conce(1,1,1),0.0_rp,0.0_rp,0.0_rp) 
        end if

        bvess_nsa(ndime+1,ipoin,1) = densi(ipoin,icomp)                     
        if(kfl_fixno_nsa(ndime+1,ipoin) == 2 ) then           ! pressure-prescribed node
           bvess_nsa(ndime+1,ipoin,1) = press(ipoin,icomp) 
        end if

        bvess_nsa(ndime+2,ipoin,1) = tempe(ipoin,icomp)

        umome(1:ndime,ipoin,icomp) =  densi(ipoin,icomp) &
             * veloc(1:ndime,ipoin,icomp)
        energ(ipoin,icomp) =  densi(ipoin,icomp) &
             * (cvcoe_nsa * tempe(ipoin,icomp)  + 0.5*velmi*velmi)        

     end do

  else if (itask > 0) then

     dummy = 0.0_rp
     !
     ! Initialize variables from reference values
     !        
     icomp=min(3_ip,ncomp_nsa)
     stapa_nsa = 0.0_rp

     xfact    = - grnor_nsa / rgasc_nsa / tempe_nsa

     if (kfl_inico_nsa==0) then
        do ipoin=1,npoin
           veloc(1:ndime,ipoin,icomp) = 0.0_rp
        end do
     else if (kfl_inico_nsa==1) then 
        !
        ! kfl_inico_nsa = 1 when 
        !   INITIAL_CONDS:     REFERENCE_VALUES
        !
        velmi = 0.0_rp
        do idime=1,ndime
           do ipoin=1,npoin
              veloc(idime,ipoin,icomp) = veloc_nsa(idime)        
           end do
           velmi = velmi + veloc_nsa(idime) * veloc_nsa(idime)
        end do
        velmi = sqrt(velmi)
     end if

     !
     ! Initial computed in nsa_initial_conditions, 
     !
     call nsa_initial_conditions(outvar)
     
     !
     ! Store the variables into press, densi, tempe, veloc, umome, energ
     !
     do ipoin = 1,npoin
        !
        ! Perturbation Velocity: OK
        !
        do idime=1,ndime
           veloc(idime,ipoin,icomp) = bvess_nsa(idime,ipoin,1)
        end do

        !
        ! Total variables:
        !
        ! bvess in initial_conditions is total variables
        ttot = bvess_nsa(ndofn_nsa, ipoin,1)
        if(outvar == 1) then
           dtot = bvess_nsa(ndofn_nsa-1, ipoin,1)
           call nsa_stalaw(2,0,dtot,ptot,ttot,coord(ivert_nsa,ipoin),dummy,dummy,0.0_rp,0.0_rp)
        else if(outvar == 2) then
           ptot = bvess_nsa(ndofn_nsa-1, ipoin,1)
           call nsa_stalaw(1,0,dtot,ptot,ttot,coord(ivert_nsa,ipoin),dummy,dummy,0.0_rp,0.0_rp)
        end if

        !
        ! Reference (HS) variables: !OK
        !
        d_ref = rekee_nsa(ndime+1, ipoin) !ok
        t_ref = rekee_nsa(ndime+2, ipoin) !ok
        p_ref = rekee_nsa(ndime+3, ipoin) !ok

        !Populate pert press(), tempe(), densi():
        press(ipoin,icomp) = ptot
        tempe(ipoin,icomp) = ttot
        densi(ipoin,icomp) = dtot

        !Momentum and energy are computed from total variables
        !(see equations for meteorology to understand why):
        umome(1:ndime,ipoin,icomp) =  dtot &
             * veloc(1:ndime,ipoin,icomp)

        !eneg is not really needed for our problems
        !but we initialize it just in case
        energ(ipoin,icomp) =  dtot &
             * (cvcoe_nsa * ttot*(pbaro_nsa/ptot)**(-rgasc_nsa/cpcoe_nsa) &
             + 0.5*velmi*velmi)
     end do
     
     !if(ndime < 3) then
     !   call nsa_outvtk('ALYA_INITIAL_dynamics_SVN.vtk')
     !else
     !   call nsa_outvtk_3d('ALYA_INITIAL_dynamics_SVN_3D.vtk')
     !end if

  end if

end subroutine nsa_inimet

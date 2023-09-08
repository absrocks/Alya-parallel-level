subroutine tem_iniunk()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_iniunk
  ! NAME 
  !    tem_iniunk
  ! DESCRIPTION
  !    This routine sets up the initial condition for the temperature.
  ! USED BY
  !    tem_begste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use def_kermod
  use mod_ker_proper 
  use mod_ker_space_time_function
  use mod_chktyp, only : check_type
  use mod_communications, only: PAR_MAX, PAR_MIN
  use mod_ADR, only : ADR_manufactured_nodal_solution
  use mod_ADR, only : ADR_manufactured_error

  implicit none
  integer(ip) :: ipoin,icomp,ifunc,ivalu,itime,ielem,inode
  integer(ip) :: kfl_advec_old,kfl_timei_old, kfl_value
  real(rp)    :: dtinv_tmp, cploc(6,2), dummr, tenew

  if( kfl_rstar == 0 ) then 
 
     !-------------------------------------------------------------------
     !
     ! Normal run
     !
     !-------------------------------------------------------------------

     if( INOTMASTER ) then

        icomp = min(3_ip,ncomp_tem)
        if( kfl_inico_tem < - 100 ) then 
           !
           ! Take initial condition from a value function field
           ! 
           kfl_value = -  kfl_inico_tem - 100
           call check_type(xfiel,kfl_value,1_ip,nunkn_tem) ! Check if value function exist
           do ipoin = 1,nunkn_tem
              therm(ipoin,icomp) = xfiel(kfl_value) % a(1,ipoin,1)
           end do

        else if( kfl_inico_tem < 0 ) then 
           ! 
           ! Space time function
           !
           ifunc = -kfl_inico_tem
           do ipoin = 1,nunkn_tem
              call ker_space_time_function(&
                   ifunc,coord(1:ndime,ipoin),cutim,therm(ipoin,icomp))
           end do
           
        else if( kfl_inico_tem == 1 ) then
           
           do ipoin = 1,nunkn_tem
              therm(ipoin,icomp) = initial_tem
           end do

        else

           do ipoin = 1,nunkn_tem              
              therm(ipoin,icomp) = bvess_tem(1,ipoin,1)  
           end do
           
        end if
        !
        ! Smooth Dirichlet boundary conditions
        !
        call tem_smobcs()  !Currently it's not doing anything 
        !
        ! Functions and specific boundary conditions
        !
        call tem_updbcs(-1_ip)
        !
        ! When solving an exact solution, impose bc
        !
        if( kfl_exacs_tem /= 0 ) then
           !
           ! Manufactured solution
           !
           icomp = min(3_ip,ncomp_tem)
           !call ADR_manufactured_nodal_solution(ADR_tem,cutim,therm(:,icomp))
           !call ADR_manufactured_error(ADR_tem,ittim,cutim,therm(:,icomp))

        else
           !
           ! Dirichlet conditions for the temperature
           !
           if( kfl_discr_tem == NODAL_SCHEME ) then
              posit_tem = -10e10_rp
              negat_tem =  10e10_rp 
              do ipoin = 1,nunkn_tem
                 if( kfl_fixno_tem(1,ipoin) > 0 .and. kfl_regim_tem==4) then
                    do ivalu = 1,6
                       cploc(ivalu,1) = sphec(ipoin,ivalu,1)
                       cploc(ivalu,2) = sphec(ipoin,ivalu,2)
                    end do
                    dummr = 0.0_rp 
                    call tem_comput(2_ip,bvess_tem(1,ipoin,2),dummr,cploc,tenew)
                    bvess_tem(1,ipoin,1) = tenew
                    posit_tem = max( posit_tem, tenew)
                    negat_tem = min( negat_tem, tenew)
                    therm(ipoin,icomp) =  bvess_tem(1,ipoin,1)
                 else if( kfl_fixno_tem(1,ipoin) > 0 ) then  
                    therm(ipoin,icomp) = bvess_tem(1,ipoin,1)  
                 end if
              end do
              if (kfl_posit_tem == 1) call PAR_MAX(posit_tem,'IN MY CODE')
              if (kfl_negat_tem == 1) call PAR_MIN(negat_tem,'IN MY CODE')
           end if

        end if
        !
        ! Initial value is current value
        !
        do ipoin = 1,nunkn_tem
           therm(ipoin,1) = therm(ipoin,icomp)
        end do

        if (kfl_regim_tem == 4_ip) then
           ! 
           ! Compute current temperature from enthalpy 
           !
           call tem_comtem
           ! 
           ! Store current temperature 
           !
           do ipoin=1,nunkn_tem
              tempe(ipoin,icomp) = tempe(ipoin,1)
           end do

        end if
        !
        ! Latex output format
        !
        call tem_outlat(one)

     else 

        if (kfl_posit_tem == 1) call PAR_MAX(posit_tem,'IN MY CODE')
        if (kfl_negat_tem == 1) call PAR_MIN(negat_tem,'IN MY CODE')

     end if
     !
     ! Solve diffusion problem
     !
     if( kfl_inico_tem == 3 ) then
        kfl_advec_old = kfl_advec_tem
        kfl_timei_old = kfl_timei_tem
        dtinv_tmp     = dtinv_tem
        kfl_advec_tem = 0
        kfl_timei_tem = 0
        dtinv_tem     = 0.0_rp
        call tem_updunk(8_ip)
        if( kfl_prope /= 0 ) call ker_updpro() ! Force update of kernel properties
        call tem_solite()
        call tem_updunk(9_ip)
        kfl_inidi_tem = 0
        kfl_advec_tem = kfl_advec_old 
        kfl_timei_tem = kfl_timei_old 
        dtinv_tem     = dtinv_tmp       
     end if

  else

     !-------------------------------------------------------------------
     !
     ! Read restart file
     !
     !-------------------------------------------------------------------

     call tem_restar(1_ip)     
     call tem_updunk(6_ip)
     !
     ! Functions and specific boundary conditions
     !
     call tem_updbcs(-1_ip)

     if (kfl_regim_tem == 4) then
        ! 
        ! Compute current temperature from enthalpy 
        !
        call tem_comtem
        ! 
        ! Store current temperature 
        !
        do ipoin=1,nunkn_tem
           tempe(ipoin,3) = tempe(ipoin,1)
        end do
        ! 
        ! High-order temporal schemes 
        !
        if(kfl_tisch_tem==2) then
           !
           ! BDF scheme
           !
           do ipoin=1,nunkn_tem
              do itime=2+kfl_tiaor_tem,4,-1
                 tempe(ipoin,itime) = tempe(ipoin,itime-1)
              end do
           end do
        end if

     end if

     avtim_tem = cutim ! Accumulated time for time-averaging variables 

     if( kfl_discr_tem == NODAL_SCHEME ) then
        !
        ! Dicichlet conditions for the temperature
        !
        icomp     = min(3_ip,ncomp_tem)
        posit_tem = -10e10_rp
        negat_tem =  10e10_rp
        do ipoin = 1,nunkn_tem
           if( kfl_fixno_tem(1,ipoin) > 0 .and. kfl_regim_tem==4) then
              do ivalu = 1,6
                 cploc(ivalu,1) = sphec(ipoin,ivalu,1)
                 cploc(ivalu,2) = sphec(ipoin,ivalu,2)
              end do
              dummr = 0.0_rp
              call tem_comput(2_ip,bvess_tem(1,ipoin,2),dummr,cploc,tenew)
              bvess_tem(1,ipoin,1) = tenew
              posit_tem = max( posit_tem, tenew)
              negat_tem = min( negat_tem, tenew)
              therm(ipoin,icomp) =  bvess_tem(1,ipoin,1)
           else if( kfl_fixno_tem(1,ipoin) > 0 ) then
              therm(ipoin,icomp) = bvess_tem(1,ipoin,1)
           else
              therm(ipoin,icomp) = bvess_tem(1,ipoin,1)              
           end if
        end do
        if (kfl_posit_tem == 1) call PAR_MAX(posit_tem,'IN MY CODE')
        if (kfl_negat_tem == 1) call PAR_MIN(negat_tem,'IN MY CODE')

     end if

  end if
  !-------------------------------------------------------------------
  !
  ! Interpolation from coarse to fine mesh
  !
  !-------------------------------------------------------------------

  if(kfl_meshi_tem /= 0_ip) call tem_coarfine(1_ip)

  call tem_coupli(ITASK_INIUNK)

end subroutine tem_iniunk

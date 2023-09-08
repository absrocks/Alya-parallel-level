subroutine tem_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_inivar
  ! NAME 
  !    tem_inicar
  ! DESCRIPTION
  !    This routine initializes some variables
  !    ITASK=1 ... When starting the run (from Turnon)
  !    ITASK=2 ... First time step. This is needed as some variables 
  !                are not initialized before
  ! USES
  ! USED BY
  !    tem_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use def_solver
  implicit none
  integer(ip), intent(in) :: itask

  select case(itask)

  case(0_ip)
     !
     ! Postprocess
     !
     postp(1) % wopos (1, 1) = 'TEMPE'
     postp(1) % wopos (1, 2) = 'HEATF'
     postp(1) % wopos (1, 3) = 'TESGS'
     postp(1) % wopos (1, 4) = 'ERROR'
     postp(1) % wopos (1, 5) = 'AVTEM' ! average temp
     postp(1) % wopos (1, 6) = 'VELOC'
     postp(1) % wopos (1, 7) = 'TURVI'
     postp(1) % wopos (1, 8) = 'RESID'
     postp(1) % wopos (1, 9) = 'GROUP'
     postp(1) % wopos (1,10) = 'PROJ1'
     postp(1) % wopos (1,11) = 'LIMIT'
     postp(1) % wopos (1,12) = 'LINTE'
     postp(1) % wopos (1,13) = 'TESTS'
     postp(1) % wopos (1,14) = 'WATVA'
     postp(1) % wopos (1,15) = 'FIXTE' ! I have changed the name otherwise it will overlap with the one for vel
     postp(1) % wopos (1,16) = 'TEMP2'
     postp(1) % wopos (1,17) = 'TESG2'    
     postp(1) % wopos (1,18) = 'AVTE2' ! average tempe*tempe
     postp(1) % wopos (1,19) = 'AVTEV' ! average veloc*tempe
     postp(1) % wopos (1,20) = 'AVDEN' ! average density
     postp(1) % wopos (1,21) = 'FVVEL' ! Favre average velocity
     postp(1) % wopos (1,22) = 'HEATN' ! Favre average velocity
     postp(1) % wopos (1,23) = 'GRATE' ! Temperature gradient
     postp(1) % wopos (1,24) = 'CHARA' 
     postp(1) % wopos (1,25) = 'BVNAT' 
     postp(1) % wopos (1,26) = 'ENTHA' ! Enthalpy
     postp(1) % wopos (1,27) = 'REACT' !solve_sol(1) % lpoin_reaction
     postp(1) % wopos (1,28) = 'TFLUX'
     postp(1) % wopos (1,29) = 'RESID'
     postp(1) % wopos (1,30) = 'PROJ2'
     postp(1) % wopos (1,31) = 'RESHE' ! Heat flux using residuals
     postp(1) % wopos (1,32) = 'AVRES' ! average Heat flux using residuals

     postp(1) % wopos (2, 1) = 'SCALA'
     postp(1) % wopos (2, 2) = 'SCALA'
     postp(1) % wopos (2, 3) = 'SCALA'
     postp(1) % wopos (2, 4) = 'SCALA'
     postp(1) % wopos (2, 5) = 'SCALA'
     postp(1) % wopos (2, 6) = 'VECTO'
     postp(1) % wopos (2, 7) = 'SCALA'
     postp(1) % wopos (2, 8) = 'SCALA'
     postp(1) % wopos (2, 9) = 'SCALA'
     postp(1) % wopos (2,10) = 'SCALA'
     postp(1) % wopos (2,11) = 'SCALA' 
     postp(1) % wopos (2,12) = 'SCALA' 
     postp(1) % wopos (2,13) = 'VECTO' 
     postp(1) % wopos (2,14) = 'SCALA' 
     postp(1) % wopos (2,15) = 'SCALA' 
     postp(1) % wopos (2,16) = 'SCALA' 
     postp(1) % wopos (2,17) = 'R3PVE'    
     postp(1) % wopos (2,18) = 'SCALA'
     postp(1) % wopos (2,19) = 'VECTO'
     postp(1) % wopos (2,20) = 'SCALA'
     postp(1) % wopos (2,21) = 'VECTO'
     postp(1) % wopos (2,22) = 'SCALA'
     postp(1) % wopos (2,23) = 'VECTO'
     postp(1) % wopos (2,24) = 'VECTO'
     postp(1) % wopos (2,25) = 'SCALA'
     postp(1) % wopos (2,26) = 'SCALA'
     postp(1) % wopos (2,27) = 'SCALA' !solve_sol(1) % lpoin_reaction
     postp(1) % wopos (2,28) = 'SCALA'
     postp(1) % wopos (2,29) = 'SCALA'
     postp(1) % wopos (2,30) = 'SCALA'
     postp(1) % wopos (2,31) = 'SCALA'
     postp(1) % wopos (2,32) = 'SCALA'

     postp(1) % wopos (1,33) = 'NFIXN'
     postp(1) % wopos (2,33) = 'SCALA'

     postp(1) % wopos (1,34) = 'TOUCH'
     postp(1) % wopos (2,34) = 'SCALA'
     !
     ! Set and witness variables 
     !
     postp(1) % woese(1)     = 'MEANC'
     postp(1) % wobse(1)     = 'MEANT'  
     postp(1) % wobse(2)     = 'MEANH'  
     postp(1) % wobse(3)     = 'MHEAN'  
     postp(1) % wonse(1)     = 'TEMPE'  
     postp(1) % wowit(1)     = 'TEMPE'
     postp(1) % wowit(2)     = 'FLUXX'
     postp(1) % wowit(3)     = 'FLUXY'
     postp(1) % wowit(4)     = 'FLUXZ'
     !
     ! Solver
     !     
     call soldef(-1_ip)
     solve(1) % wprob     = 'TEMPERATURE'          ! Equation name
     solve(1) % kfl_solve = 1                      ! Output flag
     !
     ! Nullify pointers
     !
     nullify(kfl_funty_tem)   
     nullify(funpa_tem)     
     nullify(tncod_tem)       
     nullify(tgcod_tem)       
     nullify(tbcod_tem)       
     nullify(kfl_fixno_tem)
     nullify(kfl_fixbo_tem)  
     nullify(kfl_funno_tem)  
     nullify(kfl_funbo_tem)   
     nullify(bvess_tem)  
     nullify(bvnat_tem)   
     nullify(fixnb_tem)        
     nullify(lmatb_tem)      
     nullify(gradc_tem)
     nullify(power_tem)      
     nullify(avtem_tem)      
     nullify(avres_tem)
     nullify(grtem_tem)    
     nullify(teold_tem)     
     nullify(avte2_tem) 
     nullify(avtev_tem)     
     nullify(avden_tem)      
     nullify(fvvel_tem)     


  case(1_ip)

     if(kfl_timei_tem==0) then                     ! Time integration
        dtinv_tem=1.0_rp
     else
        kfl_timei=1
     end if
     kfl_stead_tem=0

     if(kfl_cotur_tem /=0) then  ! Recommended not to use mu_t gradients 
        kfl_grdif_tem=0
     else 
        kfl_grdif_tem=1          ! It is recommnded for laminar flows
     end if

     kfl_tiaor_tem=kfl_tiacc_tem                    ! Time accuracy: save original value
     if(kfl_advec_tem==0) bemol_tem=0.0_rp          ! There is no advection: BEMOL_TEM=0
     if(kfl_timei_tem==1) then                      ! Number of velocity components
        if(kfl_tisch_tem==1) then
           ncomp_tem=3                              ! Trapezoidal rule
        else if(kfl_tisch_tem==2) then
           ncomp_tem=2+kfl_tiacc_tem                ! BDF scheme
        end if
     else
        ncomp_tem = 2     
     end if

     ittot_tem = 0                                  ! Others
     dpthe     = 0.0_rp                             ! Low-Mach: dp0/dt
     !
     ! Solver fixity
     !
     if( INOTMASTER ) then        
        solve(1) % bvess     => bvess_tem(:,:,1)
        solve(1) % kfl_fixno => kfl_fixno_tem
     end if
     xmass_tem =0.0_rp                              ! Low-Mach: mass
     kfl_rstar_two= .false.                         ! restarting bdf file?     
     !
     ! Finite volume
     !
     if(      kfl_discr_tem == 0 ) then
        solve(1) % kfl_where = SOL_NODES
        nunkn_tem            = npoin
     else if( kfl_discr_tem == 1 ) then
        solve(1) % kfl_where = SOL_ELEMENTS
        nunkn_tem            = nelem
     end if
     !
     ! Heat flux field
     !
     if (kfl_flux_tem > 0_ip) then
        heat_flux => xfiel(kfl_flux_tem) % a(:,:,1) 
     endif

  case(2_ip)
     !
     ! Velocity subgrid scale
     !
     if( associated(vesgs) .and. kfl_advec_tem == 1 ) then
        kfl_sgsve_tem = 1
     else
        kfl_sgsve_tem = 0
     end if

  end select

end subroutine tem_inivar

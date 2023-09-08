subroutine got_reanut()
  !-----------------------------------------------------------------------
  !****f* Gotita/got_reanut
  ! NAME 
  !    got_reanut
  ! DESCRIPTION
  !    This routine reads the numerical treatment for the incomcdropible
  !    NAvier-Stokes equations.
  ! INPUT
  ! OUTPUT
  !    staco_got ....... Stabilization parameters (default=4,2,1,1)
  !    kfl_preco_got ... Preconditioner
  !                      1= LEFT: sqrt(A)^{-1/2}, RIGHT: sqrt(A)^{1/2} 
  !                      2= LEFT: diag(A)^{-1} constricted by the solver
  !                      3= LEFT: (tau U,V) approx (L^{-1}(U),V)  
  !                      4= Same as 2 but it is constructed explicitly
  !                         by GOTITA  
  ! USES
  !    ecoute
  ! USED BY
  !    got_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_solver
  use def_gotita
  use def_domain
  use mod_memchk
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: istab

  if(kfl_paral<=0) then
     !
     !  Initializations (defaults)
     !
     kfl_algor_got = 1               ! Monolithic
     kfl_artif_got = 0               ! No artificial viscosity
     kfl_coupl_got = 1               ! Coupling on
     kfl_dttyp_got = 1               ! Module time step strategy
     kfl_ellen_got = 0               ! Minimum element length
     kfl_linea_got = 1               ! Linearization (RHS=0, Picard=1, Newton=2)
     kfl_normc_got = 2               ! L2 norm for convergence
     kfl_penal_got = 0               ! Penalization
     kfl_sgsco_got = 0               ! Stabilization strategy
     kfl_sgsti_got = 0               ! Stabilization strategy
     kfl_shocc_got = 0               ! Shock capturing continuity off
     kfl_shocm_got = 0               ! Shock capturing momentum off
     kfl_staty_got = 1               ! Stabilizatiom type
     kfl_taust_got = 2               ! Tau strategy
     kfl_tiacc_got = 1               ! First order time integ.
     kfl_weigc_got = 0               ! Weight continuity residual 
     kfl_weigm_got = 0               ! Weight momentum residual 
     itart_got     = 1e6             ! Stop artificial viscosity
     itshc_got     = 1e6             ! Stop shock capturing
     itshm_got     = 1e6             ! Stop shock capturing
     mibgs_got     = 1               ! # iterations Block Gauss-Seidel
     miinn_got     = 1               ! One internal iteration
     misgs_got     = 1               ! Subgrid scale iterations
     neule_got     = 1               ! Euler iterations     
     npica_got     = 1               ! Picard iterations
     artif_got(1)  = 0.0_rp          ! Artificial viscosity (momentum)
     artif_got(2)  = 0.0_rp          ! Artificial viscosity (continuity)
     cotol_got     = 1.0_rp          ! Internal tolerance
     cutof_got     = 0.0_rp          ! Alpha relative cut-off
     penal_got     = 0.0_rp          ! Penalization
     relax_got     = 1.0_rp          ! GOTITA relaxation
     relgs_got     = 1.0             ! Relaxation Block Gauss-Seidel
     relsg_got     = 1.0_rp          ! Subgrid scale relaxation
     safet_got     = 1.0e10          ! Safety factor
     shock_got     = 0.0_rp          ! Shock capturing value
     sstol_got     = 1.0e-5          ! Steady-satate tolerance
     staco_got(1)  = 1.0_rp          ! Diffusive term
     staco_got(2)  = 1.0_rp          ! Convective term
     staco_got(3)  = 1.0_rp          ! Reactive term
     tobgs_got     = 1e-6            ! Tolerance Block Gauss-Seidel
     tosgs_got     = 1.0e-2          ! Subgrid scale tolerance
     xmaxi_got(1)  =  1.0_rp         ! Airfoil bounding box
     xmaxi_got(2)  =  1.0_rp         ! Airfoil bounding box
     xmaxi_got(3)  =  1.0_rp         ! Airfoil bounding box
     xmini_got(1)  = -1.0_rp         ! Airfoil bounding box
     xmini_got(2)  = -1.0_rp         ! Airfoil bounding box
     xmini_got(3)  = -1.0_rp         ! Airfoil bounding box
     !
     ! Reach the section
     !
     ivari_got = 1              
     call ecoute('got_reanut')
     do while(words(1)/='NUMER')
        call ecoute('got_reanut')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDNU')
        call ecoute('got_reanut')

        if(words(1)=='TYPEO') then 
           if(words(2)=='SUPG ') then
              kfl_staty_got=1
           else if(words(2)=='ASGS '.or.words(2)=='ASGS1') then
              kfl_staty_got=2
           else if(words(2)=='ASGS2') then
              kfl_staty_got=3
           else if(words(2)=='OFF  ') then
              kfl_staty_got=0
           end if

        else if(words(1)=='TAUST') then
           if(words(2)=='OFF  ') then
              kfl_taust_got=0 
           else if(words(2)=='CODIN') then
              kfl_taust_got=1
           else if(words(2)=='EXACT') then
              kfl_taust_got=2
           else if(words(2)=='SHAKI') then
              kfl_taust_got=3
           else if(words(2)=='HAUKE') then
              kfl_taust_got=4
           end if

        else if(words(1)=='COUPL') then 
           if(words(2)=='ON   ') then
              kfl_coupl_got=1
              kfl_penal_got=1
              penal_got=getrea('PENAL',1.0e-3_rp,'#Penalization factor')
           else
              kfl_coupl_got=0
           end if

        else if(words(1)=='TRACK') then 
           if(exists('CONVE')) then
              kfl_sgsco_got=1 
              misgs_got=getint('ITERA',1_ip,  '#Subgrid scale iterations')
              tosgs_got=getrea('TOLER',1.0e-4_rp,'#Subgrid scale Tolerance')
              relsg_got=getrea('RELAX',1.0_rp,'#Subgrid scale Relaxation') 
           end if
           if(exists('TIME ')) kfl_sgsti_got=1 

        else if(words(1)=='BOUND') then 
           xmini_got(1)=getrea('XMIN ',-1.0_rp,'#Bounding box xmin')
           xmini_got(2)=getrea('YMIN ',-1.0_rp,'#Bounding box ymin')
           xmini_got(3)=getrea('ZMIN ',-1.0_rp,'#Bounding box zmin')
           xmaxi_got(1)=getrea('XMAX ', 1.0_rp,'#Bounding box xmax')
           xmaxi_got(2)=getrea('YMAX ', 1.0_rp,'#Bounding box ymax')
           xmaxi_got(3)=getrea('ZMAX ', 1.0_rp,'#Bounding box zmax')

        else if(words(1)=='LINEA') then
           if(words(2)=='PICAR') then
              kfl_linea_got = 1
           else if(words(2)=='NEWTO') then
              kfl_linea_got = 2
              npica_got     = getint('PICAR',0_ip,'#Number of Picard iterations')
           end if

        else if(words(1)=='PENAL') then 
           if(words(2)=='ON   ') then
              kfl_penal_got=1
              penal_got=getrea('VALUE',1.0e-3_rp,'#Penalization factor')
           else if(words(2)=='AUTOM') then
              kfl_penal_got=2
           end if

        else if(words(1)=='CUTOF') then 
           cutof_got=getrea('CUTOF',1.0e-3_rp,'#Alpha relative cut-off value')

        else if(words(1)=='ARTIF') then 
           if(words(2)=='ON   '.or.words(2)=='CLASS') then
              kfl_artif_got=1
           else if(words(2)=='ITERA') then
              kfl_artif_got=2
           end if
           if(exists('STOPA')) then
              itart_got=getint('STOPA',100000_ip,'#Stop artificial viscosity')
           end if
           artif_got(1)=getrea('MOMEN',1.0e-4_rp,'#Momentum artificial viscosity')
           artif_got(2)=getrea('CONTI',1.0e-4_rp,'#Continuity artificial viscosity')

        else if(words(1)=='TIMES') then
           if(words(2)=='LOCAL') then
              kfl_dttyp_got=1
           end if

        else if(words(1)=='SHOCK') then
           if(exists('MOMEN')) then
              if(exists('ISOTR').or.exists('ON   ')) then
                 kfl_shocm_got = 1
                 shock_got     = getrea('VALUE',0.0_rp,'#Shock capturing parameter')
                 if(getcha('RESID','ALL  ','#Residual')=='ALL  ') kfl_weigm_got=1
              else if(exists('ANISO')) then
                 kfl_shocm_got = 2
                 shock_got     = getrea('VALUE',0.0_rp,'#Shock capturing parameter')
                 if(getcha('RESID','ALL  ','#Residual')=='ALL  ') kfl_weigm_got=1
              end if
              if(exists('STOPA')) then
                 itshm_got=getint('STOPA',100000_ip,'#Stop shock capturing')
              end if
           else
              if(exists('ISOTR').or.exists('ON   ')) then
                 kfl_shocc_got = 1
                 shock_got     = getrea('VALUE',0.0_rp,'#Shock capturing parameter')
                 if(getcha('RESID','ALL  ','#Residual')=='ALL  ') kfl_weigc_got=1
              else if(exists('ANISO')) then
                 kfl_shocc_got = 2
                 shock_got     = getrea('VALUE',0.0_rp,'#Shock capturing parameter')
                 if(getcha('RESID','ALL  ','#Residual')=='ALL  ') kfl_weigc_got=1
              end if
              if(exists('STOPA')) then
                 itshc_got=getint('STOPA',100000_ip,'#Stop shock capturing')
              end if
           end if

        else if(words(1)=='ELEME') then
           call realen(kfl_ellen_got)

        else if(words(1)=='TIMEA') then
           kfl_tiacc_got = getint('TIMEA',1_ip,'#Time integration order')
           neule_got     = getint('EULER',0_ip,'#EULER TIME STEPS')

        else if(words(1)=='EULER') then
           neule_got=getint('EULER',0_ip,'#EULER TIME STEPS')

        else if(words(1)=='TIMEI') then
           kfl_tiacc_got = getint('ORDER',1_ip,'#Time integration order')
           neule_got     = getint('EULER',0_ip,'#EULER TIME STEPS')

        else if(words(1)=='SAFET') then
           safet_got = param(1)

        else if(words(1)=='STEAD') then
           sstol_got = param(1)

        else if(words(1)=='NORMO') then
           if(exists('L1   ')) then
              kfl_normc_got = 1
           else if(exists('L-inf')) then
              kfl_normc_got = 0
           end if

        else if(words(1)=='MAXIM') then
           miinn_got = int(param(1))

        else if(words(1)=='CONVE') then
           cotol_got = param(1)

        else if(words(1)=='MOMEN') then
           ivari_got=1
        else if(words(1)=='CONTI') then
           ivari_got=2

        else if(words(1)=='ALGEB') then
           solve_sol => solve(ivari_got:)
           call reasol(1_ip)

        else if(words(1)=='ALGOR') then
           if(exists('MONOL')) kfl_algor_got =  1
           if(exists('BLOCK')) kfl_algor_got =  4
           if(kfl_algor_got==4) then
              mibgs_got = getint('ITERA',1_ip,      '#Number of Block Gauss-Seidel iterations')
              tobgs_got = getrea('TOLER',0.01_rp,'#Tolerance Block Gauss-Seidel')
              relgs_got = getrea('RELAX',1.0_rp, '#Relaxation Block Gauss-Seidel')
           end if

        else if(words(1)=='PRECO') then 
           solve_sol => solve(ivari_got:)
           call reasol(2_ip)

        else if(words(1)=='RELAX') then
           relax_got = param(1)

        end if
     end do

  end if

end subroutine got_reanut

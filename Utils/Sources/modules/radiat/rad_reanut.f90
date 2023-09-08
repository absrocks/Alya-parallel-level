subroutine rad_reanut
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_reanut
  ! NAME 
  !    rad_reanut
  ! DESCRIPTION
  !    This routine reads the numerical treatment for RADIAT module
  ! USES
  !    ecoute
  ! USED BY
  !    rad_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_solver
  use def_radiat
  use def_domain
  use mod_ecoute, only :  ecoute
  implicit none

  if(kfl_paral<=0) then
     ! 
     !  Initializations (defaults)
     !
     kfl_ellen_rad = 0                                ! Minimum element length

     kfl_sgsno_rad = 0                                ! Subgrid scale non-lniear tracking
     kfl_taust_rad = 1                                ! Tau calculation option
     kfl_ortho_rad = -2                               ! Stabilization default OFF (Galerkin)
     kfl_limit_rad = 0                                ! Lmiter
     kfl_bubbl_rad = 0                                ! Bubble function
     kfl_assem_rad = 0                                ! Assembly type
     kfl_normc_rad = 2                                ! L2 norm for convergence
     miinn_rad     = 10                               ! Internal iterations
     misgs_rad     = 1                                ! Max # of SGS iterations

     staco_rad(1)  = 1.0_rp                           ! Diffusive term
     staco_rad(2)  = 1.0_rp                           ! Convective term
     staco_rad(3)  = 1.0_rp                           ! Reactive term
     cotol_rad     = 0.05_rp                           ! Internal tolerance
     relsg_rad     = 1.0_rp                           ! Relaxation parameter of subgrid scale
     tosgs_rad     = 0.01_rp                          ! Subgrid scale tolerance
     solve_sol     => solve                       ! Solver type
     !
     ! Reach the section
     !
     call ecoute('rad_reanut')
     do while(words(1)/='NUMER')
        call ecoute('rad_reanut')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDNU')
        call ecoute('rad_reanut')

        if(words(1)=='TRACK') then
           if(exists('NONLI')) then
              kfl_sgsno_rad=1
              misgs_rad=getint('ITERA',1_ip,   '#Subgrid scale iterations')
              tosgs_rad=getrea('TOLER',0.01_rp,'#Subgrid scale Tolerance')
              relsg_rad=getrea('RELAX',1.0_rp, '#Subgrid scale Relaxation') 
           end if

        else if(words(1)=='TAUST') then
           call reatau(kfl_taust_rad)

        else if(words(1)=='STRAT' .or. words(1) == 'STABI' ) then
           if( words(2) == 'SU   ' .or. words(2) == 'FIRST' ) then
              kfl_ortho_rad = -1
           else if( words(2) == 'ASGS ' ) then
              kfl_ortho_rad =  0
           else if( words(2) == 'FULLO' ) then
              kfl_ortho_rad =  1
           else if( words(2) == 'OSS  ' ) then
              kfl_ortho_rad =  2
              if( exists('NOLIM') ) kfl_limit_rad =  0  ! No limiter
              if( exists('SOTO ') ) kfl_limit_rad =  1  ! Soto limiter
              if( exists('DIFFU') ) kfl_limit_rad =  2  ! Very diffusive limiter
              if( exists('FIRST') ) kfl_limit_rad = -1  ! First order
           else if( words(2) == 'OFF  ' .or.  words(2) == 'GALER' ) then
              kfl_ortho_rad =  -2
           else if( words(2) == 'MARGA' ) then
              kfl_ortho_rad =  4
           end if

        else if(words(1)=='BUBBL') then
           if(words(2)=='ON   ') kfl_bubbl_rad=1

        else if(words(1)=='ASSEM') then
           if(words(2)=='ADR  ') then
              kfl_assem_rad=1
           else if(words(2)=='CELL ') then
              kfl_assem_rad=2
           end if

        else if(words(1)=='ELEME') then
           call realen(kfl_ellen_rad)
           if(words(2)=='NEW  ') kfl_ellen_rad=-1

        else if(words(1)=='NORMO') then  !!F ??
           if(exists('L1   ')) then
              kfl_normc_rad = 1
           else if(exists('L2   ')) then
              kfl_normc_rad = 2
           else if(exists('L-inf')) then
              kfl_normc_rad = 0
           else if(exists('ALGEB')) then
              kfl_normc_rad = 3
           end if

        else if(words(1)=='MAXIM') then  ! Maximum number of internal iterations
           miinn_rad = int(param(1))

        else if(words(1)=='CONVE') then  ! Tolerance for internal iterations
           cotol_rad = param(1)

        else if(words(1)=='ALGEB') then !!F ???
           call reasol(1_ip)

        else if(words(1)=='PRECO') then !!F ???
           call reasol(2_ip)

        end if
     end do
  end if

end subroutine rad_reanut

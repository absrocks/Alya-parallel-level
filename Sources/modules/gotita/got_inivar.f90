subroutine got_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_inivar
  ! NAME 
  !    got_inicar
  ! DESCRIPTION
  !    This routine initializes some variables
  !    ITASK=1 ... When starting the run (from Turnon)
  !    ITASK=2 ... First time step. This is needed as some variables 
  !                are not initialized before
  !    ITASK=3 ... When starting a time step (from got_begste)
  ! USES
  ! USED BY
  !    got_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_gotita
  use def_domain
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat

  select case(itask)

  case(0)   
     !
     ! Postprocess Variable names and types
     !
     postp(1) % wopos (1, 1) = 'CDROP'
     postp(1) % wopos (1, 2) = 'VDROP'
     postp(1) % wopos (1, 3) = 'MDROP'
     postp(1) % wopos (1, 4) = 'CEXAC'
     postp(1) % wopos (1, 5) = 'VEXAC'
     postp(1) % wopos (1, 6) = 'DRAG '
     postp(1) % wopos (1, 7) = 'AIRVE'
     postp(1) % wopos (1, 8) = 'TERM '
     postp(1) % wopos (1, 9) = 'CRESI'
     postp(1) % wopos (1,10) = 'VRESI'
     postp(1) % wopos (1,11) = 'FIXIT'
     postp(1) % wopos (1,12) = 'TIMES'
     postp(1) % wopos (1,13) = 'DIFFU'
     postp(1) % wopos (1,14) = 'BETA '

     postp(1) % wopos (2, 1) = 'SCALA'
     postp(1) % wopos (2, 2) = 'VECTO'
     postp(1) % wopos (2, 3) = 'VECTO'
     postp(1) % wopos (2, 4) = 'SCALA'
     postp(1) % wopos (2, 5) = 'VECTO'
     postp(1) % wopos (2, 6) = 'SCALA'
     postp(1) % wopos (2, 7) = 'VECTO'
     postp(1) % wopos (2, 8) = 'SCALA'
     postp(1) % wopos (2, 9) = 'SCALA'
     postp(1) % wopos (2,10) = 'SCALA'
     postp(1) % wopos (2,11) = 'SCALA'
     postp(1) % wopos (2,12) = 'SCALA'
     postp(1) % wopos (2,13) = 'SCALA'
     postp(1) % wopos (2,14) = 'VECTO'
     !
     ! Solver
     !     
     call soldef(-3_ip)
     solve(1) % wprob     = 'P1'        ! Equation name
     solve(1) % kfl_solve = 1           ! Output flag
     solve(2) % wprob     = 'P2'        ! Equation name
     solve(2) % kfl_solve = 1           ! Output flag
     solve(3) % wprob     = 'P3'        ! Equation name
     solve(3) % kfl_solve = 1           ! Output flag

  case(1)
     !
     ! Dimensions
     !
     if(kfl_probl_got==1) then
        ndofn_got(3) = ndime+1
        nevat_got    = ndofn_got(3)*mnode
        ndof2_got(3) = ndofn_got(3)*ndofn_got(3)
        if(kfl_algor_got==4) then
           !
           ! Block Gauss-Seidel scheme
           !
           ndofn_got(1) = ndime
           ndofn_got(2) = 1
           ndof2_got(1) = ndofn_got(1)*ndofn_got(1)
           ndof2_got(2) = ndofn_got(2)*ndofn_got(2)
           solve(1)%wprob = 'MOMENTUM'
           solve(2)%wprob = 'CONTINUITY'

        else
           !
           ! Monolithic scheme
           !
           ndofn_got(1) = ndime+1
           ndof2_got(1) = ndofn_got(1)*ndofn_got(1)
           solve(1)%wprob = 'MOMENTUM_CONTINUITY'

        end if
     else if(kfl_probl_got==2) then
        ndofn_got(1) = ndime
        nevat_got    = ndofn_got(1)*mnode
        ndof2_got(1) = ndofn_got(1)*ndofn_got(1)
     else if(kfl_probl_got==3) then
        ndofn_got(1) = 1
        nevat_got    = ndofn_got(1)*mnode
        ndof2_got(1) = ndofn_got(1)*ndofn_got(1)
     end if
     !
     ! Solver
     !
     solve(1) % ndofn = ndofn_got(1)
     solve(2) % ndofn = ndofn_got(2)
     solve(3) % ndofn = ndofn_got(3)
     solve(1) % ndof2 = ndof2_got(1)
     solve(2) % ndof2 = ndof2_got(2)
     solve(3) % ndof2 = ndof2_got(3)

     kfact_got=densi_got*ddrop_got*ddrop_got*veair_got&
          /(18.0_rp*leinf_got*muair_got)

  case(2)   
     !
     ! Number of droplet velocity components 
     !
     if(kfl_timei_got==1) then
        ncomp_got=3
     else
        ncomp_got=2     
     end if
     !
     ! Time variables
     !
     if(kfl_timei_got==0) then
        dtinv_got=1.0_rp
     else
        kfl_timei=1
     end if
     kfl_stead_got=0 
     !
     ! Time accuracy: save original value
     !
     kfl_tiaor_got=kfl_tiacc_got
     !
     ! Subgrid scale number of iterations
     !
     if(kfl_sgsco_got==0) misgs_got=1
     !
     ! Solver
     !
     itsol_got = 0
     !
     ! Artificial viscosity
     !
     if(kfl_artif_got==0) then
        artif_got(1)=0.0_rp
        artif_got(2)=0.0_rp
     end if
     !
     ! Penalization
     !
     if(kfl_penal_got==2) then
        penal_got=almax_got
     else if(kfl_penal_got==0) then 
        penal_got=0.0_rp
     end if

  end select

end subroutine got_inivar

subroutine chm_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* partis/chm_inivar
  ! NAME 
  !    chm_inivar
  ! DESCRIPTION
  !    This routine initializes some variables
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_chemic
  use def_solver
  use def_kintyp

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,iclas,kpoin,iboun,inodb,ispec
  integer(ip), pointer    :: limpo(:)
  integer(ip), save       :: nbou_set_vars
  character(2)            :: ww

  select case(itask)

  case(0_ip)
     !
     ! Postprocess
     !
     postp(1) % wopos (1,1) = 'CONCE'
     postp(1) % wopos (1,2) = 'CLUST'
     postp(1) % wopos (1,3) = 'TOTAL'
     postp(1) % wopos (1,4) = 'PROJE'
     postp(1) % wopos (1,5) = 'VELOC'
     postp(1) % wopos (1,6) = 'COSGS'
     postp(1) % wopos (1,7) = 'LIMIT'
     postp(1) % wopos (1,8) = 'COPRO' 
     postp(1) % wopos (1,9) = 'TOTCO'   ! Total concentration * density
     postp(1) % wopos (1,10) = 'ACCUM'  ! accumulation
     postp(1) % wopos (1,11) = 'SOURC'  ! source term
     postp(1) % wopos (1,12) = 'PROAD'  ! P
     postp(1) % wopos (1,13) = 'VISCO'  ! Viscosity
     postp(1) % wopos (1,14) = 'SPHEA'  ! Specific heat
     postp(1) % wopos (1,15) = 'CONDU'  ! heat conductivity
     postp(1) % wopos (1,16) = 'ENTHA'  ! Enthalpy
     postp(1) % wopos (1,17) = 'HEATS'  ! Enthalpy transport source term
     postp(1) % wopos (1,18) = 'CHEMI'  ! Chemical only Heat source term
     postp(1) % wopos (1,19) = 'SUMCO'  ! Sum of concentration
     postp(1) % wopos (1,20) = 'EQUIV'  ! Equivalence ratio of fuel
     postp(1) % wopos (1,21) = 'MOLEC'  ! Molecular weight
     postp(1) % wopos (1,22) = 'SENSO'  ! Flame front sensor for TFLES
     postp(1) % wopos (1,23) = 'SPEED'  ! Laminar flame speed for TFLES
     postp(1) % wopos (1,24) = 'THICK'  ! Laminar flame thickness for TFLES
     postp(1) % wopos (1,25) = 'WRINK'  ! Subgrid scale wrinkling factor for TFLES
     postp(1) % wopos (1,26) = 'FACTO'  ! Dynamic thickening factor F for DTFLES
     postp(1) % wopos (1,27) = 'TEMPE'  ! Temperature for low-mach CFI combustion model
     postp(1) % wopos (1,28) = 'AVTEM'  ! Averaged temperature
     postp(1) % wopos (1,29) = 'IMEAN'  ! Enthalpy scalar CFI model
     postp(1) % wopos (1,30) = 'AVCON'  ! Averaged concentration 
     postp(1) % wopos (1,31) = 'AVCO2'  ! Averaged squared of concentration
     postp(1) % wopos (1,32) = 'AVVAR'  ! Averaged VRPV
     postp(1) % wopos (1,33) = 'AVIME'  ! Averaged enthalpy scalar
     postp(1) % wopos (1,34) = 'AVCHM'  ! Averaged chemical heat
     postp(1) % wopos (1,35) = 'AVMIX'  ! Averaged mixture fraction
     postp(1) % wopos (1,36) = 'AVMI2'  ! Averaged squared of mixture fraction
     postp(1) % wopos (1,37) = 'SPEC1'  ! Species post-processing CFI model, radiation
     postp(1) % wopos (1,38) = 'SPEC2'  ! Species post-processing CFI model, radiation
     postp(1) % wopos (1,39) = 'RADIA'  ! Radiation source term CFI model
     postp(1) % wopos (1,40) = 'YSCAL'  ! Scaled reaction progress for partially premixed conditions
     postp(1) % wopos (1,41) = 'CVARI'  ! Variance of the scaled reaction progress varibale 

     postp(1) % wopos (2,1) = 'SCALA'
     postp(1) % wopos (2,2) = 'SCALA'
     postp(1) % wopos (2,3) = 'SCALA'
     postp(1) % wopos (2,4) = 'SCALA'
     postp(1) % wopos (2,5) = 'VECTO'
     postp(1) % wopos (2,6) = 'SCALA'
     postp(1) % wopos (2,7) = 'SCALA'
     postp(1) % wopos (2,8) = 'SCALA'
     postp(1) % wopos (2,9) = 'SCALA'
     postp(1) % wopos (2,10) = 'SCALA'
     postp(1) % wopos (2,11) = 'SCALA'
     postp(1) % wopos (2,12) = 'SCALA'
     postp(1) % wopos (2,13) = 'SCALA'
     postp(1) % wopos (2,14) = 'SCALA'
     postp(1) % wopos (2,15) = 'SCALA'
     postp(1) % wopos (2,16) = 'SCALA'
     postp(1) % wopos (2,17) = 'SCALA'
     postp(1) % wopos (2,18) = 'SCALA'
     postp(1) % wopos (2,19) = 'SCALA'
     postp(1) % wopos (2,20) = 'SCALA'
     postp(1) % wopos (2,21) = 'SCALA'
     postp(1) % wopos (2,22) = 'SCALA'
     postp(1) % wopos (2,23) = 'SCALA'
     postp(1) % wopos (2,24) = 'SCALA'     
     postp(1) % wopos (2,25) = 'SCALA'
     postp(1) % wopos (2,26) = 'SCALA'
     postp(1) % wopos (2,27) = 'SCALA'
     postp(1) % wopos (2,28) = 'SCALA'
     postp(1) % wopos (2,29) = 'SCALA'
     postp(1) % wopos (2,30) = 'SCALA'
     postp(1) % wopos (2,31) = 'SCALA'
     postp(1) % wopos (2,32) = 'SCALA'
     postp(1) % wopos (2,33) = 'SCALA'
     postp(1) % wopos (2,34) = 'SCALA'
     postp(1) % wopos (2,35) = 'SCALA'
     postp(1) % wopos (2,36) = 'SCALA'
     postp(1) % wopos (2,37) = 'SCALA'
     postp(1) % wopos (2,38) = 'SCALA'
     postp(1) % wopos (2,39) = 'SCALA'
     postp(1) % wopos (2,40) = 'SCALA'
     postp(1) % wopos (2,41) = 'SCALA'
     !
     ! Set variables 
     !
     postp(1) % wonse (1)     = 'CONCE'  
     postp(1) % woese (1)     = 'IPDE1' ! int_W M1 dw 
     postp(1) % woese (2)     = 'IODE1' ! int_W Q1 dw 
     postp(1) % woese (3)     = 'IODEN' ! int_W Qn dw 
     postp(1) % woese (4)     = 'TOTAL' ! int_W (q1+q2+...) dw 
     postp(1) % wobse (1)     = 'MASS ' ! int_S rho*u.n ds
     postp(1) % wobse (2:9)   = 'CONCE' ! <Yk> = int_S rho*u* Y_k dS / int_S rho*u dS, for k = 1,...,8
     nbou_set_vars = 1                ! Need to set this for later reentry when we know nspec
     !
     ! Witness variables
     !
     postp(1) % wowit (1:8) = 'CONCE' ! Species 1 to 8 (update the maximum number of species here)     
     !
     ! Solver
     !     
     call soldef(-1_ip)
     solve(1) % wprob       = 'CONCENTRATION'
     solve(1) % kfl_solve   = 1

     !
     ! Others
     !
     cputi_chm = 0.0_rp  ! CPU times
     dtmat_chm = 0.0_rp  ! Matrix-based time step

     !
     ! Nullify pointers
     !
     nullify(avtem_chm)
     nullify(avcon_chm)
     nullify(avco2_chm)
     nullify(avvar_chm)
     nullify(avime_chm)
     nullify(avchm_chm)
     nullify(avmix_chm)
     nullify(avmi2_chm)
     nullify(yscale_chm)
     nullify(cvar_chm)

  case(1_ip)

     if( nbou_set_vars+nclas_chm > nvars ) call runend('CHM_INIVAR')
     kfl_rsta2_chm= .false.                         ! restarting bdf file?

  case(3_ip)

     if(kfl_timei_chm==0) then                      ! Time integration
        dtinv_chm=1.0_rp
     else
        kfl_timei=1
     end if
     momod(modul) % kfl_stead = 0
     kfl_grdif_chm = 0
     nspec_chm     = nclas_chm+nodes_chm
     nspec = nspec_chm
     if( ISLAVE ) then
       do ispec = 1,nspec
        speci(ispec)%name = ''
       enddo
     endif
     kfl_tiaor_chm = kfl_tiacc_chm                  ! Time accuracy: save original value
     if( kfl_timei_chm == 1 ) then                  ! Number of velocity components
        if( kfl_tisch_chm == 1 ) then
           ncomp_chm = 3                            ! Trapezoidal rule
        else if( kfl_tisch_chm == 2 ) then
           ncomp_chm = 2+kfl_tiacc_chm              ! BDF scheme
        end if
        if( kfl_dttyp_chm == 2 ) then
           ncomp_chm = max(5_ip,ncomp_chm)
        else if( kfl_dttyp_chm == 3 ) then
           ncomp_chm = max(4_ip,ncomp_chm)
        end if
     else
        ncomp_chm = 2     
     end if
     ittot_chm = 0                                  ! Others
     iclai_chm = 1                                  ! Jacobi starting class
     iclaf_chm = nclas_chm                          ! Jacobi final class
     kfl_dttar_chm = min(kfl_dttar_chm,nspec_chm)   ! Time step target
     !
     ! METEO model variables
     !
     tmete_chm = timei
     tsour_chm = timei

  end select

end subroutine chm_inivar

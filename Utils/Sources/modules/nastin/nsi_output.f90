subroutine nsi_output()
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_output
  ! NAME 
  !    nsi_output
  ! DESCRIPTION
  !    End of a NASTIN time step 
  !    ITASK=0 ... When timemarching is true. There is output or
  !                post-process of results if required.
  !    ITASK=1 ... When timemarching is false. Output and/or post-process
  !                of results is forced if they have not been written 
  !                previously
  ! USES
  !    nsi_outvar
  !    postpr
  !    nsi_outbcs
  !    nsi_outset
  !    nsi_exaerr
  ! USED BY
  !    nsi_endste (itask=1)
  !    nsi_turnof (itask=2)
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_postpr
  use def_kermod,           only: ndivi,nsteps_ensemble
  use mod_output_postprocess

  implicit none
  integer(ip) :: ivari,ivarp,ii,devis(10)
  integer(ip) :: kfl_devis
  real(rp)    :: time1,time2
  logical(lg),    save       :: laux_ensemble=.true.

  
#ifdef EVENT
  call mpitrace_user_function(1)
#endif
  !
  ! Check if density and/or viscosity are needed
  !  
  kfl_devis = 0
  devis     = 0
  devis(1)  = 7  ! rho
  devis(2)  = 11 ! Tw
  devis(3)  = 17 ! mu
  devis(4)  = 25 ! y+
  devis(5)  = 30 ! Pe
  ii        = 1

  do while( ii < 10 .and. devis(ii) /= 0 )
     ivari = devis(ii)
     call posdef(11_ip,ivari)
     if( ivari /= 0 ) kfl_devis = 1
     ii = ii + 1
  end do

  if( kfl_devis == 1 ) then
     call nsi_memphy(4_ip)
     call nsi_denvis()
  end if

  !
  ! Postprocess variables
  !
! call cputim(time1)
!  do ivari = 1,nvarp
!     if( output_postprocess_check_variable_postprocess_now(ivari) ) &
!          call nsi_outvar(ivari)
!  end do
  do ivarp = 1,nvarp
     ivari = ivarp
     call posdef(11_ip,ivari)
     call nsi_outvar(ivari)
  end do! call cputim(time2)
! print*,'NSI_OUTPUT:',time2-time1

  if( postp(1) % npp_stepi(21) /= 0 ) then
     if ( nsteps_ensemble >= postp(1) % npp_stepi(21) ) call runend(' nsteps_ensemble >= postp(1) % npp_stepi(21) does not make sense')
     if( ( mod(ittim, postp(1) % npp_stepi(21) ) == nsteps_ensemble ) .and. laux_ensemble ) then
        entim_nsi = cutim    ! ensemble time
        laux_ensemble=.false.
     end if
     if( mod(ittim, postp(1) % npp_stepi(21) ) == 0 ) then     ! AVVEL frequency
        avtim_nsi = cutim  ! Update reference time for time-averaging
        laux_ensemble=.true.
     endif
  endif

  if( kfl_devis == 1 ) then
     call nsi_memphy(-4_ip)
  end if

  if( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then
     !
     ! Calculations on sets
     !
     call nsi_outset()
     !
     ! Postprocess on witness points
     !
     call nsi_outwit()
     !
     ! End of the run
     !
     call nsi_exaerr(2_ip)

  else if( ittyp == ITASK_ENDRUN ) then

     !if( INOTMASTER ) then
     !   call pspltm(&
     !        npoin,npoin,1_ip,0_ip,c_dom,r_dom,amatr(poapp_nsi),&
     !        trim(title)//': '//namod(modul),0_ip,18.0_rp,'cm',&
     !        0_ip,0_ip,2_ip,90+kfl_paral)
     !end if

  end if

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine nsi_output

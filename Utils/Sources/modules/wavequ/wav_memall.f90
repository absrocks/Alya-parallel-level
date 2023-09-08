subroutine wav_memall()
  !-----------------------------------------------------------------------
  !****f* wavequ/wav_memall
  ! NAME 
  !    wav_memall
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed to solve the
  !    wave equation equation
  ! USES
  ! USED BY
  !    wav_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_solver
  use def_wavequ
  use mod_memchk
  implicit none
  integer(ip) :: lodof,ielem,pelty,pgaus
  integer(4)  :: istat

  if(kfl_paral/=0) then
     !
     ! wave equation unknown: amplitude WAVAM
     ! 
     allocate(wavam(npoin,ncomp_wav),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'WAVAM','wav_memall',wavam)
     if(kfl_subgs_wav==1) then
        allocate(wasgs(nelem),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'WASGS','wav_memall',wasgs)
        do ielem=1,nelem
           pelty=ltype(ielem)
           pgaus=ngaus(pelty)
           allocate(wasgs(ielem)%a(mgaus,3),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'WASGS','wav_memall',wasgs(ielem)%a)
        end do
     end if
     if(kfl_timet_wav==1.and.kfl_absor_wav==1) then
        allocate(vmass_wav(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'VMASS_WAV','wav_memall',vmass_wav)        
     end if
     if(kfl_tisch_wav==3) then
        !
        ! Newmark: wave velocity and acceleration
        !
        allocate(wavac_wav(npoin,2),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'WAVAC','wav_memall',wavac_wav)        
        allocate(wavve_wav(npoin,2),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'WAVVE','wav_memall',wavve_wav)     
     end if
     !
     ! Solver
     !
     solve_sol => solve(1:)
     call soldef(4_ip)

  else

     allocate(wavam(1,3),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'WAVAM','wav_memall',wavam)

  end if

end subroutine wav_memall


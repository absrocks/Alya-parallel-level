subroutine wav_memphy()
  !------------------------------------------------------------------------
  !****f* Wavequ/wav_memphy
  ! NAME 
  !    wav_memphy
  ! DESCRIPTION
  !    Allocate memory for physical problem
  ! USES
  ! USED BY
  !    wav_sendat
  !    wav_reaphy
  !------------------------------------------------------------------------
  use def_master
  use def_domain
  use def_wavequ
  use mod_memchk
  implicit none
  integer(4) :: istat

  allocate(lmate_wav(nelem),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'LMATE_TEM','wav_reaphy',lmate_wav)
  lmate_wav=1
  allocate(densi_wav(ncoef_wav,nmate_wav),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'DENSI_TEM','wav_reaphy',densi_wav)
  allocate(kappa_wav(ncoef_wav,nmate_wav),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'KAPPA_TEM','wav_reaphy',kappa_wav)
  allocate(lawde_wav(nmate_wav),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'LAWDE_TEM','wav_reaphy',lawde_wav)
  allocate(lawka_wav(nmate_wav),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'LAWKA_TEM','wav_reaphy',lawka_wav)

end subroutine wav_memphy

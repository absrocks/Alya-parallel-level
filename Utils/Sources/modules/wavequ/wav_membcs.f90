subroutine wav_membcs(itask)
  !------------------------------------------------------------------------
  !****f* Wavequ/wav_membcs
  ! NAME 
  !    wav_membcs
  ! DESCRIPTION
  !    Allocate memory for boundary conditions
  ! USES
  ! USED BY
  !    wav_sendat
  !    wav_reabcs
  !------------------------------------------------------------------------
  use def_master
  use def_domain
  use def_wavequ
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat

  select case(itask)
     
  case(1)
     !
     ! Boundary conditions on nodes
     !
     call runend('WAV_MEMBCS: HAY QUE CAMBIAR FIXNO PARA QUE TENGA UNA PRIMERA DIMENSION')
     allocate(kfl_fixno_wav(npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXNO_WAV','wav_membcs',kfl_fixno_wav)
     allocate(bvess_wav(npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_WAV','wav_membcs',bvess_wav)

  case(2)
     !
     ! Boundary conditions on boundaries
     !
     if(kfl_onbou_wav/=1) then
        allocate(kfl_fixbo_wav(nboun),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXBO_WAV','wav_membcs',kfl_fixbo_wav)
     end if

  end select

end subroutine wav_membcs

subroutine wav_outvar(ivari)
  !------------------------------------------------------------------------
  !****f* Wavequ/wav_output
  ! NAME 
  !    wav_output
  ! DESCRIPTION
  !    Output a postprocess variable
  ! USES
  !    postpr
  !    memgen
  ! USED BY
  !    wav_output
  !***
  !------------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_wavequ
  use      mod_postpr
  use      mod_memchk
  implicit none
  integer(ip), intent(in) :: ivari
  integer(ip)             :: ipoin
  real(rp)                :: dummr,rutim
  !
  ! Define postprocess variable
  !
  rutim = cutim

  select case (ivari)  

  case(0)

     return

  case(1)
     !
     ! Wave amplitude
     !
     gesca => wavam(:,1) 

  case(2)
     !
     ! Wave velocity
     !
     gesca => wavve_wav(:,1) 

  case(3)
     !
     ! Wave acceleration
     !
     gesca => wavac_wav(:,1)

  end select

  call outvar(&
       ivari,&
       ittim,rutim,postp(1) % wopos(1,ivari))

end subroutine wav_outvar

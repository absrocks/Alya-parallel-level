subroutine wav_parall(itask)
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_parall
  ! NAME
  !    wav_parall
  ! DESCRIPTION
  !    This routine is a bridge to Parall service  
  ! USED BY
  !    wav_turnon
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_wavequ
  implicit none
  integer(ip), intent(in) :: itask

  if(kfl_paral>=0) then

     select case(itask)

     case(1)
        !
        ! Exchange data read in wav_reaphy, wav_reanut and wav_reaous
        ! always using MPI, even if this is a partition restart
        !
        call wav_sendat(1_ip)

     case(2)
        !
        ! Exchange data read in wav_reabcs
        !
        call wav_sendat(2_ip)

     case(3) 
        !
        ! Sum up residual contribution of slave neighbors
        !
        if(kfl_paral>0) then
           party =  3                     ! nodes
           pardi =  1                     ! 1-dimensions
           parki =  5                     ! Real number (2 dimensions in one)
           pard1 =  1
           parr1 => unkno(1:npoin)
           call Parall(400_ip)
           parr1 => null()
        end if

     case(4)
        !
        ! Sum up residual contribution of slave neighbors
        !
        if(kfl_paral>0) then
           party =  3                     ! nodes
           pardi =  1                     ! 1-dimensions
           parki =  5                     ! Real number (2 dimensions in one)
           pard1 =  1
           parr1 => rhsid(1:npoin)
           call Parall(400_ip)
           parr1 => null()
        end if

     case(5)
        !
        ! Modified mass matrix for absorbing b.c.'s
        !
        if(kfl_paral>0) then
           party =  3                     ! nodes
           pardi =  1                     ! 1-dimensions
           parki =  5                     ! Real number (2 dimensions in one)
           pard1 =  1
           parr1 => vmass_wav(1:npoin)
           call Parall(400_ip)
           parr1 => null()
        end if

     end select

  end if

end subroutine wav_parall

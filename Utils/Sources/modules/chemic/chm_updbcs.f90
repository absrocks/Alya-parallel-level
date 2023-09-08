subroutine chm_updbcs(itask)
  !-----------------------------------------------------------------------
  !****f* partis/chm_updbcs
  ! NAME 
  !    chm_updbcs
  ! DESCRIPTION
  !    This routine prepares a new time step
  ! USES
  !    chm_updunk
  ! USED BY
  !    partis
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  implicit none
  integer(ip), intent(in) :: itask
  real(rp)                :: udotn
  integer(ip)             :: ibopo,ipoin,idime,iffix,iclas,kfl_inout,inorm
  integer(ip), save       :: ipass=0
  real(rp)                :: a,b,c
  real(rp),    pointer    :: velol(:,:)

  if( IMASTER ) return

  select case ( itask )

  case( 1_ip )

     !-------------------------------------------------------------------
     !  
     ! Before a time step
     !     
     !-------------------------------------------------------------------

        if( ipass == 0 ) then
           ipass = 1
           do ipoin = 1,npoin
              do iclas = 1,nclas_chm

                 if(  kfl_fixno_chm(iclas,ipoin) == 2 .or.&
                      kfl_fixno_chm(iclas,ipoin) == 3 ) then
                    ibopo = lpoty(ipoin)
                    if( ibopo == 0 ) kfl_fixno_chm(iclas,ipoin) = 0
                 end if

              end do
           end do
        end if

        do ipoin = 1,npoin
           if( lpoty(ipoin) > 0 ) then
              do iclas = 1,nclas_chm

                 if(  abs(kfl_fixno_chm(iclas,ipoin)) == 2 .or. &
                      abs(kfl_fixno_chm(iclas,ipoin)) == 3 ) then

                    iffix = abs(kfl_fixno_chm(iclas,ipoin))
                    ibopo = lpoty(ipoin)
                    udotn = 0.0_rp
                    do idime = 1,ndime
                       udotn = udotn + veloc_chm(idime,ipoin)*exnor(idime,1,ibopo)
                    end do

                    if( lawvt_chm /= 0 ) then
                       udotn = udotn - vterm_chm(ipoin,iclas)*exnor(ndime,1,ibopo)
                    end if

                    if( udotn > 0.0_rp ) then
                       !
                       ! Outflow u.n > 0
                       !
                       kfl_fixno_chm(iclas,ipoin) = -iffix              
                    else
                       !
                       ! Inflow u.n <= 0
                       !
                       kfl_fixno_chm(iclas,ipoin) =  iffix
                       bvess_chm(iclas,ipoin)     =  0.0_rp
                    end if

                 end if
              end do
           end if
        end do

  case( 2_ip )

     !-------------------------------------------------------------------
     !
     ! Before a global iteration
     !  
     !-------------------------------------------------------------------
     !
     ! a*T^2 + b*T + c : Saturation fraction of water in air
     !
     a = 2.0e-5_rp
     b = 0.0003_rp
     c = 0.0025_rp     
     inorm = 0
     do ipoin = 1,npoin
        do iclas = 1,nclas_chm
           if ( kfl_fixno_chm(iclas,ipoin) == 2 ) then
              inorm = 1
              bvess_chm(iclas,ipoin) = a * tempe(ipoin,1)*tempe(ipoin,1) + b * tempe(ipoin,1) + c
              conce(ipoin,iclas,1)   = bvess_chm(iclas,ipoin)
           end if
        end do
     end do
     if (kfl_norma_chm > 0 .and. inorm /= 0 ) call chm_updunk(9_ip) ! Normalize

  end select

end subroutine chm_updbcs

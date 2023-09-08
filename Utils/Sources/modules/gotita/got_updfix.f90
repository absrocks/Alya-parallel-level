subroutine got_updfix()
  !-----------------------------------------------------------------------
  !****f* Gotita/got_updfix
  ! NAME 
  !    got_updfix
  ! DESCRIPTION
  !    This routine updates the boundary conditions of alpha
  ! USES
  ! USED BY
  !    got_solite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_gotita
  implicit none
  integer(ip) :: ipoin,idime,ibopo
  real(rp)    :: udotn

  if(kfl_paral/=0) then
     if(kfl_exacs_got==0.and.kfl_probl_got/=2) then
        do ipoin=1,npoin
           if(  .not.(coord(1,ipoin)<xmini_got(1).or.coord(1,ipoin)>xmaxi_got(1).or.&
                &     coord(2,ipoin)<xmini_got(2).or.coord(2,ipoin)>xmaxi_got(2).or.&
                &     coord(ndime,ipoin)<xmini_got(ndime).or.coord(ndime,ipoin)>xmaxi_got(ndime))) then
              ibopo=lpoty(ipoin)
              if(ibopo>0.and.kfl_fixno_got(ipoin)/=1) then
                 udotn=0.0_rp
                 do idime=1,ndime
                    udotn=udotn+vdrop(idime,ipoin,1)*exnor(idime,1,ibopo)
                 end do
                 if(udotn<=0.0_rp) then
                    kfl_fixno_got(ipoin)=2
                 else
                    kfl_fixno_got(ipoin)=0
                 end if
              end if
           end if
        end do
     end if
  end if

end subroutine got_updfix

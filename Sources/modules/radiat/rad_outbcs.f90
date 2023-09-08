subroutine rad_outbcs
!-----------------------------------------------------------------------
!****f* Radiat/rad_outbcs
! NAME 
!    rad_outbcs
! DESCRIPTION
!    Postprocess boundary conditions. This could be useful to include
!    this new file when running the same problem over again.
! USES
! USED BY
!    rad_output
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_radiat
  use      mod_iofile
  implicit none
  integer(ip) :: ipoin

  do ipoin=1,npoin
     write(lun_bound_rad,1,advance='no') ipoin
     write(lun_bound_rad,2,advance='no') kfl_fixno_rad(1,ipoin)        
     write(lun_bound_rad,5,advance='no') radav_rad(ipoin,1)
     if(kfl_conbc_rad==0) then
        write(lun_bound_rad,4,advance='no') kfl_funno_rad(ipoin)
     end if
  end do
 
  close(lun_bound_rad)

1 format(1x,i7)
2 format(1x,i1)
4 format(1x,i2)
5 format(1x,e16.8E3)

end subroutine rad_outbcs

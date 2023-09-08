subroutine ada_chkpoi(wopad)
!-----------------------------------------------------------------------
!****f* adapti/ada_chkpoi
! NAME 
!    ada_chkpoi
! DESCRIPTION
!    This routine writes the checkpoint-like file with the projected
!    values
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_inpout
  use      mod_postpr
  use      def_adapti
  use      mod_iofile
  implicit none
  integer(ip)     :: lunau,npoau,nelau,kflau
  character(5)    :: wopad(2)
  
  

  kflau     = kfl_outfo
  kfl_outfo = -2
  lunau     = lun_postp
  lun_postp = lun_chpoi_ada  
  npoau= npoin
  nelau= nelem
  npoin= npnew_ada
  nelem= nenew_ada
  
  if (wopad(1)=='START') then 
     write(lun_postp) ittim,cutim,dtime     
  else if (wopad(2) == 'SCALA') then
     call postpr(scnew_ada,wopad,ittim,cutim)
  else if (wopad(2) == 'VECTO') then
     call postpr(venew_ada,wopad,ittim,cutim)
  end if
  
  npoin     = npoau
  nelem     = nelau
  lun_postp = lunau
  kfl_outfo = kflau
  

end subroutine ada_chkpoi

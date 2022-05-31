subroutine ibm_begste()
  !-----------------------------------------------------------------------
  !****f* immbou/ibm_begste
  ! NAME
  !    ibm_begste
  ! DESCRIPTION
  !    This routines begins a time step
  ! USED BY
  !    Immbou
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_immbou
  implicit none
  integer(ip) :: iimbo,idime
  !
  ! Put force to zero
  ! 
  do iimbo = 1,nimbo
     do idime = 1,ndime
        imbou(iimbo)%force(idime,1) = 0.0_rp
     end do
  end do
  !
  ! Obtain xrota_ibm and xline_ibm from nstro_ibm and nstli_ibm
  !
  xrota_ibm = 0.0_rp
  xline_ibm = 0.0_rp
  do idime =1,3_ip
     if ( ittim >= nstro_ibm(idime) )  xrota_ibm(idime) = 1.0_rp
     if ( ittim >= nstli_ibm(idime) )  xline_ibm(idime) = 1.0_rp
  end do 

end subroutine ibm_begste

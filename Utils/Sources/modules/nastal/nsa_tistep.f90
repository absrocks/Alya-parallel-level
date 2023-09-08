subroutine nsa_tistep
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_tistep
  ! NAME 
  !    nsa_tistep
  ! DESCRIPTION
  !    This routine sets the time step already computed in nsa_updtss
  !
  ! USED BY
  !    nsa_begste
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_parame
  use      def_domain

  use      def_nastal
  use mod_outfor, only : outfor

  implicit none 

  real(rp)    :: dtaux(10),dtlof,dtmii,dtmix
  integer(ip) :: itimi,ipoin

  dtmii= dtcri_nsa(1)
  dtmix= dtcri_nsa(1)

  do itimi=1,ndtdf_nsa
     if (dtcri_nsa(itimi) < dtmii) dtmii= dtcri_nsa(itimi)
     if (dtcri_nsa(itimi) > dtmix) dtmix= dtcri_nsa(itimi)
  end do

  if(kfl_timco/=2) then
     dtaux = dtime
  else if(kfl_timco==2) then
     dtaux = dtcri_nsa
  end if
  
  !
  ! Compute BDF parameters
  !

  if (kfl_tisch_nsa == 1 .or. kfl_tisch_nsa == 3) then
     !
     ! Euler scheme (for or back) and CN scheme
     !
     call parbdf(kfl_tiacc_nsa,pabdf_nsa)

  else if(kfl_tisch_nsa==2) then
     !
     ! BDF higher-order scheme: increase integration order at each time step
     !
     if(ittim<=neule_nsa) then
        kfl_tiacc_nsa=1
     else
        kfl_tiacc_nsa=min(kfl_tiaor_nsa,ittim)
     end if
     call parbdf(kfl_tiacc_nsa,pabdf_nsa)
     
  end if

  nbdfp_nsa = kfl_tiacc_nsa + 1

  routp(1)=dtmii
  routp(2)=0.0_rp
  routp(3)=0.0_rp
  ioutp(1)=kfl_timei_nsa
  ioutp(2)=kfl_stead_nsa
  call outfor(8_ip,lun_outpu,' ')

end subroutine nsa_tistep

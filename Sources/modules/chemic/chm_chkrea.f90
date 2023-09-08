subroutine chm_chkrea()
  !------------------------------------------------------------------------
  !****f* partis/chm_chkrea
  ! NAME 
  !    chm_chkrea
  ! DESCRIPTION
  !    Defines the chemical reactions
  ! USES
  ! USED BY
  !    chm_reaphy
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_chemic
  use def_domain
  use mod_memchk
  implicit none
  integer(ip)   :: ireac,jreac,ncoef,iequa,icoef
  integer(ip)   :: lirea(100),ljrea(100)
  character(20) :: wirea,wjrea

  do ireac = 1,nreac_chm
     ncoef = size(lreac_chm(ireac)%l)
     do icoef = 1,ncoef
        lirea(icoef) = lreac_chm(ireac)%l(icoef)
     end do
     call heapsorti1(2_ip,ncoef,lirea)
     do jreac = ireac+1,nreac_chm
        iequa = 0
        if( ncoef == size(lreac_chm(jreac)%l) ) then
           do icoef = 1,ncoef
              ljrea(icoef) = lreac_chm(jreac)%l(icoef)
           end do
           call heapsorti1(2_ip,ncoef,ljrea)
           do icoef = 1,ncoef
              if( lirea(icoef) == ljrea(icoef) ) then
                 iequa = iequa + 1
              end if
           end do
           if( iequa == ncoef ) then
              wirea = intost(ireac)
              wjrea = intost(jreac)
              call runend('CHM_CHKREA: REACTIONS '//trim(wirea)//' AND '//trim(wjrea)//' ARE IDENTICAL')
           end if
        end if
     end do
  end do

end subroutine chm_chkrea

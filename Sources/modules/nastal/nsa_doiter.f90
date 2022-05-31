subroutine nsa_doiter
!-----------------------------------------------------------------------
!****f* Nastal/nsa_doiter
! NAME 
!    nsa_doiter
! DESCRIPTION
!    This routine solves an iteration of the NS equations
! USES
!    nsa_begite
!    nsa_solite
!    nsa_endite
! USED BY
!    Nastal
!***
!-----------------------------------------------------------------------
use def_parame
use def_master
use def_solver
use def_nastal

implicit none
   
  if(kfl_stead_nsa==0) then
     call nsa_begite
     do while(kfl_goite_nsa==1)
        call nsa_solite
        call nsa_endite(one)
     end do
     call nsa_endite(two)
  end if

end subroutine nsa_doiter


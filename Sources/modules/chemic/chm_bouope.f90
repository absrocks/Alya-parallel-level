subroutine chm_bouope()
  !-------------------------------------------------------------------------
  !****f* emic/chm_bouope
  ! NAME 
  !    chm_bouope
  ! DESCRIPTION
  !    Boundary operations
  ! USES
  ! USED BY
  !    chm_matrix 
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_ADR, only : ADR_add_sgs_or_bubble
  implicit none
  real(rp)    :: elmat(mnode,mnode)
  real(rp)    :: elrhs(mnode)
  real(rp)    :: gbdif(mgaub) 
  real(rp)    :: baloc(9),bocod(ndime,mnodb)
  integer(ip) :: ielem,iboun,inodb,pblty,kfl_robin
  integer(ip) :: pnodb,pnode,pelty,pgaus,pgaub
  integer(ip) :: iclas,izmat,izrhs,idime,ipoin
  real(rp)    :: gbsur(mgaub)
  real(rp)    :: dummr(ndime*mnode)

  kfl_robin = 0

  if( INOTMASTER .and. kfl_robin_chm /= 0 ) then
     iboun     = 0
     do while( iboun < nboun )
        iboun = iboun + 1
        iclas = iclai_chm - 1
        do while( iclas < iclaf_chm )
           iclas = iclas + 1
           if( kfl_fixbo_chm(iclas,iboun) == 2 ) then
              kfl_robin = 1
              iclas     = iclaf_chm
              iboun     = nboun
           end if
        end do
     end do

     if( kfl_robin == 1 ) then
       call runend('ROBIN-TYPE BCs ARE NOT CODED IN CHEMIC') 
     end if

  else  !! Master

     if( kfl_robin == 1 ) then
       call runend('ROBIN-TYPE BCs ARE NOT CODED IN CHEMIC') 
     end if

  end if

end subroutine chm_bouope

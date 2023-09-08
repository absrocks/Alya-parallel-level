subroutine dod_order_boundaries(&
     itask,ipoin,mnodb,nboun_local,npoin_global,nedge_ipoin,ltypb_global,&
     lnodb_global,lboch_global,ledbo_global,ledge_npoin_global,&
     list_ordered_boundaries,lperm,iview,kboun)
  use def_elmtyp
  use def_kintyp
  implicit none  
  integer(ip) , intent(in)    :: itask
  integer(ip) , intent(in)    :: ipoin
  integer(ip) , intent(in)    :: mnodb
  integer(ip) , intent(inout) :: nboun_local
  integer(ip) , intent(in)    :: npoin_global
  integer(ip) , intent(in)    :: nedge_ipoin
  integer(ip) , intent(in)    :: ltypb_global(*)
  integer(ip) , intent(in)    :: lnodb_global(mnodb,*)
  integer(ip) , intent(in)    :: lboch_global(*)
  integer(ip) , intent(inout)    :: ledbo_global(4,*)
  integer(ip) , intent(in)    :: ledge_npoin_global(nedge_ipoin)
  integer(ip) , intent(inout) :: list_ordered_boundaries(nboun_local)
  integer(ip) , intent(inout) :: lperm(3,4)
  integer(ip) , intent(out)   :: iview
  integer(ip) , intent(inout) :: kboun
  integer(ip)                 :: boun1,boun2
  integer(ip)                 :: kont,cont,flag,cpoin,nbopa,inodb
  integer(ip)                 :: iboun,jboun,ii,jj,ibopo,eboun,aux(4),salto
  integer(ip)                 :: iboun_last,iedge,kedge,inodb_1,inodb_2,ipoin_1,ipoin_2
  integer(ip)                 :: tab(2,4,2),itab
  integer(ip)                 :: nedge,pblty,nnode
  !logical(lg)                 :: ifoun
  integer(ip)                 :: ifoun
  !
  ! TRI03
  !
  tab(1,1,1) = 1
  tab(2,1,1) = 2
  tab(1,2,1) = 2
  tab(2,2,1) = 3
  tab(1,3,1) = 3
  tab(2,3,1) = 1
  !
  ! QUA04
  !
  tab(1,1,2) = 1
  tab(2,1,2) = 2
  tab(1,2,2) = 2
  tab(2,2,2) = 3
  tab(1,3,2) = 3
  tab(2,3,2) = 4
  tab(1,4,2) = 4
  tab(2,4,2) = 1

  if( itask == 1 ) then
     !
     ! Choose first edge and order first two boundaries
     !
     nedge   = nedge_ipoin
     iedge   = ledge_npoin_global(1)
     boun1   = ledbo_global(1,iedge)
     pblty   = abs(ltypb_global(boun1))
     if( pblty == TRI03 ) then
        itab = 1
     else
        itab = 2
     end if

     inodb_1 = ledbo_global(2,iedge)
     inodb_2 = ledbo_global(3,iedge)
     ipoin_1 = lnodb_global(inodb_1,boun1)
     ipoin_2 = lnodb_global(inodb_2,boun1)
     if( ipoin == ipoin_1 ) then
        if( inodb_2 /= tab(2,inodb_1,itab) ) then
           boun1 = ledbo_global(4,iedge)
           boun2 = ledbo_global(1,iedge)           
        else
           boun2 = ledbo_global(4,iedge)
        end if
     else
        if( inodb_1 /= tab(2,inodb_2,itab) ) then
           boun1 = ledbo_global(4,iedge)
           boun2 = ledbo_global(1,iedge)           
        else
           boun2 = ledbo_global(4,iedge)
        end if
     end if
     if( boun1 == 0 .or. boun2 == 0 ) then
        call runend('DOD_ORDENA: COULD NOT ORDER BOUNDARIES')
     end if
     !
     ! BOUN1 is first boundary
     ! BOUN2 is next  boundary
     !
     list_ordered_boundaries(1) =  boun1
     list_ordered_boundaries(2) =  boun2
     ledbo_global(1,iedge)      = -ledbo_global(1,iedge)
     ledbo_global(4,iedge)      = -ledbo_global(4,iedge)
     !
     ! Look for next boundaries
     !
     do iboun = 3,nboun_local
        iboun_last = list_ordered_boundaries(iboun-1)
        ifoun      = 0
        iedge_loop: do kedge = 2,nedge
           iedge = ledge_npoin_global(kedge)
           if( ledbo_global(1,iedge) == iboun_last ) then
              ifoun = 4
              exit iedge_loop
           else if( ledbo_global(4,iedge) == iboun_last ) then
              ifoun = 1
              exit iedge_loop
           end if
        end do iedge_loop
        if( ifoun /= 0 ) then         
           list_ordered_boundaries(iboun) =  ledbo_global(ifoun,iedge)
           ledbo_global(ifoun,iedge)      = -ledbo_global(ifoun,iedge)
        else
           call runend('DOD_ORDENA: COULD NOT FIND NEXT BOUDNARY')
        end if
     end do
     do kedge = 1,nedge
        iedge = ledge_npoin_global(kedge)
        ledbo_global(1,iedge) = abs(ledbo_global(1,iedge))
        ledbo_global(4,iedge) = abs(ledbo_global(4,iedge))
     end do
     !
     ! Take off BOFEM boundaries
     !
     nbopa = 0
     do iboun = 1,nboun_local
        eboun = list_ordered_boundaries(iboun)
        if( lboch_global(eboun) /= BOFEM ) then      
           nbopa = nbopa + 1
           list_ordered_boundaries(nbopa) = eboun
        end if
     end do
     nboun_local = nbopa


  else if( itask == 2 ) then

     !-------------------------------------------------------------------
     !
     ! Compute permutation
     !
     !-------------------------------------------------------------------

     iboun = list_ordered_boundaries(kboun)
     if( ltypb_global(iboun) == TRI03 ) then
       if( lnodb_global(1,iboun) == ipoin ) then
           iview          = 1           
           lperm(1,iview) = 2
           lperm(2,iview) = 3
        else if( lnodb_global(2,iboun) == ipoin ) then
           iview          = 2
           lperm(1,iview) = 3
           lperm(2,iview) = 1
        else
           iview          = 3
           lperm(1,iview) = 1
           lperm(2,iview) = 2
        end if

     else if( ltypb_global(iboun) == QUA04 ) then

        if(      lnodb_global(1,iboun) == ipoin ) then
           iview          = 1
           lperm(1,iview) = 2
           lperm(2,iview) = 3
           lperm(3,iview) = 4
        else if( lnodb_global(2,iboun) == ipoin ) then
           iview          = 2
           lperm(1,iview) = 3
           lperm(2,iview) = 4
           lperm(3,iview) = 1
        else if( lnodb_global(3,iboun) == ipoin ) then
           iview          = 3
           lperm(1,iview) = 4
           lperm(2,iview) = 1
           lperm(3,iview) = 2
        else
           iview          = 4
           lperm(1,iview) = 1
           lperm(2,iview) = 2
           lperm(3,iview) = 3
        end if

     end if
  end if

end subroutine dod_order_boundaries

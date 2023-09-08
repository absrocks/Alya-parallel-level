subroutine nsa_bounod(kfl_bours_nsa,kfl_boufu_nsa)
!-----------------------------------------------------------------------
!****f* Nastal/nsa_bounod
! NAME 
!    nsa_bounod
! DESCRIPTION
!    This routine passes conditions on boundaries to nodal conditions.
! USES
! USED BY
!    nsa_reabcs
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastal
  implicit none
  integer(ip), intent(in) :: kfl_bours_nsa(nboun),kfl_boufu_nsa(nboun)
  integer(ip)             :: iboun,pnodb,inodb,ipoin,idime,ibopo
!
! Loop over boundaries
!
  do iboun=1,nboun
     !
     ! Dirichlet 
     !
     if(kfl_fixbo_nsa(iboun)==21) then
        pnodb=nnode(ltypb(iboun))
        do inodb=1,pnodb
           ipoin=lnodb(inodb,iboun) 
           if(  kfl_fixno_nsa(    1,ipoin)==-1.and.&
                kfl_fixno_nsa(    2,ipoin)==-1.and.&
                kfl_fixno_nsa(ndime,ipoin)==-1) then
              ibopo=lpoty(ipoin)
              kfl_fixno_nsa(1:ndime,ipoin)=1
              do idime=1,ndime
                 bvess_nsa(idime,ipoin,1)=bvnat_nsa(idime,iboun)
              end do 
              kfl_fixrs_nsa(ibopo)=kfl_bours_nsa(iboun)
              if(kfl_conbc_nsa==0) then
                 kfl_funno_nsa(ipoin)=kfl_boufu_nsa(iboun)
                 do idime=1,ndime
                    bvess_nsa(idime,ipoin,2)=bvess_nsa(idime,ipoin,1)
                 end do                 
              end if
           end if
        end do
     !
     ! No slip wall
     !
     else if(kfl_fixbo_nsa(iboun)==27) then
        pnodb=nnode(ltypb(iboun))
        do inodb=1,pnodb
           ipoin=lnodb(inodb,iboun)
           if(  kfl_fixno_nsa(    1,ipoin)==-1.and.&
                kfl_fixno_nsa(    2,ipoin)==-1.and.&
                kfl_fixno_nsa(ndime,ipoin)==-1) then
              ibopo=lpoty(ipoin)
              kfl_fixno_nsa(1:ndime,ipoin)=1
              bvess_nsa(1:ndime,ipoin,1)=0.0_rp
              kfl_fixrs_nsa(ibopo)=0
              if(kfl_conbc_nsa==0) kfl_funno_nsa(ipoin)=0
           end if
        end do        
     end if
  end do

end subroutine nsa_bounod

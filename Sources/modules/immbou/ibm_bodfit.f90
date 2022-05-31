subroutine ibm_bodfit()
  !-----------------------------------------------------------------------
  !****f* ibm_bodfit/ibm_bodfit
  ! NAME
  !    ibm_bodfit
  ! DESCRIPTION
  !    This routines computes the body fitted contribution 
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_immbou
  use mod_kdtree
  use mod_memchk
  implicit none
  integer(ip)          :: pblty,iimbo
  integer(ip)          :: idime,isets,kpoin,kboun,inodb,iboun,ipoin

  if( IMASTER ) then
     !
     ! No boundaries nor nodes
     !
     do iimbo = 1,nimbo
        if( imbou(iimbo) % kfl_typeb >= 1 ) then
           imbou(iimbo) % nboib = 0
           imbou(iimbo) % npoib = 0
        end if
     end do

  else

     call memgen(1_ip,npoin,0_ip)

     do iimbo = 1,nimbo

        if( imbou(iimbo) % kfl_typeb >= 1 ) then

           igene = iimbo
           mnoib = max(mnoib,mnodb)
           mgaib = max(mgaib,mgaub)
           imbou(iimbo) % nboib = 0
           imbou(iimbo) % npoib = 0
           isets = imbou(iimbo) % kfl_typeb
           !
           ! Count number of boundaries and mark nodes
           !
           kboun = 0
           do iboun = 1,nboun
              if( lbset(iboun) == isets ) then
                 lexib(ltypb(iboun)) = 2
                 kboun = kboun + 1
                 do inodb = 1,nnode(ltypb(iboun))
                    ipoin = lnodb(inodb,iboun)
                    gisca(ipoin) = 1
                 end do
              end if
           end do
           imbou(iimbo) % nboib = kboun
           !
           ! Count number of nodes
           !
           kpoin = 0
           do ipoin = 1,npoin
              if( gisca(ipoin) == 1 ) then
                 kpoin = kpoin + 1
                 gisca(ipoin) = kpoin
              end if
           end do
           imbou(iimbo) % npoib = kpoin
           !
           ! Allocate memory
           !
           call ibm_memall(3_ip)
           !
           ! Coordinates
           !
           do ipoin = 1,npoin
              if( gisca(ipoin) /= 0 ) then
                 kpoin = gisca(ipoin)
                 imbou(iimbo) % lninv(kpoin) = ipoin
                 do idime = 1,ndime
                    imbou(iimbo) % cooib(idime,kpoin) = coord(idime,ipoin)
                 end do
              end if
           end do
           !
           ! Boundaries
           !
           kboun = 0
           do iboun = 1,nboun
              if( lbset(iboun) == isets ) then
                 kboun = kboun + 1
                 imbou(iimbo) % ltyib(kboun) = ltypb(iboun)
                 imbou(iimbo) % lbinv(kboun) = iboun
                 do inodb = 1,nnode(ltypb(iboun))
                    imbou(iimbo) % lnoib(inodb,kboun) = gisca(lnodb(inodb,iboun))
                 end do
              end if
           end do
           !
           ! Put GISCA to 0
           !
           do ipoin = 1,npoin
              gisca(ipoin) = 0
           end do
           !
           ! Gauss points and integration rule
           !
           do pblty = ibsta_dom,ibsto_dom
              if( lexib(pblty) == 2 ) then
                 ngaib(pblty) = ngaus(pblty)
                 lruib(pblty) = lrule(pblty)
              end if
           end do

        end if

     end do
     !
     ! Shape functions
     !
     call ibm_cshder(2_ip)
     do pblty = ibsta_dom,ibsto_dom
        if( lexib(pblty) == 2 ) lexib(pblty) = 1
     end do
     !
     ! GISCA: Deallocate
     !
     call memgen(3_ip,npoin,0_ip)

  end if

end subroutine ibm_bodfit

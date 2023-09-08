subroutine ibm_wallib()
  !-----------------------------------------------------------------------
  !****f* immbou/ibm_wallib
  ! NAME
  !    ibm_wallib
  ! DESCRIPTION
  !    This routines computes the host elements of the Gauss points for
  !    the IBM. In parallel, if several subdomains find a host element
  !    for the same Gauss point (due to finite tolerance), only the
  !    subdomain with the lowest number keep it.
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_immbou 
  use mod_memchk
  use mod_htable
  use mod_kdtree
  implicit none
  integer(ip)          :: inodb,ipoin,iboun,mbpoi,izbou,msize,pnodb
  integer(ip)          :: jboun,inode,kboun,mbobo,pblty,jnodb,iwaib
  integer(ip)          :: idime,kpoin,nboun_tmp,ii,dummi
  integer(4)           :: istat 
  integer(ip), pointer :: lbopo(:),pbopo(:),nbpoi(:)   ! Node to boundary graph
  integer(ip), pointer :: lbobo(:),pbobo(:)            ! Boundary to boundary graph
  integer(ip), pointer :: lista(:)
  integer(ip), pointer :: lconv(:)
  type(hash_t)         :: ht
  integer(ip), pointer :: lwaib(:,:)
  integer(ip), pointer :: lnodb_tmp(:,:)
  integer(ip), pointer :: lboel_tmp(:,:)
  integer(ip), pointer :: lelbo_tmp(:)
  integer(ip), pointer :: ltypb_tmp(:)
  integer(ip), pointer :: nwaib_tmp(:)
  real(rp),    pointer :: bouno_tmp(:,:)
  real(rp)             :: pladi,inorm(3),jnorm(3),xfact
  integer(ip)          :: ipoib,jpoib,sboun,sboun_tmp,nconv,mblty,iconv,conve

  if( INOTMASTER ) then

     !----------------------------------------------------------------------
     !
     ! Indentify walls and count wall boundaries: MWAIB
     !
     !----------------------------------------------------------------------

     allocate( lwaib(2,nboun) , stat = istat ) 
     call memchk(zero,istat,memor_dom,'LWAIB','ibm_wallib',lwaib)

     nboun_tmp = 0

     if(  kfl_geome == 1 ) then
        !
        ! Use geometrical boundary conditions
        !
        do iboun = 1,nboun
           if(       kfl_geobo(iboun) == 30 &
                .or. kfl_geobo(iboun) == 40 &
                .or. kfl_geobo(iboun) == 50 ) then
              nboun_tmp = nboun_tmp + 1
              lwaib(1,nboun_tmp) = iboun
           end if
        end do

     else if( 1 == 1 ) then 
        !
        ! All boundaries are walls
        !
        do iboun = 1,nboun
           nboun_tmp = nboun_tmp + 1
           lwaib(1,nboun_tmp) = iboun
        end do

     end if

     !----------------------------------------------------------------------
     !
     ! Allocate memory
     !
     !----------------------------------------------------------------------

     allocate( ltypb_tmp(nboun_tmp) , stat = istat ) 
     call memchk(zero,istat,memor_dom,'LTYPB_TMP','ibm_wallib',ltypb_tmp)
     allocate( lnodb_tmp(mnodb,nboun_tmp) , stat = istat ) 
     call memchk(zero,istat,memor_dom,'LNODB_TMP','ibm_wallib',lnodb_tmp)
     allocate( lboel_tmp(mnodb,nboun_tmp) , stat = istat ) 
     call memchk(zero,istat,memor_dom,'LBOEL_TMP','ibm_wallib',lboel_tmp)
     allocate( bouno_tmp(ndime,nboun_tmp) , stat = istat ) 
     call memchk(zero,istat,memor_dom,'BOUNO_TMP','ibm_wallib',bouno_tmp)
     allocate( lelbo_tmp(nboun_tmp) , stat = istat ) 
     call memchk(zero,istat,memor_dom,'LELBO_TMP','ibm_wallib',lelbo_tmp)

     do iwaib = 1,nboun_tmp
        iboun            = lwaib(1,iwaib)
        ltypb_tmp(iwaib) = ltypb(iboun)
        pnodb            = nnode(ltypb(iboun))
        do inodb = 1,pnodb
           lboel_tmp(inodb,iwaib) = lboel(inodb,iboun)
           lnodb_tmp(inodb,iwaib) = lnodb(inodb,iboun)
        end do
        lelbo_tmp(iwaib)         = lelbo(iboun)
     end do
     if( nboun_tmp > 0 ) then
        call bounor(nboun_tmp,lnodb_tmp,ltypb_tmp,lelbo_tmp,dummi,bouno_tmp)
     end if

     call memchk(two,istat,memor_dom,'LBOEL_TMP','ibm_wallib',lboel_tmp)
     deallocate(lboel_tmp,stat=istat)
     if(istat/=0) call memerr(two,'LBOEL_TMP','ibm_wallib',0_ip)

     call memchk(two,istat,memor_dom,'LELBO_TMP','ibm_wallib',lelbo_tmp)
     deallocate(lelbo_tmp,stat=istat)
     if(istat/=0) call memerr(two,'LELBO_TMP','ibm_wallib',0_ip)

     !----------------------------------------------------------------------
     !
     ! Compute boundary graph
     !
     !----------------------------------------------------------------------
     !
     ! Allocate memory for NBPOI and compute it
     !
     allocate(nbpoi(npoin),stat=istat)
     call memchk(zero,istat,memor_dom,'NBPOI','ibm_wallib',nbpoi)

     do iboun = 1,nboun_tmp
        do inode = 1,nnode(ltypb_tmp(iboun))
           ipoin = lnodb_tmp(inode,iboun)
           nbpoi(ipoin) = nbpoi(ipoin) + 1
        end do
     end do
     !
     ! Allocate memory for PBOPO and compute it
     !
     allocate(pbopo(npoin+1),stat=istat)
     call memchk(zero,istat,memor_dom,'PBOPO','ibm_wallib',pbopo)
     pbopo(1) = 1
     do ipoin = 1,npoin
        pbopo(ipoin+1) = pbopo(ipoin) + nbpoi(ipoin)
     end do
     !
     ! Allocate memory for LBOPO and construct the list
     !
     allocate(lbopo(pbopo(npoin+1)),stat=istat)
     call memchk(zero,istat,memor_dom,'LBOPO','ibm_wallib',lbopo)
     do iboun = 1,nboun_tmp
        do inode = 1,nnode(ltypb_tmp(iboun))
           ipoin = lnodb_tmp(inode,iboun)
           lbopo(pbopo(ipoin)) = iboun
           pbopo(ipoin) = pbopo(ipoin) + 1
        end do
     end do
     !
     ! Recompute PBOPO and maximum number of element neighbors MBPOI
     !
     pbopo(1) = 1
     mbpoi   = -1
     do ipoin = 1,npoin
        pbopo(ipoin+1) = pbopo(ipoin) + nbpoi(ipoin)
        mbpoi = max(mbpoi,nbpoi(ipoin))
     end do
     !
     ! Deallocate memory for temporary node/element connectivity
     !
     call memchk(two,istat,memor_dom,'NBPOI','ibm_wallib',nbpoi)
     deallocate(nbpoi,stat=istat)
     if(istat/=0) call memerr(two,'NBPOI','ibm_wallib',0_ip)

     allocate(pbobo(nboun_tmp+1),stat=istat)
     call memchk(zero,istat,mem_servi(1:2,servi),'PBOBO','ibm_wallib',pbobo)
     !
     ! Compute Hash table (initialize, reset, add and destroy)
     !
     mbobo = mbpoi * mnodb
     msize = mbobo * nboun_tmp
     allocate(lista(msize),stat=istat)
     call memchk(zero,istat,mem_servi(1:2,servi),'LISTA','ibm_wallib',lista)
     call htaini( ht, mbobo )
     pbobo(1) = 1
     do iboun = 1, nboun_tmp
        call htares( ht, lista(pbobo(iboun):) )
        pblty = abs(ltypb_tmp(iboun))
        do inodb= 1, nnode(pblty)
           jnodb = lnodb_tmp(inodb,iboun)
           do jboun = pbopo(jnodb), pbopo(jnodb+1)-1
              kboun = lbopo(jboun)
              if( kboun /= iboun ) call htaadd( ht, kboun )
           end do
        end do
        pbobo(iboun+1) = pbobo(iboun) + ht%nelem
     end do
     call htades( ht )
     nedgb = pbobo(nboun_tmp+1)-1
     !
     ! Allocate memory and compute list of adjacancies LBOBO
     !
     allocate(lbobo(nedgb),stat=istat)
     call memchk(zero,istat,mem_servi(1:2,servi),'LBOBO','ibm_wallib',lbobo)
     do iboun = 1, nedgb
        lbobo(iboun) = lista(iboun)
     end do

     call memchk(two,istat,memor_dom,'LISTA','ibm_wallib',lista)
     deallocate(lista,stat=istat)
     if(istat/=0) call memerr(two,'LISTA','ibm_wallib',0_ip)

     call memchk(two,istat,memor_dom,'PBOPO','ibm_wallib',pbopo)
     deallocate(pbopo,stat=istat)
     if(istat/=0) call memerr(two,'PBOPO','ibm_wallib',0_ip)

     call memchk(two,istat,memor_dom,'LBOPO','ibm_wallib',lbopo)
     deallocate(lbopo,stat=istat)
     if(istat/=0) call memerr(two,'LBOPO','ibm_wallib',0_ip)

     !----------------------------------------------------------------------
     !
     ! Identify convex sets
     !
     !----------------------------------------------------------------------
     !
     ! Allocate memory
     !
     if ( nboun_tmp > 0 ) then
        allocate(lconv(nboun_tmp),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,servi),'LCONV','ibm_wallib',lconv)
        
        do iboun = 1,nboun_tmp       
           lconv(iboun) = -1_ip  
        end do
        !
        ! Loop over boundaries, iboun 
        !
        lconv(1)  = 1_ip
        sboun_tmp = 0_ip
        sboun     = 1_ip
        nconv     = 1_ip
        
        do while ( sboun < nboun_tmp )
           !
           ! Loop while the current convex set find at least a new face
           !        
           do while ( sboun > sboun_tmp)
              sboun_tmp = sboun
              !
              ! Loop over all the walls
              !        
              do iboun = 1,nboun_tmp       
                 if (lconv(iboun) == nconv) then
                    call runend(' ibm:wallib: BUG a corregir!!!!!!!!!!!')        
                    mblty =  ltypb(iboun)
                    pnodb =  nnode(mblty)
                    inorm(3) = 0.0_rp
                    call extbou(2_ip,pnodb,lnodb(1,iboun),coord,inorm)
                    !
                    ! Loop over neighbor boundaries, jboun 
                    !
                    do izbou = pbobo(iboun),pbobo(iboun+1)-1
                       
                       jboun = lbobo(izbou)
                       if (lconv(jboun) == -1) then                          
                          !
                          ! Test all nodes in jboun to check if this face is convex to iboun
                          !
                          pblty = ltypb(jboun)
                          pnodb = nnode(pblty)
                          ipoib = lnodb(1,iboun)                                                 
                          conve = 0_ip
                          do jnodb = 1,pnodb
                             jpoib = lnodb(jnodb,jboun)  
                             iconv = 0_ip                                                    
                             do idime = 1,ndime
                                jnorm(idime) = coord(idime,jpoib) - coord(idime,ipoib)
                             end do
                             call vecuni(ndime,jnorm,xfact)
                             pladi = 0.0_rp
                             do idime = 1,ndime
                                pladi = pladi + inorm(idime) * jnorm(idime)
                             end do
                             !
                             ! This tolerancy of -1.0e-5_rp is required to minimize the number of convex sets
                             ! 
                             if (pladi > -1.0e-5_rp ) Then
                                conve = conve + 1_ip
                             end if
                          end do
                          !
                          ! If all nodes in jboun are convex to iboun, then iboun and jboun are convex
                          !                            
                          if (conve == pnodb) then
                             sboun = sboun + 1
                             lconv(jboun) = nconv                       
                          end if
                          !pblty = ltypb(jboun)
                          !pnodb = nnode(pblty)                          
                          !call extbou(2_ip,pnodb,lnodb(1,jboun),coord,jnorm)
                          !conve = 0_ip
                          !do idime = 1,ndime
                          !   if (jnorm(idime) > inorm(idime)-1.0e-2_rp .and. jnorm(idime) < inorm(idime)+1.0e-2_rp) then
                          !      conve = conve + 1_ip
                          !   end if
                          !end do
                          !if (conve == ndime) then
                          !   sboun = sboun + 1
                          !   lconv(jboun) = nconv                       
                          !end if
                       end  if
                    end do
                 end if
              end do
           end do
           !
           ! Find a new set of convex faces
           !                            
           iboun = 1
           do while ( iboun <  nboun_tmp )
              iboun = iboun + 1
              if (lconv(iboun) == -1) then
                 sboun = sboun + 1
                 nconv = nconv + 1
                 lconv(iboun) = nconv
                 iboun = nboun_tmp
              end if
           end do
        end do

        nwaib = nconv
        allocate( twall_ibm(nwaib), stat = istat )
        


        do iwaib = 1,nwaib

           sboun = 0_ip
           do iboun = 1,nboun_tmp
              if (lconv(iboun) == iwaib) then
                 sboun = sboun + 1_ip
              end if
           end do
           
           twall_ibm(iwaib) % nboun = sboun
           allocate( twall_ibm(iwaib) % ltypb(sboun) , stat = istat ) 
           call memchk(zero,istat,memor_dom,'LTYPB','ibm_wallib',twall_ibm(iwaib) % ltypb)
           allocate( twall_ibm(iwaib) % lnodb(mnodb,sboun) , stat = istat ) 
           call memchk(zero,istat,memor_dom,'LNODB','ibm_wallib',twall_ibm(iwaib) % lnodb)
           allocate( twall_ibm(iwaib) % bouno(ndime,sboun) , stat = istat ) 
           call memchk(zero,istat,memor_dom,'BOUNO','ibm_wallib',twall_ibm(iwaib) % bouno)
           sboun = 0_ip
           do iboun = 1,nboun_tmp
              if (lconv(iboun) == iwaib) then
                 sboun = sboun + 1_ip
                 pnodb = nnode(ltypb_tmp(iboun))
                 twall_ibm(iwaib) % ltypb(sboun) = ltypb_tmp(iboun)
                 do inodb = 1,pnodb
                    twall_ibm(iwaib) % lnodb(inodb,sboun) = lnodb_tmp(inodb,iboun)
                 end do
                 do idime = 1,ndime
                    twall_ibm(iwaib) % bouno(idime,sboun) = bouno_tmp(idime,iboun)
                 end do
              end if
           end do
        end do
        !----------------------------------------------------------------------
        !
        ! Coordinates
        !
        !----------------------------------------------------------------------        
        do iwaib = 1,nwaib
           call memgen(1_ip,npoin,0_ip)
           sboun = 0_ip           
           do iboun = 1,nboun_tmp
              if (lconv(iboun) == iwaib) then
                 sboun = sboun + 1_ip
                 pnodb = nnode(twall_ibm(iwaib) % ltypb(sboun))
                 do inodb = 1,pnodb
                    ipoin = twall_ibm(iwaib) % lnodb(inodb,sboun)
                    gisca(ipoin) = 1
                 end do
              end if
           end do
           
           kpoin = 0
           do ipoin = 1,npoin
              if( gisca(ipoin) /= 0 ) then
                 kpoin = kpoin + 1
                 gisca(ipoin) = kpoin
              end if
           end do
           twall_ibm(iwaib) % npoin = kpoin
           allocate( twall_ibm(iwaib) % coord(ndime,kpoin) , stat = istat ) 
           call memchk(zero,istat,memor_dom,'COORD','ibm_wallib',twall_ibm(iwaib) % coord)
           allocate( twall_ibm(iwaib) % lninv(kpoin) , stat = istat ) 
           call memchk(zero,istat,memor_dom,'LNINV','ibm_wallib',twall_ibm(iwaib) % lninv)
           do ipoin = 1,npoin
              if( gisca(ipoin) /= 0 ) then
                 kpoin = gisca(ipoin)
                 twall_ibm(iwaib) % lninv(kpoin) = ipoin
                 do idime = 1,ndime
                    twall_ibm(iwaib) % coord(idime,kpoin) = coord(idime,ipoin)
                 end do
              end if
           end do
           
           sboun = 0_ip
           do iboun = 1,nboun_tmp
              if (lconv(iboun) == iwaib) then
                 sboun = sboun + 1_ip
                 pnodb = nnode(twall_ibm(iwaib) % ltypb(sboun))
                 do inodb = 1,pnodb
                    ipoin = gisca( twall_ibm(iwaib) % lnodb(inodb,sboun) )
                    twall_ibm(iwaib) % lnodb(inodb,sboun) = ipoin
                 end do
              end if
           end do
           call memgen(3_ip,npoin,0_ip)
        end do
        !
        ! Deallocate memory
        !
        call memchk(two,istat,memor_dom,'LCONVE','ibm_wallib',lconv)
        deallocate(lbopo,stat=istat)
        
        !----------------------------------------------------------------------
        !
        ! Tree structure
        !
        !----------------------------------------------------------------------
        
        if( INOTMASTER ) then
           do iwaib = 1,nwaib
              call kdtree(&
                   1_ip,mnodb,twall_ibm(iwaib) % npoin,twall_ibm(iwaib) % nboun,&
                   twall_ibm(iwaib) % coord,twall_ibm(iwaib) % lnodb,twall_ibm(iwaib) % ltypb,&
                   twall_ibm(iwaib) % fabox,twall_ibm(iwaib) % bobox,twall_ibm(iwaib) % sabox,&
                   twall_ibm(iwaib) % blink,twall_ibm(iwaib) % stru2,twall_ibm(iwaib) % ldist,&
                   twall_ibm(iwaib) % lnele )
           end do
        end if
        
        !----------------------------------------------------------------------
        !
        ! Deallocate memory
        !
        !----------------------------------------------------------------------
        
        call memchk(two,istat,memor_dom,'PBOBO','ibm_wallib',pbobo)
        deallocate(pbobo,stat=istat)
        if(istat/=0) call memerr(two,'PBOBO','ibm_wallib',0_ip)
        
        call memchk(two,istat,memor_dom,'LBOBO','ibm_wallib',lbobo)
        deallocate(lbobo,stat=istat)
        if(istat/=0) call memerr(two,'LBOBO','ibm_wallib',0_ip)
        
        call memchk(two,istat,memor_dom,'LNODB_TMP','ibm_wallib',lnodb_tmp)
        deallocate(lnodb_tmp,stat=istat)
        if(istat/=0) call memerr(two,'LNODB_TMP','ibm_wallib',0_ip)
        
        call memchk(two,istat,memor_dom,'LTYPB_TMP','ibm_wallib',ltypb_tmp)
        deallocate(ltypb_tmp,stat=istat)
        if(istat/=0) call memerr(two,'LTYPB_TMP','ibm_wallib',0_ip)
        
        call memchk(two,istat,memor_dom,'LWAIB','ibm_wallib',lwaib)
        deallocate(lwaib,stat=istat)
        if(istat/=0) call memerr(two,'LWAIB','ibm_wallib',0_ip)
     else
        nwaib = 0_ip
     end if
  end if

  !----------------------------------------------------------------------
  !
  ! Unify wall convexization numbering
  !
  !----------------------------------------------------------------------

  !allocate( nwaib_tmp(npart+1) )
  !do ii = 1,npart+1
  !   nwaib_tmp(ii) = 0
  !end do
  !nwaib_tmp(kfl_paral+1) = nwaib
  !call parari('SUM',0_ip,npart+1,nwaib_tmp)
  !do ii = 1,npart
  !   nwaib_tmp(ii+1) = nwaib_tmp(ii+1) + nwaib_tmp(ii) 
  !end do
  !if( INOTMASTER ) then
  !   do iwaib = 1,nwaib
  !      twall_ibm(iwaib) % numbe = iwaib + nwaib_tmp(kfl_paral)
  !   end do

  !   do iboun = 1,nboun_2
  !      iwaib = lconv(iboun)
  !      lconv(iboun) = twall_ibm(iwaib) % numbe
  !   end do

     !do ineig = 1,nneig
        !nbous = frcom(ineig) % nbcos
        !do ibous = 1,nbous
        !   frcom(ineig) % lbcos
        !end do
     !end do

  !end if
  !call runend('POPO')

end subroutine ibm_wallib

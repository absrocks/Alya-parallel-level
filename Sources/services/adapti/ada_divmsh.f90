subroutine ada_divmshNEW(itask)
!-----------------------------------------------------------------------
!****f* adapti/ada_divmsh
! NAME 
!    ada_divmsh
! DESCRIPTION
!    This routine performs adaptivity mesh division 
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_elmtyp

  use      def_adapti
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: &
       ielem,jelem,kelem,inode,jnode,knode,pelty,pnode,ipoin,&
       iface,jface,kface,jjfac,ifacs,pnodb,&
       idime,ifoun,icoun,ncoun,ipone,jj,ii,imark,ifada
  integer(ip)             :: &
       lnomo(2*mnode),lnodi(mnode),lnodj(mnode),lnoei(mnode),lnoej(mnode),lnaux(mnode)
  real(rp)                :: &
       elcod(ndime,mnode),coaux(ndime)

 
  if (kfl_redgr_ada >= 1) then

     if (itask == ITASK_TURNON) then           
        !
        ! kfl_redgr_ada = 1 --> Non-defined error estimator: divide all of the elements
        !
        do ielem=1,nelem
           pelty=ltype(ielem)
           pnode=nnode(pelty)
           lnseg_ada(ielem)                   = 0      ! segmented element mark
           do inode= 1,pnode
              lnoed_ada(inode,inode  ,ielem)  = 0
              lnoed_ada(inode,inode  ,ielem)  = 0
              lnoed_ada(inode,mnode+1,ielem)  = 0
              lnoed_ada(inode,mnode+2,ielem)  = 0
              lnoed_ada(inode,mnode+3,ielem)  = 4 * (ndime - 1)
           end do
        end do
     end if

     if (itask == ITASK_ENDSTE) then
        
        call ada_livinf('Red-green adaptivity partition -> Proceed.')
        call ada_livinf('Red partition -> Proceed.')

        if (ndime == 2) then
           if (pelty /= TRI03) call runend('ADA_DIVMSH: REDGREEN ONLY FOR SIMPLEX ELEMENTS!')
        else if (ndime == 3) then
           if (pelty /= TET04) call runend('ADA_DIVMSH: REDGREEN ONLY FOR SIMPLEX ELEMENTS!')
        end if
        
        ipone = 0
        do ifacs= 1,nfacs
           pnodb= lfacs(mnodb+3,ifacs)

           ! identify red faces (i.e. faces to be segmented)
           imark = 1           
           call ada_markel('FACES',ifacs,imark)
           


        end do
        
        npnew_ada = ipone        
        
        call ada_livinf('Red partition -> Done.')

     end if

  end if


end subroutine ada_divmshNEW


!----------- vieja divmsh:

subroutine ada_divmsh(itask)
!-----------------------------------------------------------------------
!****f* adapti/ada_divmsh
! NAME 
!    ada_divmsh
! DESCRIPTION
!    This routine performs adaptivity mesh division 
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_elmtyp

  use      def_adapti
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: &
       ielem,jelem,kelem,inode,jnode,knode,pelty,pnode,ipoin,&
       iface,jface,kface,jjfac,&
       idime,ifoun,icoun,ncoun,ipone,jj,ii,imark,ifada
  integer(ip)             :: &
       lnomo(2*mnode),lnodi(mnode),lnodj(mnode),lnoei(mnode),lnoej(mnode),lnaux(mnode)
  real(rp)                :: &
       elcod(ndime,mnode),coaux(ndime)

 
  if (kfl_redgr_ada >= 1) then

     if (itask == ITASK_TURNON) then           
        !
        ! kfl_redgr_ada = 1 --> Non-defined error estimator: divide all of the elements
        !
        do ielem=1,nelem
           lnseg_ada(ielem)                   = 0      ! segmented element mark
           do inode= 1,pnode
              lnoed_ada(inode,inode  ,ielem)  = 0
              lnoed_ada(inode,inode  ,ielem)  = 0
              lnoed_ada(inode,mnode+1,ielem)  = 0
              lnoed_ada(inode,mnode+2,ielem)  = 0
              lnoed_ada(inode,mnode+3,ielem)  = 4 * (ndime - 1)
           end do
        end do
     end if

     if (itask == ITASK_ENDSTE) then
        
        call ada_livinf('Red-green adaptivity partition -> Proceed.')
        call ada_livinf('Red partition -> Proceed.')

        !
        ! Set the coordinates of the new nodes
        !

        ipone = 0

        do ielem=1,nelem
           pelty=ltype(ielem)
           pnode=nnode(pelty)
           if (ndime == 2) then
              if (pelty /= TRI03) call runend('ADA_DIVMSH: REDGREEN ONLY FOR SIMPLEX ELEMENTS!')
           else if (ndime == 3) then
              if (pelty /= TET04) call runend('ADA_DIVMSH: REDGREEN ONLY FOR SIMPLEX ELEMENTS!')
           end if
           
           ! identify red elements (i.e. elements to be segmented)
           imark = 1           
           call ada_markel('CELLS',ielem,imark)

           if (imark == 1) then

              lnseg_ada(ielem)= lnseg_ada(ielem) + 1              

              do inode=1,pnode
                 lnomo(inode      )= inode
                 lnomo(inode+pnode)= inode
                 ipoin= lnods(inode,ielem)
                 do idime=1,ndime
                    elcod(idime,inode)= coord(idime,ipoin)          
                 end do
                 lnodi(inode)= ipoin
              end do
              
              jface= 0
              kface= 0
              ii= pelel(ielem) - 1 
              
              iface= 1
              do while (kface == 0)
                 ii= ii + 1
                 if (ii == (pelel(ielem+1)-1)) kface = -1
                 jelem= lelel(ii)                 
                 
                 ifoun= 0              
                 ncoun= 0
                 if (lnoed_ada(iface,mnode+1,ielem) > 0) then
                    ncoun = -1      ! do not search, already found
                 end if
                 
                 if (ncoun == 0) then
                    do jnode= 1,pnode
                       lnodj(jnode)= lnods(jnode,jelem)
                       lnoei(jnode)= 0
                       lnoej(jnode)= 0
                    end do
                    
                    inode = 1
                    do while (ncoun < pnode-1)
                       ifoun= 0
                       jnode= 1
                       do while (ifoun==0)
                          if (lnodi(lnomo(inode)) == lnodj(lnomo(jnode))) then
                             ifoun= 1
                             ncoun= ncoun+1
                             lnoei(ncoun)= inode
                             lnoej(ncoun)= jnode
                             jnode= pnode
                          end if
                          jnode = jnode + 1
                          if (jnode > pnode) exit
                       end do
                       inode = inode + 1
                       if (inode > pnode) exit
                    end do
                    
                    jface= 1
                    do jj= pelel(jelem) , pelel(jelem+1)-1
                       if (lelel(jj) == ielem) exit
                       jface= jface+1
                    end do
                    
                    if ( lnoed_ada(iface,mnode+1,ielem) == 0 .and. &
                         lnoed_ada(jface,mnode+1,jelem) == 0 ) then
                       
                       ipone = ipone + 1
                       !
                       ! sorting lnoej and lnoei
                       !
                       if (ncoun == 2) then
                          if ((lnoej(1) - lnoej(2) == 1).or.(lnoej(2) - lnoej(1) == 2)) then 
                             jnode    = lnoej(1)
                             lnoej(1) = lnoej(2)
                             lnoej(2) = jnode
                          end if
                          if ((lnoei(1) - lnoei(2) == 1).or.(lnoei(2) - lnoei(1) == 2)) then 
                             jnode    = lnoei(1)
                             lnoei(1) = lnoei(2)
                             lnoei(2) = jnode
                          end if
                       else if (ncoun == 3) then
                          call runend('ADA_DIVMSH: 3D NOT YET PROGRAMMED!')
                       end if
                                              
                       coaux = 0.0_rp                                           
                       do icoun= 1,ncoun
                          lnoed_ada(iface,icoun,ielem) = lnodi(lnoei(icoun)) 
                          lnoed_ada(jface,icoun,jelem) = lnodj(lnoej(icoun)) 
                          do idime=1,ndime
                             coaux(idime) = coaux(idime) + 0.5_rp * elcod(idime,lnomo(lnoei(icoun)))
                          end do
                       end do
                       
                       lnoed_ada(iface,mnode+1,ielem) = jelem
                       lnoed_ada(iface,mnode+2,ielem) = ipone+npoin
                       
                       lnoed_ada(jface,mnode+1,jelem) = ielem
                       lnoed_ada(jface,mnode+2,jelem) = ipone+npoin
                       
                       conew_ada(1:ndime,ipone) =  coaux(1:ndime)
                       
                    end if
                 end if
                 
                 iface= iface+1
                 
              end do
              !
              ! Check if the element has boundary faces (no neighboring element) 
              !
              kface = iface - 1
              
              if (kface < pnode) then
                 
                 ncoun = pnode - 1
                 
                 ifoun= 0
                 iface= 0
                 lnaux= 0
                 do while (ifoun==0)
                    iface= iface+1
                    lnaux(iface) = lnoed_ada(iface,1,ielem)
                    if (lnoed_ada(iface,mnode+1,ielem)==0) then                       
                       
                       lnaux(iface) = 0
                       ipone = ipone + 1
                       coaux = 0.0_rp                    
                       
                       if (lnaux(iface) == 0) then
                          do inode= 1,pnode
                             ! ojo con el 3d!!!!!
                             knode = 0
                             do jnode= 1,pnode
                                if (lnodi(inode) == lnaux(jnode)) knode= jnode
                             end do
                             if (knode == 0) then
                                knode= inode 
                                exit
                             end if
                          end do
                       else
                          inode= 0
                          do jnode= 1,pnode
                             ! ojo con el 3d!!!!!
                             if (lnoed_ada(iface,1,ielem) == lnods(jnode,ielem)) inode= jnode
                          end do
                       end if
                       
                       do icoun= 1,ncoun
                          lnoei(icoun) = lnomo(inode+icoun-1)
                       end do
                       
                       do icoun= 1,ncoun
                          lnoed_ada(iface,icoun,ielem) = lnodi(lnoei(icoun)) 
                          do idime=1,ndime
                             coaux(idime) = coaux(idime) + 0.5_rp * elcod(idime,lnoei(icoun))
                          end do
                       end do
                       
                       lnoed_ada(iface,mnode+1,ielem) = -1
                       lnoed_ada(iface,mnode+2,ielem) = ipone+npoin
                       lnaux(iface) = lnoed_ada(iface,1,ielem)
                       
                       conew_ada(1:ndime,ipone) =  coaux(1:ndime)
                       
                    end if
                    if (iface==pnode) ifoun= -1
                 end do
                 
                 
              end if

           end if
           
        end do
        
        npnew_ada = ipone        
        
        call ada_livinf('Red partition -> Done.')

     end if

  end if


end subroutine ada_divmsh

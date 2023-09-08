subroutine ada_newcon(itask)
!-----------------------------------------------------------------------
!****f* adapti/ada_newcon
! NAME 
!    ada_newcon
! DESCRIPTION
!    This routine computes the new connectivity
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
       ielem,jelem,kelem,ielne,kelne,inode,jnode,knode,pelty,pnode,&
       ipoin,jpoin,kpoin,iface,jface,kface,jjfac,&
       idime,ifoun,icoun,ncoun,ipone,jj,ii,imark,iechi,kechi
  integer(ip)             :: &
       lnomo(2*mnode),lnodi(mnode),lnaux(mnode)
  real(rp)                :: &
       elcod(ndime,mnode),coaux(ndime),rpnod
  
  if (kfl_redgr_ada >= 1) then

     if (itask == ITASK_ENDSTE) then
        
        call ada_livinf('Green partition -> Proceed.')
        !
        ! Set the new connectivities
        !
        ipone = npnew_ada
        ielne = 0
        kelne = 0
        do ielem= 1,nelem
           pelty=ltype(ielem)
           pnode=nnode(pelty)           

           do inode=1,pnode
              lnomo(inode      )= inode
              lnomo(inode+pnode)= inode
              lnaux(inode      )= 0
           end do
           
           kechi = 0
           if (lnseg_ada(ielem) == 0) then
              !
              ! Non-segmented element
              !
              ncoun= 0
              if (kfl_redgr_ada >= 10) then
                 do iface= 1,pnode
                    if (lnoed_ada(iface,1,ielem) > 0 ) then
                       ncoun= ncoun+1
                       lnaux(ncoun)= iface
                    end if
                 end do
              end if

              if (ncoun==0) then
                 !
                 ! No segmented neighbours ->  keep the element as it is
                 !
                 ielne= ielne + 1
                 do inode= 1,pnode
                    lnnew_ada(inode,ielne) =  lnods(inode,ielem)
                 end do
                 !
                 ! Put to zero lfaso_ada
                 !
                 do iechi= 1,nechi_ada
                    lfaso_ada(iechi,ielem) = 0
                 end do
              else

                 !
                 ! Segmented neighbours -> conformize the element
                 !
                 ! 1. Add a central node                                  
                 ipone = ipone + 1
                 rpnod = 1.0_rp / real(pnode)
                 coaux = 0.0_rp                                                            
                 do inode= 1,pnode
                    ipoin= lnods(inode,ielem)
                    do idime=1,ndime
                       coaux(idime) = coaux(idime) + rpnod * coord(idime,ipoin)
                    end do
                 end do
                 conew_ada(1:ndime,ipone) =  coaux(1:ndime)

                 lnseg_ada(ielem) = - (ipone+npoin)     ! label for central-noded elements

                 if (ndime == 2) then
                    ! 2. Connect the segmented faces to the central node
                    do icoun= 1,ncoun
                       iface= lnaux(icoun)
                       ielne= ielne+1
                       kechi= kechi+1
                       lfaso_ada(kechi,ielem) = ielne
                       lsofa_ada(ielne      ) = ielem
                       lnnew_ada(1,ielne) = lnoed_ada(iface,      1,ielem)
                       lnnew_ada(2,ielne) = lnoed_ada(iface,mnode+2,ielem)
                       lnnew_ada(3,ielne) = ipone + npoin 
                       ielne= ielne+1
                       kechi= kechi+1
                       lfaso_ada(kechi,ielem) = ielne
                       lsofa_ada(ielne      ) = ielem
                       lnnew_ada(1,ielne) = lnoed_ada(iface,      2,ielem)
                       lnnew_ada(2,ielne) = ipone + npoin                               
                       lnnew_ada(3,ielne) = lnoed_ada(iface,mnode+2,ielem)
                       lnaux(icoun)       = lnoed_ada(iface,      1,ielem)
                    end do
                    
                    ! 3. Connect the non-segmented faces to the central node
                    
                    do inode= 1,pnode                       
                       ifoun= 0
                       do icoun= 1,ncoun
                          if (lnaux(icoun)==lnods(inode,ielem)) then
                             ifoun= 1
                          end if
                       end do
                       if (ifoun == 0) then
                          ielne= ielne+1
                          kechi= kechi+1
                          lfaso_ada(kechi,ielem) = ielne
                          lsofa_ada(ielne      ) = ielem
                          jnode= lnomo(inode+1)
                          lnnew_ada(1,ielne) = lnods(inode,ielem)
                          lnnew_ada(2,ielne) = lnods(jnode,ielem)
                          lnnew_ada(3,ielne) = ipone + npoin 
                       end if
                    end do

                 else if (ndime==3) then
                    

                 end if

              end if

           else if (lnseg_ada(ielem) == 1) then
              !
              ! Segmented element
              !
              do inode= 1,pnode
                 ipoin = lnods(inode,ielem)
                 do iface= 1,pnode
                    if (lnoed_ada(iface,1,ielem) == ipoin) then                    
                       lnaux(inode) = iface
                    end if
                 end do
              end do
              
              if (ndime == 2) then
                 kelne = ielne + 4
                 do inode= 1,pnode
                    jnode= lnomo(inode+pnode-1)
                    ielne= ielne + 1
                    lnnew_ada(    1,ielne) =  lnods(inode,ielem)
                    lnnew_ada(    2,ielne) =  lnoed_ada(lnaux(inode),mnode+2,ielem)
                    lnnew_ada(    3,ielne) =  lnoed_ada(lnaux(jnode),mnode+2,ielem)              
                    lnnew_ada(inode,kelne) =  lnoed_ada(lnaux(inode),mnode+2,ielem)
                 end do
                 ielne = ielne + 1
                 kechi= kechi+1
                 lfaso_ada(kechi,ielem) = ielne
                 lsofa_ada(ielne      ) = ielem
              else 
                 call runend('ADA_NEWCON: 3D NOT YET PROGRAMMED!')
              end if
              
           end if

        end do
        
        ! total value of elements and points in the adapted mesh
        nenew_ada = ielne
        npnew_ada = ipone + npoin

        ! store current lnseg to lsori
        do ielem=1,nelem
           lsori_ada(ielem)= lnseg_ada(ielem)
        end do
        !
        ! Make the announcement of success
        !
        call ada_livinf('Red-green adaptivity partition done.')
        call ada_livinf('Mesh adapted... leaving.')
        
     end if

  end if


end subroutine ada_newcon

subroutine ibm_fielib(ipoin,propo,lnode,shapl,limit,dista)   
  !--------------------------------------------------------------------------
  !****f* ibm_fielib
  ! NAME
  !    ibm_fielib
  ! DESCRIPTION
  !    Find the best element formed by a given node and their neiborgs that 
  !    contains the point of intersection with the particle surface    
  !
  !    Returns the point of prejection, a list of elemental nodes 
  !    and the value of the shape functions
  ! USED BY
  !    nastin/nsi_embedd
  !--------------------------------------------------------------------------
  use def_kintyp,  only    :  ip,rp
  use def_immbou
  use def_master
  use def_kermod
  use def_domain
  use mod_kdtree
  use def_elsest, only     :  lmini,lmaxi
  use mod_elmgeo, only  :  elmgeo_natural_coordinates 
  implicit none


  integer(ip), intent(in)     :: ipoin
  integer(ip), intent(out)    :: lnode(ndime+1),limit
  integer(ip)                 :: jpoin,iimbo,idime,izdom,ilist,nlist
  integer(ip)                 :: ielem,jelem,inode
  integer(ip)                 :: nlele,iter,count,ifoun,chang,tempo,pnode
  integer(ip)                 :: kfl_lelem,updat,kface
  integer(ip)                 :: lnod1(ndime),lnod2(ndime)
  integer(ip), pointer        :: list(:),lelem(:,:)
  integer(4)                  :: istat
  real(rp),    intent(in)     :: propo(ndime),dista
  real(rp),    intent(out)    :: shapl(ndime+1)
  real(rp)                    :: coloc(3),Q,bestQ
  real(rp)                    :: Qmax,total,dist2,facto
  real(rp)                    :: deriv(ndime,mnode),shapf(mnode),proje(ndime)
  real(rp)                    :: elcod(ndime,mnode)


  !
  ! Initialize
  !
  Qmax = 1.0e9_rp
  do idime = 1,ndime+1
     lnode(idime) = 0_ip
     shapl(idime) = 0.0_rp
  end do

  pnode     = ndime+1
  kfl_lelem = 0
  limit     = ndime+1

  !
  ! Neighbor nodes (only one level)
  !     
  ilist = 0
  nlist = r_dom(ipoin+1) - r_dom(ipoin) 
  if (npoin_2 > npoin .and. ipoin >= npoi1+1 .and. ipoin <= npoin) then
     nlist = nlist + r_dom_2(ipoin-npoi1+1) - r_dom_2(ipoin-npoi1) 
  end if
  allocate( list(nlist), stat = istat )

  do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
     jpoin = c_dom(izdom)
     if ( lntib(jpoin) == 0 ) then        
        ilist       = ilist + 1
        list(ilist) = jpoin
     end if
  end do
  if (npoin_2 > npoin .and. ipoin >= npoi1+1 .and. ipoin <= npoin) then
     do izdom = r_dom_2(ipoin-npoi1),r_dom_2(ipoin-npoi1+1)-1
        jpoin = c_dom_2(izdom)           
        if ( lntib(jpoin) == 0 ) then
           ilist       = ilist + 1
           list(ilist) = jpoin
        end if
     end do
  end if
  nlist = ilist
  !
  ! Possible elements
  !
  nlele = 0_ip
  if( ndime == 2 ) then
     do iter = 1,nlist-1
        nlele = nlele + iter
     end do
     if( nlele > 0 ) then
        allocate(lelem(3,nlele))
        kfl_lelem = 1
        count = 0
        do ielem = 1,nlist-1
           do jelem = ielem+1,nlist
              count          = count + 1
              lelem(1,count) = ipoin
              lelem(2,count) = list(ielem)
              lelem(3,count) = list(jelem)
           end do
        end do
     else
        nlele = 0
     end if
  else if( ndime == 3 ) then
     do iter = 1,nlist-2
        nlele = nlele + iter
     end do
     if (nlele > 0) then
        allocate(lelem(4,nlele))
        kfl_lelem = 1
        count = 0
        do ielem = 1,nlist-2
           do jelem = ielem+2,nlist
              count          = count + 1
              lelem(1,count) = ipoin
              lelem(2,count) = list(ielem)
              lelem(3,count) = list(ielem+1)
              lelem(4,count) = list(jelem)
           end do
        end do
     else
        nlele = 0
     end if
  end if

  bestQ = -1.0_rp
  do ielem = 1,nlele

     if (dista > 0.0_rp ) then     
        do idime = 1,ndime
           elcod(idime,1) = propo(idime)
        end do
        do inode = 2,pnode
           do idime = 1,ndime
              elcod(idime,inode) = coord(idime,lelem(inode,ielem))
           end do
        end do
        call ibm_qualib(elcod,chang,Q,Qmax)  
        !
        ! If the element is well formed
        ! 
        if( Q < Qmax ) then
           if( chang == 1 ) then
              tempo          = lelem(2,ielem)
              lelem(2,ielem) = lelem(3,ielem)
              lelem(3,ielem) = tempo
           end if
           !call elmgeo_natural_coordinates(      &
           !     ndime,pelty,pnode,elcod,shapf,   &
           !     deriv,xforc_material_nsi(12:,pmate),coloc,ifoun)
           call runend('IBM NOT CODED')
           !call elsest_chkelm(&
           !     ndime,1_ip,pnode,elcod,shapf,deriv, &
           !     coord(1:ndime,ipoin),coloc,ifoun,lmini,lmaxi)
        else
           ifoun = 0_ip
        end if
     else
        do inode = 1,pnode
           do idime = 1,ndime
              elcod(idime,inode) = coord(idime,lelem(inode,ielem))
           end do
        end do
        call ibm_qualib(elcod,chang,Q,Qmax)
        !
        ! If the element is well formed
        ! 
        if( Q < Qmax ) then
           if (chang == 1) then
              tempo          = lelem(2,ielem)
              lelem(2,ielem) = lelem(3,ielem)
              lelem(3,ielem) = tempo
           end if
           call runend('IBM NOT CODED')
           !call elsest_chkelm(&
           !     ndime,1_ip,pnode,elcod,shapf,deriv, &
           !     propo,coloc,ifoun,lmini,lmaxi)
        else
           ifoun = 0_ip
        end if

     end if

     if ( shapf(1) >= bestQ-1.0e-3_rp .and. ifoun == 1 ) then             
        updat = 1_ip
        if ( shapf(1) <= bestQ + 1.0e-3_rp) then
           updat = 0_ip
           do inode = 2,pnode
              lnod1(inode-1) = lninv_loc( lnode(inode) )
              lnod2(inode-1) = lninv_loc( lelem(inode,ielem) )
           end do
           call sortin(pnode-1,lnod1)
           call sortin(pnode-1,lnod2)                 

           inode = 0_ip
           do while ( inode < pnode-1_ip)
              inode = inode + 1_ip
              if (lnod2(inode) < lnod1(inode)) then
                 updat = 1_ip
                 inode = pnode
              elseif (lnod2(inode) > lnod1(inode)) then
                 inode = pnode
              end if
           end do
        end if
        if (updat == 1) then
           do inode = 1,pnode
              lnode(inode) = lelem(inode,ielem)
              shapl(inode) = shapf(inode)
           end do
           bestQ = shapl(1)
        end if
     end if
  end do



  if ( nlele == 0 .or. bestQ <= 0.0_rp ) then    
     !
     ! If there is not any element that contains the point of intersection
     !         
     bestQ    =  -1.0_rp
     iimbo = abs(lntib(ipoin))


     !
     ! Find a line with the shortest distance from its own surface point of intersection
     ! to the original one.
     !                 
     do ilist = 1,nlist
        facto = 0.0_rp
        if (dista > 0.0) then
           facto = -3.0_rp ! Use to find the interseccion when the segment doesn't intersect the particle (i.e.: Lohner method)
        end if
        call faceli(&
             imbou(iimbo) % sabox,imbou(iimbo) % blink,imbou(iimbo) %ltyib,imbou(iimbo) %lnoib,imbou(iimbo) %cooib, & 
             coord(1,ipoin),coord(1,list(ilist)),mnoib,imbou(iimbo) % nboib,facto,proje,kface,dist2)

        if (dista > 0.0) then
           total = 0.0_rp
           dist2 = 0.0_rp
           do idime = 1,ndime
              total = total + (coord(idime,list(ilist)) - proje(idime))** 2.0_rp
              dist2 = dist2 + (coord(idime,list(ilist)) - coord(idime,ipoin))** 2.0_rp
           end do
           total = sqrt(total)
           dist2 = sqrt(dist2)
           dist2 = dist2/total
        end if
        dist2 = 1.0_rp-dist2
   
        if ( dist2 > (bestQ-1.0e-3_rp) .and. kface /= -100 ) then
           updat = 1_ip
           if ( dist2 <= (bestQ+1.0e-3_rp) ) then
              updat = 0_ip
              lnod1(1) = lninv_loc( lnode(2)    )
              lnod2(1) = lninv_loc( list(ilist) )                 
              if (lnod2(1) < lnod1(1)) then
                 updat = 1_ip
              end if
           end if
           if (updat == 1) then
              lnode(1)       = ipoin
              lnode(2)       = list(ilist)
              shapl(1)       = dist2
              bestQ          = dist2
              shapl(2)       = 1.0_rp-dist2
              ! The intersection point in the line is the new point of intersection
              !do idime = 1,ndime
              !   propo(idime) = proje(idime)
              !end do
           end if
        end if
        limit = 2
     end do
  end if

  deallocate(list,stat=istat)
  if( kfl_lelem == 1 ) deallocate(lelem,stat=istat)


end subroutine ibm_fielib


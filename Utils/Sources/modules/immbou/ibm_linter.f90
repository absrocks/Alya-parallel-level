subroutine ibm_linter()
  !-----------------------------------------------------------------------
  !****f* ibm_linter/ibm_linter
  ! NAME
  !    ibm_linter
  ! DESCRIPTION
  !    This routines find the neighbor nodes for all the fringe nodes
  !    and also determines the kriging coefficients to interpolate  
  !    the value of the particle velocity using all these nodes.
  ! USED BY
  !    ibm_solite
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_immbou
  use mod_kdtree
  use mod_intpol
  use mod_messages, only : livinf

  implicit none
  integer(ip)         :: idime,ipoin,jpoin
  integer(ip)         :: izdom
  integer(ip)         :: iinte,ninte
  integer(ip)         :: iimbo,iboib,dummi
  integer(ip)         :: inode,ihole,isout,pnode
  integer(ip)         :: ielem,ielpo,minim,mini2
  integer(ip)         :: nele1,nele2
  real(rp)            :: dista,norfa(ndime),propo(ndime),total
  real(rp)            :: coor1(ndime),dumma(ndime),dummr

  if( kfl_diric_ibm /= 0 ) then

     if( ittim /= 0 ) call livinf(165_ip,'INT,',0_ip) 

     if( INOTMASTER ) then
        !
        ! Deallocate the interpolation arrays from the previous step
        !
        do ipoin = 1,npoin
           if (lnint(ipoin) % limit > 0_ip) then
              deallocate(lnint(ipoin) % lnode)
              deallocate(lnint(ipoin) % shapl)
           end if
        end do

        do ipoin = 1,npoin
           if (lnin2(ipoin) % limit > 0_ip) then
              deallocate(lnin2(ipoin) % lnode)
              deallocate(lnin2(ipoin) % shapl)
           end if
        end do


        do ipoin = 1,npoin
           lnint(ipoin) % limit = 0_ip
           lnin2(ipoin) % limit = 0_ip           
           if( lntib(ipoin) < 0 ) then            
              !
              ! Projection point on the surface mesh
              !     
              iimbo     = abs(lntib(ipoin))
              call dpopar(1_ip,coord(1:ndime,ipoin),&
                   imbou(iimbo) % npoib,mnoib,imbou(iimbo) % nboib,1.0e10_rp,&
                   imbou(iimbo) % ltyib,imbou(iimbo) % lnoib,imbou(iimbo) % cooib,&
                   dista,norfa,propo,iboib,imbou(iimbo) % fabox,imbou(iimbo) % sabox,&
                   imbou(iimbo) % blink,imbou(iimbo) % struc,imbou(iimbo) % ldist,&
                   imbou(iimbo) % lnele)

              !----------------------------------------------------------------------
              !
              ! Interpolation with elements
              !
              !---------------------------------------------------------------------- 
 
              if (kfl_inter_ibm == 1) then
                 allocate( lnint(ipoin) % lnode(ndime+1) )
                 allocate( lnint(ipoin) % shapl(ndime+1) )                 
                 call ibm_fielib(ipoin,propo,lnint(ipoin) % lnode,lnint(ipoin) % shapl,lnint(ipoin) % limit,dista)

                 if (dista < 0.0_rp) then
                    if (lnint(ipoin) % shapl(1) < 0.2_rp) then
                       total = 1.0_rp - (0.2_rp - lnint(ipoin) % shapl(1))
                       lnint(ipoin) % shapl(1) = 0.2_rp
                       do iinte = 2,lnint(ipoin) % limit
                          lnint(ipoin) % shapl(iinte) = lnint(ipoin) % shapl(iinte)*total
                       end do
                    end if

                    do iinte = 2,lnint(ipoin) % limit
                       lnint(ipoin) % shapl(iinte) = -lnint(ipoin) % shapl(iinte)/lnint(ipoin) % shapl(1)
                    end do
                    lnint(ipoin) % shapl(1) = 1.0_rp/lnint(ipoin) % shapl(1)
                 end if
              !----------------------------------------------------------------------
              !
              ! Kriging interpolation
              !
              !----------------------------------------------------------------------  
              elseif (kfl_inter_ibm == 2) then   
                 !
                 ! Find the number of free neighbor nodes
                 !
                 ninte = 0_ip
                 do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                    jpoin = c_dom(izdom)
                    if ( lntib(jpoin) == 0 ) then
                       ninte = ninte + 1_ip
                    end if
                 end do
                 !
                 ! Find the number of free neighbor nodes in other subdomains
                 !
                 if (npoin_2 > npoin .and. ipoin >= npoi1+1 .and. ipoin <= npoin) then
                    do izdom = r_dom_2(ipoin-npoi1),r_dom_2(ipoin-npoi1+1)-1
                       jpoin = c_dom_2(izdom)           
                       if ( lntib(jpoin) == 0 ) then
                          ninte = ninte + 1
                       end if
                    end do
                 end if
                 !
                 ! Allocate the variables
                 !
                 allocate( lnint(ipoin) % lnode(ninte) )
                 allocate( lnint(ipoin) % shapl(ninte+1) )
                 lnint(ipoin) % limit = ninte
                 !
                 ! Find the free neighbor nodes
                 !
                 iinte = 0
                 do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                    jpoin = c_dom(izdom)                    
                    if ( lntib(jpoin) == 0 ) then
                       iinte = iinte + 1                          
                       lnint(ipoin) % lnode(iinte) = jpoin                       
                    end if
                 end do
                 !
                 ! Find the free neighbor nodes in other subdomains
                 !
                 if (npoin_2 > npoin .and. ipoin >= npoi1+1 .and. ipoin <= npoin) then
                    do izdom = r_dom_2(ipoin-npoi1),r_dom_2(ipoin-npoi1+1)-1
                       jpoin = c_dom_2(izdom)           
                       if ( lntib(jpoin) == 0 ) then
                          iinte = iinte + 1                          
                          lnint(ipoin) % lnode(iinte) = jpoin
                       end if
                    end do
                 end if
                 !
                 ! Find the kriging interpolation coefficients
                 !
                 if (lnint(ipoin) % limit > 1) then
                    iimbo    = abs(lntib(ipoin))
                    call dpopar(1_ip,coord(1:ndime,ipoin),&
                         imbou(iimbo) % npoib,mnoib,imbou(iimbo) % nboib,1.0e10_rp,&
                         imbou(iimbo) % ltyib,imbou(iimbo) % lnoib,imbou(iimbo) % cooib,&
                         dummr,dumma,propo,dummi,imbou(iimbo) % fabox,imbou(iimbo) % sabox,&
                         imbou(iimbo) % blink,imbou(iimbo) % struc,imbou(iimbo) % ldist,&
                         imbou(iimbo) % lnele)

                    call krigin(coord,coord(:,ipoin),lnint(ipoin) % limit,lnint(ipoin) % lnode,lnint(ipoin) % shapl,propo)
                 else
                    !
                    ! In this situation, we only have one free neighbor node
                    !
                    do idime = 1,ndime
                       coor1(idime) =  coord(idime,ipoin) - 2.0_rp*(coord(idime,lnint(ipoin)%lnode(1)) - coord(idime,ipoin))
                    end do

                    call faceli(&
                         imbou(iimbo) % sabox,imbou(iimbo) % blink,imbou(iimbo) %ltyib,imbou(iimbo) %lnoib,imbou(iimbo) %cooib, & 
                         coor1,coord(1,lnint(ipoin)%lnode(1)),coord(1,ipoin),mnoib,imbou(iimbo) %nboib,propo,dummi,dista)

                    lnint(ipoin) % shapl(1) = (coord(1,ipoin) -  propo(1)) / ((coord(1,lnint(ipoin)%lnode(1)) - propo(1)))
                    lnint(ipoin) % shapl(2) = 1.0_rp - lnint(ipoin) % shapl(1)
                 end if
                 !
                 ! For each fringe node, we find all the elements store in different subdomains 
                 ! that use each of their free neighbor nodes.                 
                 !
                 ! Then, we multiply the kriging coefficients of the free neighbor nodes
                 ! by the number of found elements store in my subdomain and divide 
                 ! by the sum of all found elements. 
                 !
                 ! This allow us to do a parallel reduction of each free neigbor node correctly.
                 !
                 if (npoin_2 > npoin .and. ipoin > npoi1) then 
                    pnode = nnode(ltype(1)) 
                    do iinte = 1,lnint(ipoin) % limit
                       nele1 = 0_ip; nele2 = 0_ip
                       jpoin = lnint(ipoin) % lnode(iinte)                       
                       if (jpoin > npoi1 .and. jpoin <= npoin) then
                          do ielpo = pelpo(ipoin),pelpo(ipoin+1)-1
                             ielem = lelpo(ielpo)
                             do inode = 1,pnode
                                if ( lnods(inode,ielem) == jpoin ) nele1 = nele1 + 1
                             end do
                          end do
                          do ielpo = pelpo_2(ipoin-npoi1),pelpo_2(ipoin-npoi1+1)-1 
                             ielem = lelpo_2(ielpo)
                             do inode = 1,pnode
                                if ( lnods(inode,ielem) == jpoin ) nele2 = nele2 + 1
                             end do
                          end do
                          lnint(ipoin) % shapl(iinte) = lnint(ipoin) % shapl(iinte) &
                               * real(nele1,8) * (1.0_rp/real(nele1+nele2,8))
                       end if
                    end do
                 end if
              end if
           end if
           !----------------------------------------------------------------------
           !
           ! Interpolation for ALE. Free and fringe nodes with mesh movement.
           !
           !----------------------------------------------------------------------  
           if (lntib(ipoin) <= 0 .and. lpoty(ipoin) == 0) then
              !
              ! Find the number of free and fringe neighbor nodes
              !
              ninte = 0_ip
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)
                 if ( lntib(jpoin) == 0 .or. (lntib(jpoin) < 0 .and. lnti2(jpoin) <= 0) ) then
                    ninte = ninte + 1_ip
                 end if
              end do
              !
              ! Find the number of free and fringe neighbor nodes
              !
              if (npoin_2 > npoin .and. ipoin >= npoi1+1 .and. ipoin <= npoin) then
                 do izdom = r_dom_2(ipoin-npoi1),r_dom_2(ipoin-npoi1+1)-1
                    jpoin = c_dom_2(izdom)           
                    if ( lntib(jpoin) == 0 .or. (lntib(jpoin) < 0 .and. lnti2(jpoin) <= 0) ) then
                       ninte = ninte + 1
                    end if
                 end do
              end if

              allocate( lnin2(ipoin) % lnode(ninte) )
              allocate( lnin2(ipoin) % shapl(ninte+1) )
              lnin2(ipoin) % limit = ninte
              !
              ! Find the free and fringe neighbor nodes
              !
              iinte = 0
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)
                 if ( lntib(jpoin) == 0 .or. (lntib(jpoin) < 0 .and. lnti2(jpoin) <= 0) ) then
                    iinte = iinte + 1                          
                    lnin2(ipoin) % lnode(iinte) = jpoin
                 end if
              end do
              !
              ! Find the free and fringe neighbor nodes in other subdomains
              !
              if (npoin_2 > npoin .and. ipoin >= npoi1+1 .and. ipoin <= npoin) then
                 do izdom = r_dom_2(ipoin-npoi1),r_dom_2(ipoin-npoi1+1)-1
                    jpoin = c_dom_2(izdom)           
                    if ( lntib(jpoin) == 0 .or. (lntib(jpoin) < 0 .and. lnti2(jpoin) <= 0) ) then
                       iinte = iinte + 1                          
                       lnin2(ipoin) % lnode(iinte) = jpoin
                    end if
                 end do
              end if
              !
              ! For each free node, we find all the elements store in different subdomains 
              ! that use each of their free and fringe neighbor nodes.                 
              !
              ! Then, we multiply the kriging coefficients of the free and fringe neighbor nodes
              ! by the number of found elements store in my subdomain and divide 
              ! by the sum of all the found elements. 
              !
              ! This allow us to do a parallel reduction of each free and fringe neigbor node correctly.
              !
              do iinte = 1,lnin2(ipoin) % limit
                 lnin2(ipoin) % shapl(iinte) = 1.0_rp
              end do
              if (npoin_2 > npoin .and. ipoin > npoi1) then
                 pnode = nnode(ltype(1))
                 do iinte = 1,lnin2(ipoin) % limit
                    nele1 = 0_ip; nele2 = 0_ip
                    jpoin = lnin2(ipoin) % lnode(iinte)
                    if (jpoin > npoi1 .and. jpoin <= npoin) then
                       do ielpo = pelpo(ipoin),pelpo(ipoin+1)-1
                          ielem = lelpo(ielpo)
                          do inode = 1,pnode
                             if ( lnods(inode,ielem) == jpoin ) nele1 = nele1 + 1
                          end do
                       end do
                       do ielpo = pelpo_2(ipoin-npoi1),pelpo_2(ipoin-npoi1+1)-1
                          ielem = lelpo_2(ielpo)
                          do inode = 1,pnode
                             if ( lnods(inode,ielem) == jpoin ) nele2 = nele2 + 1
                          end do
                       end do
                       lnin2(ipoin) % shapl(iinte) = real(nele1,8) * (1.0_rp/real(nele1+nele2,8))
                    end if
                 end do
              end if

           end if
        end do
     end if
  end if
end subroutine ibm_linter

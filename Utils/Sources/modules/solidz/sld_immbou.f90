!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_immbou.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   Coupling with Immbou: impose force FORCF: interpolate b.c.
!> @details Coupling with Immbou: impose force FORCF: interpolate b.c.
!> @} 
!-----------------------------------------------------------------------
subroutine sld_immbou(aa,bb,uu)
  use def_kintyp
  use def_domain
  use def_master
  use def_solidz

  use mod_gradie
  use mod_kdtree

  implicit none
  integer(ip)                :: inode,ipoin
  integer(ip)                :: idime,jdime,izdom,jpoin,idofn,iimbo
  integer(ip)                :: izdod

  real(rp),    intent(out)   :: aa(ndime,ndime,nzdom) 
  real(rp),    intent(out)   :: bb(ndime,npoin)
  real(rp)                   :: uu(ndime,npoin)
  real(rp)                   :: x(3),v(3),aad(3)

  integer(ip)                :: updat,index,kpoin,initi,limit
  real(rp)                   :: propo(ndime)
  integer(ip), pointer       :: lnode(:)
  real(rp),    pointer       :: shapl(:)

  integer(ip)                :: kdime,kboun
  real(rp)                   :: deter,temp1,temp2,temp3,invma(ndime,ndime)
  real(rp)                   :: tempo,dista,dumma(ndime),proje(ndime)

  !----------------------------------------------------------------------
  !
  ! Set pressure and velocity of hole nodes to zero
  !
  !----------------------------------------------------------------------  

  do ipoin = 1,npoin   

     if( lntib(ipoin) > 0 ) then
        !
        ! IZDOD: Diagonal
        !
        izdod = r_dom(ipoin) - 1
        jpoin = 0
        do while( jpoin /= ipoin )
           izdod = izdod + 1
           jpoin = c_dom(izdod)
        end do

        do idime = 1,ndime
           aad(idime) = aa(idime,idime,izdod)
           if( abs(aad(idime)) < zeror ) aad(idime) = 1.0_rp
        end do
        !
        ! Set line to zero
        !
        do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
           do idime = 1,ndime
              do jdime = 1,ndime
                 aa(jdime,idime,izdom) = 0.0_rp
              end do
           end do
        end do
        !
        ! Presrcibe value to zero
        !
        do idime = 1,ndime
           aa(idime,idime,izdod) = aad(idime)
        end do

        idofn = (ipoin-1)*ndime
        do idime = 1,ndime
           idofn           = idofn + 1
           bb(idime,ipoin) = aad(idime)*0.0_rp
        end do
     end if
  end do

  do iimbo = 1,nimbo    

     do ipoin = 1,npoin
        if( lntib(ipoin) == -iimbo ) then
           !
           ! Element use to interpolate the velocity
           !
           lnode => lnint(ipoin) % lnode
           shapl => lnint(ipoin) % shapl
           limit =  lnint(ipoin) % limit
           !
           ! IZDOD: Diagonal
           !
           izdod = r_dom(ipoin) - 1
           jpoin = 0
           do while( jpoin /= ipoin )
              izdod = izdod + 1
              jpoin = c_dom(izdod)
           end do
           do idime = 1,ndime
              Aad(idime) = Aa(idime,idime,izdod)
              if( abs(Aad(idime)) < zeror ) Aad(idime) = 1.0_rp
           end do
           !
           ! Set line to zero
           !
           do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
              jpoin = c_dom(izdom)
              if (jpoin <= npoi1 .or. (jpoin>= npoi2 .and. jpoin<= npoi3)) then           
                 do idime = 1,ndime
                    do jdime = 1,ndime
                       Aa(jdime,idime,izdom) = 0.0_rp
                    end do
                 end do
              end if
           end do

           !----------------------------------------------------------------------
           !
           ! Normal interpolation
           ! Set velocity line equal to shape functions of the new element
           !
           !---------------------------------------------------------------------- 

           do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
              jpoin = c_dom(izdom)
              if (jpoin <= npoi1 .or. (jpoin >= npoi2 .and. jpoin <= npoi3)) then
                 if (ipoin /= jpoin) then
                    kpoin = 0
                    index = 1
                    do while ( kpoin /= jpoin .and. index < limit )
                       index = index+1_ip
                       kpoin = lnode(index)                       
                    end do
                    if ( kpoin == jpoin ) then  
                       do idime = 1,ndime        
                          Aa(idime,idime,izdom) = -shapl(index) * Aad(idime)
                       end do
                    end if
                 else
                    do idime = 1,ndime                                  
                       Aa(idime,idime,izdom) = Aad(idime)
                    end do
                 end if
              end if
           end do

           if (ipoin <= npoi1 .or. (ipoin>= npoi2 .and. ipoin<= npoi3)) then              
              !
              ! Presrcibe velocity in the projection point: propo
              !
              v(1)     = 0.0_rp
              v(2)     = 0.0_rp
              v(3)     = 0.0_rp
              do idime = 1,ndime
                 bb(idime,ipoin) = v(idime) * shapl(1) * Aad(idime) 
              end do

              !----------------------------------------------------------------------
              !
              ! Normal interpolation using two levels of neighbors
              !
              !----------------------------------------------------------------------  

              do inode = 2,limit
                 updat = 1_ip        
                 do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                    jpoin = c_dom(izdom)                    
                    if ( lnode(inode) == jpoin ) then
                       updat = 0_ip
                    end if
                 end do
                 if (updat == 1) then
                    do idime = 1,ndime
                       bb(idime,ipoin) = bb(idime,ipoin) + du_sld(idime,lnode(inode)) * shapl(inode) * Aad(idime)
                    end do
                 end if
              end do

           end if
        end if
     end do
  end do

  !----------------------------------------------------------------------
  !
  ! Disconnect free and fringe nodes (IPOIN) from hole nodes (JPOIN)
  !
  !----------------------------------------------------------------------     

  do ipoin = 1,npoin     
     if( lntib(ipoin) <= 0 ) then  
        do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
           jpoin = c_dom(izdom)
           if( lntib(jpoin) > 0 ) then
              do idime = 1,ndime
                 do jdime = 1,ndime
                    Aa(jdime,idime,izdom) = 0.0_rp
                 end do
              end do
           end if
        end do
     end if
  end do

end subroutine sld_immbou

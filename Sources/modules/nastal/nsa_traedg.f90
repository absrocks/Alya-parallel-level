!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_traedg.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Free nodes in trailing edges (only 3D)
!> @details Free nodes in trailing edges (only 3D)
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_traedg
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use mod_memchk

  use def_nastal
  implicit none

  integer(ip)  :: &
       iboun,kboun,pblty,pnodb,ielem,nifmark,&
       ibopo,idime,ipoin,jpoin,kpoin,inodb
  real(rp)     :: venor(ndime),veta1(ndime),veta2(ndime),venoe(ndime),rauxi

  if( INOTMASTER ) then
     !
     ! Loop over boundaries to detect trailing edges comparing exnor and face normal 
     !
     boundaries_detect: do iboun = 1,nboun
        pblty = ltypb(iboun) 
        pnodb = nnode(pblty)
        ielem = lelbo(iboun)
        !
        ! Mark boundary faces with all nodes 20000
        !
        nifmark= 0
        do inodb=1,pnodb
           ipoin= lnodb(inodb,iboun)
           if (kfl_fixno_nsa(1,ipoin) == 2) then        ! slip condition
              if (kfl_fixno_nsa(3,ipoin) == 0) then     ! do not include 201 condition
                 nifmark= nifmark+1
              end if
           else if (kfl_fixno_nsa(1,ipoin) == 4) then        ! those already changed in the first pass
              nifmark= nifmark+1
           end if           
        end do

        if ( nifmark == pnodb) then  
           
           ! three nodes to define a plane
           ipoin= lnodb(1,iboun)
           jpoin= lnodb(2,iboun)
           kpoin= lnodb(3,iboun)
           do idime= 1,ndime
              veta1(idime)= coord(idime,ipoin) - coord(idime,jpoin) 
              veta2(idime)= coord(idime,jpoin) - coord(idime,kpoin)         
           end do
           ! compute venoe, the vector normal to the boundary element
           call vecpro(veta2,veta1,venoe,ndime)
           rauxi= sqrt(venoe(1)*venoe(1) + venoe(2)*venoe(2) + venoe(3)*venoe(3)) 
           venoe(1)= venoe(1) / rauxi
           venoe(2)= venoe(2) / rauxi
           venoe(3)= venoe(3) / rauxi
           
           ! check angles with each exnor 
           do inodb=1,pnodb           
              ipoin= lnodb(inodb,iboun)
              ibopo= lpoty(ipoin)
              if (ipoin > 0) then
                 venor(1:ndime)= exnor(1:ndime,1,ibopo)
                 rauxi= 0.0_rp
                 do idime= 1,ndime
                    rauxi= rauxi + venor(idime)*venoe(idime)
                 end do
                 ! the threshold is arbitrary, it must detect acute edges only. it is set in nsa_reabcs.
                 rauxi= abs(rauxi)
                 if (rauxi < angle_tredg_nsa) then
                    kfl_fixno_nsa(1,ipoin) = 4_ip  ! first pass, put a fixno with no meaning (it is 2 digit)
                 end if                 
              end if
           end do
           
        end if

     end do boundaries_detect

     !
     ! Interchange fixno values, first pass
     !  
     call nsa_parall(9_ip)

     do ipoin = 1,npoin
        if (kfl_fixno_nsa(1,ipoin)==4) kfl_fixno_nsa(1,ipoin)= 0
     end do


     !
     ! Loop over fixnos to correct
     !
!     boundaries_release: do kboun = 1,nbouz(lzone(ID_NASTAL))
!        
!        iboun = lbouz(lzone(ID_NASTAL)) % l(kboun)
!        pblty = ltypb(iboun) 
!        pnodb = nnode(pblty)
!        ielem = lelbo(iboun)
!
!        do inodb=1,pnodb
!           ipoin= lnodb(inodb,iboun)
!           if (kfl_fixno_nsa(1,ipoin) == 11) then
!              kfl_fixno_nsa(1,ipoin) = 4
!           end if
!        end do
!
!     end do boundaries_release

  end if
  

end subroutine nsa_traedg

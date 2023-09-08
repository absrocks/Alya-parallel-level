!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_noredg.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Correct the 201 boundary condition for wall-symmetry edges (only 3D)
!> @details Correct the 201 boundary condition for wall-symmetry edges (only 3D)
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_noredg
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use mod_memchk
  use mod_memory

  use def_nastal
  implicit none

  integer(ip)  :: &
       iboun,kboun,pblty,pnodb,ielem,ifmark(mnodb),nifmark,kifmark,&
       ibopo,idime,ipoin,jpoin,kpoin,inodb
  integer(4)   :: istat
  real(rp)     :: venor(ndime),veta1(ndime),veta2(ndime),venoe(ndime),rauxi

  real(rp), pointer :: veta1_global(:,:) 
  real(rp), pointer :: veta2_global(:,:) 

  nullify(veta1_global)
  nullify(veta2_global)


  !
  ! Allocate temporary vectors
  !
  if (IMASTER) then
     allocate(veta1_global(1,1),stat=istat) 
     allocate(veta2_global(1,1),stat=istat) 
  else
     allocate(veta1_global(ndime,npoin),stat=istat)
     call memchk(0_ip,istat,mem_modul(1:2,modul),'VETA1_GLOBAL','nsa_noredg',veta1_global)
     allocate(veta2_global(ndime,npoin),stat=istat)
     call memchk(0_ip,istat,mem_modul(1:2,modul),'VETA2_GLOBAL','nsa_noredg',veta2_global)
  end if


  if( INOTMASTER ) then
     !
     ! Loop over boundaries to erase tangents
     !
     kifmark= 0
     boundaries_erase: do iboun = 1,nboun
        pblty = ltypb(iboun) 
        pnodb = nnode(pblty)
        ielem = lelbo(iboun)
        ifmark= 0
        nifmark= 0
        do inodb=1,pnodb
           ipoin= lnodb(inodb,iboun)
           if (kfl_fixno_nsa(1,ipoin) == 2) then
              if (kfl_fixno_nsa(ndime,ipoin) == 1) then
                 ! node with condition 201 identified
                 ifmark(inodb)= ipoin                         
                 nifmark= nifmark+1
              else
                 ! node with condition 200, so this element is not symetry with nodes 201,
                 ! but wall with nodes 201...
                 nifmark= -(10*mnodb)         
              end if
           end if
        end do

        if ( nifmark > 1) then  
           kifmark= 1                               ! at least one 201 detected
           ! boundary with more than 1 nodes 201
           do inodb=1,pnodb           
              ipoin= ifmark(inodb)
              if (ipoin > 0) then
                 ibopo= lpoty(ipoin)
                 ! erase both tangents
                 do idime= 1,ndime
                    exnor(idime,2,ibopo) = 0.0_rp
                    exnor(idime,3,ibopo) = 0.0_rp
                 end do
              end if
           end do

        end if

     end do boundaries_erase

     if (kifmark == 1) then
        !
        ! Loop over boundaries to correct triads
        !
        boundaries_correct: do iboun = 1,nboun
           pblty = ltypb(iboun) 
           pnodb = nnode(pblty)
           ielem = lelbo(iboun)
           ifmark= 0
           nifmark= 0
           do inodb=1,pnodb
              ipoin= lnodb(inodb,iboun)
              if (kfl_fixno_nsa(1,ipoin) == 2) then
                 if (kfl_fixno_nsa(ndime,ipoin) == 1) then
                    ! node with condition 201 identified
                    ifmark(inodb)= ipoin                         
                    nifmark= nifmark+1
                 else
                    ! node with condition 200, so this element is not symetry with nodes 201,
                    ! but wall with nodes 201...
                    nifmark= -(10*mnodb)         
                 end if
              end if
           end do
           
           if ( nifmark > 1) then
              
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
              
              do inodb=1,pnodb           
                 ipoin= ifmark(inodb)
                 if (ipoin > 0) then
                    ibopo= lpoty(ipoin)
                    venor(1:ndime)= exnor(1:ndime,1,ibopo)
                    
                    ! compute veta1, the new first tangent vector
                    call vecpro(venor,venoe,veta1,ndime)
                    rauxi= sqrt(veta1(1)*veta1(1) + veta1(2)*veta1(2) + veta1(3)*veta1(3)) 
                    veta1(1)= veta1(1) / rauxi
                    veta1(2)= veta1(2) / rauxi
                    veta1(3)= veta1(3) / rauxi
                    
                    ! compute veta2, the new second tangent vector
                    call vecpro(venor,veta1,veta2,ndime)
                    rauxi= sqrt(veta2(1)*veta2(1) + veta2(2)*veta2(2) + veta2(3)*veta2(3)) 
                    veta2(1)= veta2(1) / rauxi
                    veta2(2)= veta2(2) / rauxi
                    veta2(3)= veta2(3) / rauxi
                    
                    do idime= 1,ndime
                       exnor(idime,2,ibopo)= exnor(idime,2,ibopo) + veta1(idime)
                       exnor(idime,3,ibopo)= exnor(idime,3,ibopo) + veta2(idime)
                    end do
                    
                 end if
              end do

           end if

        end do boundaries_correct
                
        !
        ! Parallel interchange of the tangent vectors
        !
        !
        !   1. Copy boundary defined exnor tangents in globally defined unkno and rhsid 
        !
        if (ISLAVE) then
           do ipoin = 1,npoin        
              do idime= 1,ndime
                 veta1_global(idime,ipoin) = 0.0_rp
                 veta2_global(idime,ipoin) = 0.0_rp
              end do
              ibopo= lpoty(ipoin)
              if (ibopo > 0) then
                 do idime= 1,ndime
                    veta1_global(idime,ipoin) = exnor(idime,2,ibopo)
                    veta2_global(idime,ipoin) = exnor(idime,3,ibopo)
                 end do
              end if
           end do
        end if
        
     end if

  end if
  
  !
  !   2. Perform the parallel interchange
  !
  call rhsmod(ndime,veta1_global)
  call rhsmod(ndime,veta2_global)
  
  !
  !   3. Copy globally defined unkno and rhsid to boundray defined exnor tangents
  !
  if (ISLAVE) then        
     if (kifmark == 1) then
        do ipoin = 1,npoin        
           ibopo= lpoty(ipoin)
           if (ibopo > 0) then
              do idime= 1,ndime
                 exnor(idime,2,ibopo) = veta1_global(idime,ipoin) 
                 exnor(idime,3,ibopo) = veta2_global(idime,ipoin) 
              end do
           end if
        end do
     end if
     
  end if

  if (INOTMASTER) then
     !
     ! Normalize new tangents
     !
     do ibopo= 1,nbopo
        venor(2)= sqrt(exnor(1,2,ibopo)*exnor(1,2,ibopo) &
             + exnor(2,2,ibopo)*exnor(2,2,ibopo) + exnor(3,2,ibopo)*exnor(3,2,ibopo))        
        venor(3)= sqrt(exnor(1,3,ibopo)*exnor(1,3,ibopo) &
             + exnor(2,3,ibopo)*exnor(2,3,ibopo) + exnor(3,3,ibopo)*exnor(3,3,ibopo))        
        do idime= 1,ndime
           exnor(idime,2,ibopo)= exnor(idime,2,ibopo) / venor(2)
           exnor(idime,3,ibopo)= exnor(idime,3,ibopo) / venor(3)
        end do
     end do

     !
     ! Free memory
     !
     if (IMASTER) then
        deallocate(veta1_global,stat=istat) 
        deallocate(veta2_global,stat=istat) 
     else
        call memchk(2_ip,istat,mem_modul(1:2,modul),'VETA1_GLOBAL','nsa_noredg',veta1_global)
        deallocate(veta1_global,stat=istat)
        call memchk(2_ip,istat,mem_modul(1:2,modul),'VETA2_GLOBAL','nsa_noredg',veta2_global)
        deallocate(veta2_global,stat=istat)
     end if

  end if


end subroutine nsa_noredg

!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_inibcs.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Initialize the boundary conditions
!> @details Initialize the boundary conditions
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_inibcs
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use mod_memchk

  use def_nastal
  implicit none

  integer(ip)  :: ipoin,kpoin,ibopo,idofn

  if( INOTMASTER ) then
     !
     ! Allocate memory
     !
     call nsa_membcs(1_ip)
     call nsa_membcs(2_ip)

     !-------------------------------------------------------------
     !
     ! Node codes
     !
     !-------------------------------------------------------------

     if( kfl_icodn > 0 ) then               !if( kfl_icodn == 0 ) CODE CONDITIONS ARE ASSIGNED BUT THERE IS NO NODAL CODE
        if( kfl_conbc_nsa == 0 ) then       !kfl_conbc_nsa == 1 by default
           iffun      =  1
           kfl_funno  => kfl_funno_nsa
        else
           iffun      =  0
        end if
        ifbop     =  0
        ifloc     =  1
        kfl_fixrs => kfl_fixrs_nsa
        kfl_fixno => kfl_fixno_nsa
        bvess     => bvess_nsa(:,:,1)
        tncod     => tncod_nsa(1:)

        call reacod(10_ip)
  
     end if
     !-------------------------------------------------------------
     !
     ! Boundary codes
     !
     !-------------------------------------------------------------

     if( kfl_icodb > 0 ) then
        if( kfl_conbc_nsa == 0 ) then
           iffun      =  1
           kfl_funbo  => kfl_funbo_nsa
        end if
        kfl_fixbo => kfl_fixbo_nsa
        bvnat     => bvnat_nsa
        tbcod     => tbcod_nsa(1:)
        !tncod     => momod(ID_NASTIN) % tncod(1:)   
        call reacod(20_ip)
     end if
     !
     !
     ! Final corrections
     !
     !
     do ipoin = 1,npoin
        !
        ! Correct errors of reacod assignement on interior nodes
        !
        ibopo= lpoty(ipoin)
        !I think the following if should go (Fer)
        if (ibopo == 0) then                     ! ibopo=0 means this is an inner node
           do idofn= 1,ndofn_nsa
!              if (kfl_fixno_nsa(idofn,ipoin) == -1) kfl_fixno_nsa(idofn,ipoin) = 0
              kfl_fixno_nsa(idofn,ipoin) = 0
           end do

        else if (ibopo > 0) then                 ! ibopo>0 means this is boundary node

           if(kfl_fixno_nsa(1,ipoin)==2 .and. kfl_fixrs_nsa(ibopo)==0) &
                kfl_fixrs_nsa(ibopo)=-1
           if(kfl_fixno_nsa(1,ipoin)==5 .and. kfl_fixrs_nsa(ibopo)==0) &
                kfl_fixrs_nsa(ibopo)=-1
           if(kfl_fixno_nsa(1,ipoin)==9 .and. kfl_fixrs_nsa(ibopo)==0) &
                kfl_fixrs_nsa(ibopo)=-1
        end if


     end do

  end if

  !
  ! Correct the famous "201 condition" for 3D problems
  !     
  if (ndime > 2) call nsa_noredg

  !
  ! Correct boundary conditions in trailing edges
  !     
  if (ndime > 2) then
     if (kfl_tredg_nsa == 1) call nsa_traedg
  end if

  if( INOTMASTER ) then
     
     !
     ! Non-constant boundary conditions: store references in bvess_nsa(:,:,2) 
     !     
     if( kfl_conbc_nsa == 0 ) then
        do ipoin = 1,npoin
           bvess_nsa(1:ndofn_nsa,ipoin,2) = bvess_nsa(1:ndofn_nsa,ipoin,1)
        end do
        !        do iboun = 1,nboun
        !           bvnat_nsa(1,iboun,2) = bvnat_nsa(1,iboun,1)
        !        end do
     end if

  end if

  
end subroutine nsa_inibcs

!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_bouimp.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Boundary conditions for the implicit case
!> @details Boundary conditions for the implicit case
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_bouimp(ktask)
!-----------------------------------------------------------------------
!****f* Nastal/nsa_bouimp
! NAME 
!    nsa_bouimp
! DESCRIPTION
!    This routine computes boundary conditions for the implicit case:
!    
!    itask=1 converts physical to implicit (on the unknowns) b.c.
!    itask=2 converts back implicit to physical (on the unknowns) b.c.
!    
!    kfl_fixno fixes u_i, rho, T
!    the following rules apply (suppose 3D):
!    
!    supersonic inflow
!      11111 fixes u, rho, T, then it will fix U, rho, E
!    supersonic outflow
!      all free
!    subsonic inflow
!      11101 fixes u and T, then it will TO BE SET
!    subsonic outflow
!      00010 fixes rho, then it will fix rho
!    slip condition
!      20000 fixes u_n to zero, then it will fix U_n to zero 
!    non-slip condition
!      11101 fixes u and T, then it will fix TO BE SET
!    others
!      TO BE SET
!    
!    
! USED BY
!    nsa_upcons
!***
!-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_parame

  use      def_nastal

  integer(ip) ,intent(in)      :: ktask   !< Conversion direction
  integer(ip)      :: ipoin,kpoin,idime,idofn,kvelo,kdens,ktemp,ktota
  real(rp)         :: vsqua

  if (ktask == 1) then
     do ipoin = 1,npoin
        kvelo= 0
        kdens= 0
        ktemp= 0
        ktota= 0
        do idime=1,ndime
           if (kfl_fixno_nsa(idime,ipoin)==1)then
              kvelo= kvelo+1 
           else if (kfl_fixno_nsa(idime,ipoin)==2)then
              kvelo= -1                                   ! usually the first dimension (normal)
           end if           
        end do
        if (kfl_fixno_nsa(ndime+1,ipoin)==1) kdens= 1
        if (kfl_fixno_nsa(ndime+2,ipoin)==1) ktemp= 1
        ktota= kvelo+kdens+ktemp

        if (ktota == ndofn_nsa) then
           ! supersonic inflow
           vsqua= 0.0_rp
           do idime= 1,ndime
              vsqua= vsqua + bvess_nsa(idime,ipoin,1)*bvess_nsa(idime,ipoin,1)
              bvess_nsa(idime,ipoin,1)= &
                   bvess_nsa(ndime+1,ipoin,1)*bvess_nsa(idime,ipoin,1)  ! u to U
           end do
           bvess_nsa(ndime+2,ipoin,1)= &
                bvess_nsa(ndime+1,ipoin,1) &
                * (cvcoe_nsa * bvess_nsa(ndime+2,ipoin,1)+ vsqua * 0.5_rp) ! T to E
           
        else

           !.... TO BE PROGRAMMED

        end if

     end do

  else if (ktask == 2) then
     do ipoin = 1,npoin
        kvelo= 0
        kdens= 0
        ktemp= 0
        ktota= 0
        do idime=1,ndime
           if (kfl_fixno_nsa(idime,ipoin)==1)then
              kvelo= kvelo+1 
           else if (kfl_fixno_nsa(idime,ipoin)==2)then
              kvelo= -1                                   ! usually the first dimension (normal)
           end if           
        end do
        if (kfl_fixno_nsa(ndime+1,ipoin)==1) kdens= 1
        if (kfl_fixno_nsa(ndime+2,ipoin)==1) ktemp= 1
        ktota= kvelo+kdens+ktemp

        if (ktota == ndofn_nsa) then
           ! supersonic inflow
           vsqua= 0.0_rp
           do idime= 1,ndime
              bvess_nsa(idime,ipoin,1)= bvess_nsa(idime,ipoin,1) / bvess_nsa(ndime+1,ipoin,1)  ! U to u
              vsqua= vsqua + bvess_nsa(idime,ipoin,1)*bvess_nsa(idime,ipoin,1)
           end do
           bvess_nsa(ndime+2,ipoin,1)= (bvess_nsa(ndime+2,ipoin,1)/bvess_nsa(ndime+1,ipoin,1) - vsqua * 0.5_rp) / cvcoe_nsa ! E to T
           
        else

           !.... TO BE PROGRAMMED

        end if

     end do


  end if


end subroutine nsa_bouimp

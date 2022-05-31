subroutine nsa_turbul(dwall,velmo,xvisc,gvelo,dvolu,xnutu,xunkn)
!-----------------------------------------------------------------------
!****f* Nastal/nsa_turbul
! NAME 
!    nsa_turbul
! DESCRIPTION
!    This routine computes the turbulent viscosity from the gradients of the velocity field 
!    at gauss points
!    Turbulence models available:
!       --> SMAGORINSKY (kfl_cotur_nsa == -1)
!           --> A Van Driest damping function is also available (Piomelli et al 1987) 
!       --> WALE        (kfl_cotur_nsa == -2)
!       --> SIGMA       (kfl_cotur_nsa == -3)
!
! USED BY
!    nsa_gauvalxy
!***
!-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_parame
  use      def_nastal

  implicit none
  real(rp),      intent(in) :: gvelo(ndime,ndime),dvolu,xunkn(ndofn_nsa),dwall,velmo,&
                               xvisc
  real(rp),      intent(out):: xnutu
  integer(ip)               :: ipoin,idime,jpoin,jdime,kdime,idumi
  real(rp)                  :: xmile,seci4,g2_ij(ndime,ndime), &
                               g2_kk,seci5,sd_ij,denom,G__ij(ndime,ndime), &
                               G_val(3),vdumm(3,3),sigma(3),ynorm,fvaDr,Smagc

  ! Compute turbulent viscosity xnutu   
    
  xmile = 0.0_rp
  xnutu = 0.0_rp

  xmile = dvolu**0.3333333_rp

  if (kfl_cotur_nsa == -1) then           
  !
  ! SMAGORINSKY 
  ! -------------------------------------------------------------
     seci4 = 0.0_rp
     do idime = 1,ndime                     ! 2 S_ij : S_ij
        do jdime = 1,ndime         
           seci4 = seci4 + gvelo(idime,jdime) * (gvelo(idime,jdime) + gvelo(jdime,idime))
        end do
     end do

     ! Van Driest dumping function
     !
     if (kfl_dampf_nsa == 1 ) then
        ynorm = 0.04_rp * dwall * (velmo / ndime) * xunkn(ndime+1)  / xvisc  
        fvaDr = 1.0_rp - exp(- (ynorm * ynorm * ynorm))
        Smagc = turbu_nsa*xmile*xmile*fvaDr
     else
        Smagc = turbu_nsa*xmile*xmile
     endif

     xnutu = Smagc*sqrt(seci4)

  else if  (kfl_cotur_nsa == -2) then
  !
  ! WALE  - Wall-adapting local eddy-viscosity - Nicoud-Ducros 99
  ! -------------------------------------------------------------
  !
     ! h calculated as vol**(1/ndime)  
     !
     if (ndime == 2) then
        call runend('nsa_lawvis: WALE only ready in 3d - 2d LES does not make sense')
     endif 
          
     ! Computation of the square of the velocity gradient g2_ij

     g2_ij = 0.0_rp
     g2_kk = 0.0_rp

     do idime = 1,ndime
        do jdime = 1,ndime
           do kdime = 1,ndime
              g2_ij(idime,jdime) = g2_ij(idime,jdime) + gvelo(idime,kdime)*gvelo(kdime,jdime)
           end do
        end do 
     end do

     g2_kk = g2_ij(1,1) + g2_ij(2,2) + g2_ij(3,3)

     seci4 = 0.0_rp
     seci5 = 0.0_rp
     sd_ij = 0.0_rp
     denom = 0.0_rp

     do idime = 1,ndime
        do jdime = 1,ndime
           seci4 = seci4 + 0.5_rp*gvelo(idime,jdime) * &          ! S_ij : S_ij
                   (gvelo(idime,jdime) + gvelo(jdime,idime)) 
           sd_ij = 0.5_rp * ( g2_ij(idime,jdime) + g2_ij(jdime,idime) ) 
           if (idime == jdime) then
              sd_ij = sd_ij - 0.3333333_rp*g2_kk    ! -1/3 delta_ij  * gkk**2
           endif
           seci5 = seci5 + sd_ij*sd_ij     ! Sd_ij:Sd_ij
        end do
     end do

     denom = max ( zeror , seci4**2.5_rp + seci5**1.25_rp )

     xnutu = turbu_nsa*xmile*xmile*seci5**1.5_rp/denom

  else if( kfl_cotur_nsa == -3 ) then     
  !
  ! Sigma model  - Baya Toda, Nicoud 2010
  ! -------------------------------------------------------------
  ! h calculated as vol**(1/ndime)  
  !
     G__ij = 0.0_rp                      ! G = g^T*g = g_ki * g_kj 
     G_val = 0.0_rp                      ! Note gvelo(kdime,idime) is g_ik = du_i/dx_k from  Nicoud

     do idime = 1,ndime
        do jdime = 1,ndime
           do kdime = 1,ndime
              G__ij(idime,jdime) = G__ij(idime,jdime) + gvelo(idime,kdime)*gvelo(jdime,kdime)
           end do
        end do 
     end do
        
     call spcdec(G__ij,G_val,vdumm,idumi,0_ip,'NSA_TURBUL')   ! eigenvalues sorted in descending order
     !
     ! sigma is in descending order
     !
     if (abs(G_val(3)) <= 2.3d-10  )then
        G_val(3) = 0.0_rp
     elseif  (G_val(3) <= -2.3d-10  )then
        write(*,*)'Sigma values',G_val(1),G_val(2),G_val(3)
        write(*,*)
        write(*,*)'Tensor G__ij',G__ij
        call runend('nsa_lawvis: negative eigenvalue in sigma-model')
     endif

     sigma(1) = sqrt(G_val(1))
     sigma(2) = sqrt(G_val(2))
     sigma(3) = sqrt(G_val(3))

     xnutu = turbu_nsa*xmile*xmile

     if ( sigma(1)*sigma(1) > 10.0_rp * zeror ) then     ! Avoid divide by zero
        xnutu = xnutu * sigma(3) * ( sigma(1) - sigma(2) ) * ( sigma(2) - sigma(3) ) / ( sigma(1) * sigma(1) )
     else
        xnutu = 0.0_rp
     end if

  end if
     
end subroutine nsa_turbul

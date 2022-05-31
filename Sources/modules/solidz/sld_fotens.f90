!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_fotens.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Compute some tensors for postprocess
!> @details Compute some tensors for postprocess
!> @}
!-----------------------------------------------------------------------

subroutine sld_fotens()

  use def_kintyp,   only : ip, rp
  use def_elmtyp
  use def_master
  use def_domain
  use def_solidz
  implicit none
  integer(ip)             :: ipoin,ivoig,ibopo,idumm,imeth, i, j
  real(rp)                :: dummr,nx,ny,vdumm(3,3),G_ij(3,3),G_val(3)
  real(rp)                :: E_ij(3,3),B_ij(3,3),B_val(3),iden(3,3)

  !
  ! Initialize
  !
  do ipoin = 1,npoin
     do ivoig = 1,nvoig_sld
        caust_sld(ivoig,ipoin) = 0.0_rp
        green_sld(ivoig,ipoin) = 0.0_rp
        lepsi_sld(ivoig,ipoin) = 0.0_rp
  !      epsee_sld(ivoig,ipoin) = 0.0_rp
     end do
  end do

  grlst_sld(:,:) = 0.0_rp
  iden(:,:) = 0.0_rp
  do i=1,3
    iden(i,i) = 1.0_rp
  end do

  imeth = 1
  if( imeth == 2 ) call opeclo(1_ip)
  call sld_elmope(3_ip)
 ! call sld_elmope(15_ip)
  if( imeth == 2 ) call opeclo(2_ip)

  call rhsmod(nvoig_sld,caust_sld)
  call rhsmod(nvoig_sld,green_sld)
  call rhsmod(nvoig_sld,lepsi_sld)
  if (kfl_plast_sld == 1_ip) call rhsmod(nvoig_sld,epsee_sld)

  do ipoin = 1,npoin
     dummr = 1.0_rp/vmass(ipoin)
     do ivoig = 1,nvoig_sld
        caust_sld(ivoig,ipoin) = dummr * caust_sld(ivoig,ipoin)
        green_sld(ivoig,ipoin) = dummr * green_sld(ivoig,ipoin)
        lepsi_sld(ivoig,ipoin) = dummr * lepsi_sld(ivoig,ipoin)
        if (kfl_plast_sld == 1_ip) epsee_sld(ivoig,ipoin) = dummr * epsee_sld(ivoig,ipoin)
     end do

     ibopo= lpoty(ipoin)
     caunn_sld(ipoin)= 0.0_rp
     if (ibopo > 0) then
        if (ndime == 2)  then
           nx = exnor(1,1,ibopo)
           ny = exnor(2,1,ibopo)
           caunn_sld(ipoin)= &
                nx * caust_sld(1,ipoin) * nx + &
                ny * caust_sld(2,ipoin) * ny + &
                2.0_rp * nx * caust_sld(3,ipoin) * ny
        else if (ndime == 3) then

        end if
     end if
  end do
  !
  ! Combined Stresses
  !
  if( ndime == 3_ip ) then

     do ipoin = 1,npoin
        !
        ! Von Mises Stress
        !
        seqvm_sld(ipoin) = &
             (sqrt(0.5_rp*((caust_sld(1,ipoin)-caust_sld(2,ipoin))**2.0_rp+&
             (caust_sld(2,ipoin)-caust_sld(3,ipoin))**2.0_rp+&
             (caust_sld(3,ipoin)-caust_sld(1,ipoin))**2.0_rp+&
             6.0_rp*(caust_sld(4,ipoin)**2.0_rp+caust_sld(5,ipoin)**2.0_rp+caust_sld(6,ipoin)**2.0_rp))))

        if (kfl_rotei_sld == 1) call sld_troloc(0_ip,G_ij)     !compute roloc_sld tensor

        !G_ij(1,1)= caust_sld(1,ipoin)
        !G_ij(2,2)= caust_sld(2,ipoin)
        !G_ij(3,3)= caust_sld(3,ipoin)
        !G_ij(1,2)= caust_sld(4,ipoin)
        !G_ij(1,3)= caust_sld(5,ipoin)
        !G_ij(2,3)= caust_sld(6,ipoin)
        !G_ij(2,1)= G_ij(1,2)
        !G_ij(3,1)= G_ij(1,3)
        !G_ij(3,2)= G_ij(2,3)

        do i = 1, ndime
           do j = 1, ndime
              ivoig = nvgij_inv_sld(i, j)
              G_ij(i, j) = caust_sld(ivoig,ipoin)
              E_ij(i, j) = green_sld(ivoig,ipoin)
           end do
        end do

        if (kfl_rotei_sld == 1) call sld_troloc(ipoin,G_ij)     !rotate and correct tensor
        if (kfl_rotei_sld == 1) call sld_troloc(ipoin,E_ij)     !rotate and correct tensor

        ! eigenvalues sorted in ascending order
        call spcdec(G_ij,G_val,vdumm,idumm,0_ip,'SLD_FOTENS, COMPUTE CAUCHY EIGENVALUES FOR SIGEI')

        G_val(1) = abs(G_val(1))
        G_val(3) = abs(G_val(3))
        if (G_val(3) .gt. G_val(1)) G_val(1)= G_val(3)

        sigei_sld(ipoin) = G_val(1)

       !
       ! Principal stretches
       !       
       B_ij = transpose(2.0_ip*E_ij + iden) 
       
       ! eigenvalues sorted in ascending order
       call spcdec(B_ij,B_val,vdumm,idumm,0_ip,'SLD_FOTENS, COMPUTE GREEN-LAGRANGE EIGENVALUES FOR EPRIN')
        
        B_val(1) = abs(B_val(1))
        B_val(3) = abs(B_val(3))
        if (B_val(3) .gt. B_val(1)) B_val(1)= B_val(3)
 
       eprin_sld(ipoin) = B_val(1)
     end do

  else if ( ndime == 2_ip ) then
     !
     ! Von Mises Stress (General Plane Stress assumption)
     !
     seqvm_sld(1:npoin) = sqrt(caust_sld(1,1:npoin)**2.0_rp + caust_sld(2,1:npoin)**2.0_rp &
          - caust_sld(1,1:npoin)*caust_sld(2,1:npoin) + 3.0_rp*caust_sld(3,1:npoin)**2.0_rp)

     sigei_sld(:) = 0.0_rp ! not coded

  end if

end subroutine sld_fotens

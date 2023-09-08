subroutine sld_stress_model_104(pgaus,pmate,gpgdi,gpstr,gpdet,flagt,gpdds,gpmof)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_elmcla
  ! NAME
  !    sld_stress_model_103
  ! DESCRIPTION
  !    Compressible Mooney-Rivlin stress model
  !    Xiao and Belytschko's formulation
  !    Compute second Piola-Kirchoff stress tensor S_{IJ}
  !
  !    GPGDI ... Deformation tensor ...................... F = grad(phi)
  !    GPCAU ... Right Cauchy-Green deformation tensor ... C = F^t x F
  !    GPCIN ... Inverce of Right Caughy-Green tensor .... C^-1
  !    GPSTR ... 2nd P-K Stress tensor ........................... S
  !    GPPIO ... 1st Piola-Kirchhoff stress tensor ....... P = F.S
  !    GPENE ... Stored energy function .................. W
  !    FLAGT ... Flag to activate GPDDS (when implicit)
  ! USES
  ! USED BY
  !    sld_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode
  use def_solidz, only       :  lawco_sld,parco_sld,densi_sld,velas_sld
  use def_master, only       :  ittim
  implicit none
  integer(ip), intent(in)    :: pgaus,pmate,flagt
  real(rp),    intent(in)    :: gpmof(pgaus)
  real(rp),    intent(in)    :: gpgdi(ndime,ndime,pgaus),gpdet(pgaus)
  real(rp),    intent(out)   :: gpstr(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: gpdds(ndime,ndime,ndime,ndime,pgaus)
  real(rp)                   :: gpcau(ndime,ndime,pgaus), gpene(pgaus)
  real(rp)                   :: gpcin(ndime,ndime,pgaus),tkron(ndime,ndime)
  integer(ip)                :: igaus,idime,jdime,kdime,i,j,k,l
  real(rp)                   :: lambda0,mu0,logj,bidon,i1,i2
  real(rp)                   :: gpcal(ndime,ndime,pgaus), dummr, gpcau2(ndime,ndime,pgaus), gppre(pgaus), Iden(ndime,ndime)
  real(rp)                   :: cstc1,cstc2,cstd1,cstm1,cstm2


  ! Mooney Rivlin's law (not implemented here)
  ! W(C) = C_1*(I_1-3) + C_2*(I_2-3) + M_1*(exp(M_2(I_1-3))-1) 

  cstd1 = parco_sld(2,pmate)  ! also called D_1, here is lambda
  cstc1 = parco_sld(3,pmate)  ! also called C_1
  cstc2 = parco_sld(4,pmate)  ! C_2=0 is the NeoHook material 
  cstm1 = parco_sld(5,pmate)  
  cstm2 = parco_sld(6,pmate)  

  gpstr = 0.0_rp

  do idime = 1,ndime
     do jdime = 1,ndime
        if (idime==jdime) then
           Iden(idime,jdime) = 1.0_rp
        else
           Iden(idime,jdime) = 0.0_rp
        end if
     enddo
  enddo

  ! Kronecker delta
  tkron = 0.0_rp
  do idime = 1, ndime
      tkron(idime,idime) = 1.0_rp
  end do

  !GPCAL: Left Cauchy tensor b = F F^T
  !
  do igaus=1,pgaus
     do idime = 1,ndime
        do jdime = 1,ndime
           gpcal(idime,jdime,igaus) = 0.0_rp
           do kdime = 1,ndime
              gpcal(idime,jdime,igaus) = gpcal(idime,jdime,igaus) + gpgdi(idime,kdime,igaus)*gpgdi(jdime,kdime,igaus)
           end do
        end do
     end do
  end do

  !Alternative computation
  !do igaus=1,pgaus
  !      gpcal(:,:,igaus) = matmul(gpgdi2(:,:,igaus),transpose(gpgdi2(:,:,igaus)))
  !      traceb = 0.0_rp
  !      do idime=1,ndime
  !         traceb = traceb + gpcal(idime,idime,igaus)
  !      end do
  !enddo



  do igaus=1,pgaus

     ! compute C^2 for the second invariant
     gpcau2 = 0.0_rp
     gpcau2(:,:,igaus) = matmul(gpcau(:,:,igaus),gpcau(:,:,igaus))

     call invmtx(gpcau(:,:,igaus),gpcin(:,:,igaus),bidon,ndime)

     ! * * * *
     !INVARIANTS
     ! * * * *

     !I_1=Trace(C)
     i1=0.0_rp
     i2=0.0_rp
     dummr=0.0_rp
     do idime = 1,ndime
        i1 = i1 + gpcau(idime,idime,pgaus)
        dummr = dummr + gpcau2(idime,idime,pgaus)
     end do
     i2 = 0.5_rp*(i1*i1 - dummr);
  

     ! Energy density function
     gpene(igaus) = cstc1*(i1-3.0_rp) + cstc2*(i2-3.0_rp) & 
                    + cstm1*(exp(cstm2*(i1-3.0_rp))-1.0_rp)
                    !& - 0.25_rp*cstd1*((log(gpdet(igaus)))**2.0_rp))


     ! Belystschko
     do idime=1,ndime
        do jdime=1,ndime
           gpstr(idime,jdime,igaus) = (cstc1 + cstc2*i1)*Iden(idime,jdime) - cstc2*gpcau(idime,jdime,igaus) &
                                     + 2.0_rp*cstm1*cstm2*exp(cstm2*(i1-3.0_rp))* Iden(idime,jdime)
                                      ! +(2.0_rp/cstd1)*(gpdet(igaus)-1)*gpdet(igaus)*gpcin(idime,jdime,igaus)
        end do
     end do

     if (flagt == 1_ip) then

       do i=1,ndime
         do j=1,ndime
           do k=1,ndime
             do l=1,ndime 
               gpdds(i,j,k,l,igaus) =  2.0_rp*cstc2*(tkron(i,k)*tkron(j,l)+tkron(i,l)*tkron(j,k)) &
                    + 2.0_rp*cstm1*cstm2*cstm2*exp(cstm2*(i1-3.0_rp))*(tkron(i,k)*tkron(j,l)+tkron(i,l)*tkron(j,k))
                    !+ 2.0_rp*cstd1*gpdet(igaus)*((2.0_rp*gpdet(igaus)-1.0_rp)&
                    !*(gpcin(i,k,igaus)*gpcin(j,l,igaus)+gpcin(i,l,igaus)*gpcin(j,k,igaus))&
                    !- 0.5_rp*(gpdet(igaus)-1.0_rp)*(tkron(i,k)*tkron(j,l)+tkron(i,l)*tkron(j,k)))
             enddo
           enddo
         enddo
       enddo


     endif

  enddo

end subroutine sld_stress_model_104

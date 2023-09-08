!---------------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_stress_model_100.f90
!> @author  
!> @date    
!>
!> @brief   Elastic isotropic material
!> @details GPGDI ... Deformation tensor .............................. F = grad(phi))
!>          GPPRE ... Green-Lagrange Strain ........................... E = 1/2(F^t.F-delta I)
!>          GPSTR ... 2nd P-K Stress tensor ........................... S
!>          GPDDS ... Tangent moduli in terms of 2nd P-K stress ....... dS/dE
!>          FLAGT ... Flag to activate GPDDS (when implicit)
!>
!> @}
!---------------------------------------------------------------------------------
subroutine sld_stress_model_100(pgaus, pmate, gpgdi, gpstr, gpcau, gpdet, gptlo, gpigd, ielem, elcod, flagt, gpdds,&
                                gpigd_eps, gpgdi_eps, gpdet_eps, gpmof)

  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode
  use def_solidz, only       :  lawco_sld,parco_sld,densi_sld,velas_sld,kfl_plane_sld,kfl_ninex_sld,epsex_sld
  use def_master, only       :  kfl_modul,cutim,ittim,postp,ittyp,ITASK_ENDRUN
  use def_master, only       :  lun_livei,coupling


  implicit none
  integer(ip), intent(in)    :: pgaus,pmate,ielem,flagt
  real(rp),    intent(in)    :: gpgdi(ndime,ndime,pgaus),gpmof(pgaus)
  real(rp),    intent(out)   :: gpstr(ndime,ndime,pgaus,2)
  real(rp),    intent(out)   :: gpdds(ndime,ndime,ndime,ndime,pgaus)
  real(rp)                   :: gpgre(ndime,ndime,pgaus,2)
  integer(ip)                :: igaus,i,j,k,l,idime,ininex
  real(rp)                   :: lame1,lame2,young,poiss,trgpg

  real(rp)                   :: i4f,i4n,i4s,lamda(3),&
                                scasta,scastv,scastr,&
                                nfibe0(ndime,pgaus),norma0(ndime,pgaus),nshet0(ndime,pgaus),&
                                nfibe(ndime,pgaus),nfibt(ndime,pgaus),&
                                nshet(ndime,pgaus),nshtt(ndime,pgaus),&
                                norma(ndime,pgaus),normt(ndime,pgaus),tkron(ndime,ndime)

  real(rp)                   :: castr_active(ndime,ndime),gpstr_active(ndime,ndime)
  real(rp),    intent(in)    :: gpcau(ndime,ndime,pgaus),gpdet(pgaus),gptlo(pgaus),gpigd(ndime,ndime,pgaus)

  real(rp)                   :: trace(pgaus,2),gphyd(ndime,ndime,pgaus),gpdev(ndime,ndime,pgaus),plane,dummr
  real(rp),    intent(in)    :: elcod(ndime,mnode)
  
  real(rp),    intent(in)        :: gpgdi_eps(ndime,ndime,pgaus)           ! Displacement Gradient F for the perturbed state
  real(rp),    intent(in)        :: gpigd_eps(ndime,ndime,pgaus)           ! Displacement Gradient Inverse F^{-1} for the perturbed state
  real(rp),    intent(in)        :: gpdet_eps(pgaus)                       ! 

  real(rp)   :: gpgdi_aux(ndime,ndime,pgaus)           

  logical                    :: debugging
  
  debugging = .false.

  !Lame Constants:

  young = parco_sld(1,pmate)
  poiss = parco_sld(2,pmate)
  lame1 = parco_sld(3,pmate) !lambda
  lame2 = parco_sld(4,pmate) !mu

  if (kfl_plane_sld == 1_ip .and. ndime == 3_ip) then
    call runend('Trying to run a plane stress case in 3D. That is not possible.')
  end if 

  if (kfl_plane_sld == 1_ip) then ! 2D plane stress assumption
    plane = (1.0_rp-2.0_rp*poiss)/(1.0_rp-poiss)
  else                            ! 2D plane strain assumption 
    plane = 1_ip
  end if
  
  gpstr = 0.0_rp
  gpgre = 0.0_rp

  gphyd = 0.0_rp
  gpdev = 0.0_rp
  !
  ! Kronecker delta
  !
  tkron = 0.0_rp
  tkron(    1,    1) = 1.0_rp
  tkron(    2,    2) = 1.0_rp
  tkron(ndime,ndime) = 1.0_rp

  !
  ! Compute the tangent moduli if required (i.e., implicit scheme)
  ! Cijkl: Elastic moduli tensor
  !
  if( flagt == 1_ip ) then

     
        ! dSdE_{ijkl} = lame1*delta^{-1}_{ij}*delta^{-1}_{kl}
        !               + lame2*[delta^{-1}_{ik}*delta^{-1}_{jl}+
        !                        delta^{-1}_{il}*delta^{-1}_{jk}]
        !forall(igaus=1:pgaus, i=1:ndime, j=1:ndime, k=1:ndime, l=1:ndime)
        !   gpdds(i,j,k,l,igaus)= &
        !                lame1*tkron(i,j)*tkron(k,l)+ &
        !                lame2*((tkron(i,k)*tkron(j,l))+(tkron(i,l)*tkron(j,k)))
        !end forall
        !
        ! Cijkl = lambda * dij * dkl  +  mu * ( dik * djl + dil * djk )
        ! C     = lambda I x I + 2 mu I
        !
        ! dSdE_{ijkl} = lame1*plane_stress_factor*delta^{-1}_{ij}*delta^{-1}_{kl}
        !               + lame2*[delta^{-1}_{ik}*delta^{-1}_{jl}+
        !                        delta^{-1}_{il}*delta^{-1}_{jk}]
        forall(igaus=1:pgaus, i=1:ndime, j=1:ndime, k=1:ndime, l=1:ndime)
           ! when no modulating fields are present, gpmof is 1.0
           gpdds(i,j,k,l,igaus)= &
                        lame1*gpmof(igaus)*plane*tkron(i,j)*tkron(k,l)+ &
                        lame2*((tkron(i,k)*tkron(j,l))+(tkron(i,l)*tkron(j,k)))
        end forall  
     
  end if

  gpgdi_aux = gpgdi
  do ininex = 1, kfl_ninex_sld + 1

     do igaus = 1,pgaus
        !
        ! E_ij = 0.5(F_ki*F_kj-Kro_ij) (Holzapfel eq. 2.69)
        !
        do i = 1,ndime
           do j = 1,ndime
              if (ininex == 2) then
                 gpgdi_aux(i,j,igaus) = gpgdi(i,j,igaus) + 0.5_rp * epsex_sld * gpgdi(i,j,igaus)
              end if
              dummr = 0.0_rp
              do k = 1,ndime
                 dummr = dummr + gpgdi_aux(k,i,igaus) * gpgdi_aux(k,j,igaus)
              end do
              gpgre(i,j,igaus,ininex) = 0.5_rp*( dummr - tkron(i,j) )
              if (ininex == 2) then
                 gpgdi_aux(i,j,igaus) = gpgdi(i,j,igaus)
              end if
           end do
        end do
        !
        ! trace E
        !
        trace(igaus,ininex) = 0.0_rp
        do idime = 1,ndime
           trace(igaus,ininex) = trace(igaus,ininex) + gpgre(idime,idime,igaus,ininex)
        end do
     end do
     
     !
     ! S directly with fourth order tensor (Belytch p. 228)
     ! Sij = Cijkl Ekl = lambda * trace(E) delta_ij + 2 mu * Eij
     !
     do igaus = 1,pgaus
        ! when no modulating fields are present, gpmof is 1.0
        gpstr(:,:,igaus,ininex) = &
             lame1 * gpmof(igaus) * plane * trace(igaus,ininex) * tkron + 2.0_rp * lame2 * gpgre(:,:,igaus,ininex)
     end do
     
  end do

end subroutine sld_stress_model_100


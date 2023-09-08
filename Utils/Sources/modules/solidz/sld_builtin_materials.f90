!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_fe2.f90
!> @author  Guido Giuntoli
!> @date    July, 2018
!>
!> @brief   Interface for calling built-in subroutines
!> @details
!>
! INPUT
!    PGAUS ... Number of Gauss points
!    PMATE ... Material number
!    TEMP0 ... Previous temperature ............................ temp(n)
!    TEMP1 ... Updated temperature ............................. temp(n+1)
!    GPGI0 ... Previous deformation gradient tensor ............ F(n)
!    GPGDI ... Updated deformation gradient tensor ............. F(n+1)
!    GPIGD ... Inverse of updated deformation gradient tensor .. F^{-1}(n+1)
!    GPCAU ... Updated right Cauchy-Green tensor ............... C(n+1)
!    GPDET ... Updated jacobian ................................ J = det(F(n+1))
!    GPRAT ... Rate of Deformation tensor ...................... Fdot = grad(phidot)
!    FLAGT ... Flag for tangent moduli calculations ............ 1 if yes; 0 if no
!    FLAGS ... Flag for strain gradient calculations ........... 1 if yes; 0 if no
! INPUT/OUTPUT
!    GPPIO ... 1st Piola-Kirchhoff stress tensor at t(n)/t(n+1)  P(n)/P(n+1)
!    GPPIO_EPS ... Perturbed 1st Piola-Kirchhoff stress tensor at t(n)/t(n+1)  P(n)/P(n+1)
!    GPSTR ... 2nd Piola-Kirchhoff stress tensor at t(n)/t(n+1)  S(n)/S(n+1)
!    GPENE ... Stored energy function .......................... W
!    GPGRE ... Green-Lagrange         .......................... E
! OUTPUT
!    GPTMO ... Tangent moduli at t(n+1) ........................ dP/dF(n+1)
!    GGSGO ... Deformation gradient gradient tensor ............ dF/dX(n+1)
!> @}
!-----------------------------------------------------------------------

subroutine sld_builtin_materials(pgaus,pmate,gpvol,temp0,temp1,&
     gpgi0,gpgdi,gpigd,gpcau,gpdet,gprat,gppio,gpstr,&
     gpene,flagt,gptmo,flags,ggsgo,ielem,gptlo,nfibe,elcod,pnode,lnods,gpsha,gpdds,&
     gpigd_eps,gpgdi_eps,gpdet_eps,gppio_eps,&
     gpmof)

  use def_kintyp,               only:  ip,rp
  use def_domain,               only:  ndime,mnode
  use def_solidz
  use mod_sld_stress_model_152, only: sld_stress_model_152
  use mod_sld_stress_model_154, only: sld_stress_model_154
  use mod_sld_plastic_model,    only: sld_plastic_model_biso
  use mod_sld_fe2
  use def_material

  implicit none
  integer(ip), intent(in)    :: pgaus,pmate,ielem,pnode,lnods(pnode)
  real(rp),    intent(in)    :: temp0(pgaus)
  real(rp),    intent(in)    :: temp1(pgaus)
  real(rp),    intent(in)    :: gpgi0(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpgdi(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpcau(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpdet(pgaus),gpmof(pgaus)
  real(rp),    intent(in)    :: gprat(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus)
  integer(ip), intent(in)    :: flagt,flags
  real(rp),    intent(in)    :: gpigd(ndime,ndime,pgaus)
  real(rp),    intent(inout) :: gppio(ndime,ndime,pgaus)
  real(rp),    intent(inout) :: gppio_eps(ndime,ndime,pgaus)
  !   real(rp),    intent(out)   :: gpgre(ndime,ndime,pgaus)
  real(rp),    intent(inout) :: gpstr(ndime,ndime,pgaus,2)
  !   real(rp),    intent(inout) :: gplep(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: nfibe(ndime,pgaus)
  real(rp),    intent(inout) :: gpene(pgaus)
  real(rp),    intent(out)   :: gptmo(ndime,ndime,ndime,ndime,pgaus)
  real(rp),    intent(out)   :: ggsgo(ndime,ndime)
  real(rp),    intent(in)    :: gptlo(pgaus)
  real(rp),    intent(out)   :: gpdds(ndime,ndime,ndime,ndime,pgaus)
  integer(ip)                :: igaus,idime,jdime,kdime,ldime,mm,nn
  real(rp)                   :: tkron(ndime,ndime)
  !   real(rp)                   :: gpgre(ndime,ndime,pgaus,2),trace(pgaus,2)
  real(rp)                   :: dgree,dstre
  !   real(rp)                   :: gpcrt(ndime,ndime),eval(ndime),evec(ndime,ndime)
  real(rp),    intent(in)    :: elcod(ndime,mnode)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)

  real(rp),    intent(in)        :: gpgdi_eps(ndime,ndime,pgaus)           ! Displacement Gradient F for the perturbed state
  real(rp),    intent(in)        :: gpigd_eps(ndime,ndime,pgaus)           ! Displacement Gradient Inverse F^{-1} for the perturbed state
  real(rp),    intent(in)        :: gpdet_eps(pgaus)

  logical                    :: debuggingMode

  ! declaration of flags for additional inputs/outputs

  integer(ip)  :: flgs1,flgs3

  flgs1 = 0_ip ! flag for 1st P-K Stress tensor from 2nd P-K in output P(n+1)
  flgs3 = 0_ip ! flag for tangent moduli from material stress tangent  dP/dF

  ! touch the values of the "output" so that the compiler does think it is used
  ! (to be erased as soon as there exists one subroutine that uses it)

  ggsgo(1,1)       = 0.0_rp
  gptmo(1,1,1,1,1) = 0.0_rp

  debuggingMode = .false. ! (ielem == 1_ip)

  ! Kronecker delta

  tkron = 0.0_rp
  do idime = 1, ndime
     tkron(idime,idime) = 1.0_rp
  end do

  ! define additional inputs/outputs when needed

  if (lawst_sld(pmate)==100) then                     ! isolinear elastic
     flgs1 = 1_ip
     flgs3 = 1_ip
  else if (lawst_sld(pmate)==101) then                ! Neo-Hookean belytschko (textbook)
     flgs1 = 1_ip
     flgs3 = 1_ip
  else if (lawst_sld(pmate)==102) then                ! Neo-Hookean Ansys
     flgs1 = 1_ip
     flgs3 = 1_ip
  else if (lawst_sld(pmate)==103) then                ! Mooney-Rivlin
     flgs1 = 1_ip
     flgs3 = 1_ip
  else if (lawst_sld(pmate)==105) then                ! constant spring
     flgs1 = 1_ip
  else if (lawst_sld(pmate)==133) then                ! Yin & Lin
     flgs1 = 1_ip
  else if (lawst_sld(pmate)==1331) then                ! Yin & Lin
     flgs1 = 1_ip
  else if (lawst_sld(pmate)==134) then                ! Holza & Ogden
     flgs1 = 1_ip
     flgs3 = 1_ip
  else if (lawst_sld(pmate)==135) then                ! Pole-Zero
     flgs1 = 1_ip
  else if (lawst_sld(pmate)==136) then                ! Guccione
     flgs1 = 1_ip
     flgs3 = 1_ip
  else if (lawst_sld(pmate)==137) then                !HGO-C for arteries with 2 fiber fields
     flgs1 = 1_ip
     flgs3 = 1_ip
  else if (lawst_sld(pmate)==151) then                ! Ortholinear elastic
     flgs1 = 1_ip
     flgs3 = 1_ip
  else if (lawst_sld(pmate)==152) then                ! Transversally isotropic damage (Maimi 2017)
     flgs1 = 1_ip
     flgs3 = 1_ip
  else if (lawst_sld(pmate)==153) then                ! Isolinear damage law (Oliver 1990)
     flgs1 = 1_ip
     flgs3 = 1_ip
  else if (lawst_sld(pmate)==154) then                ! Transversally isotropic damage (XXX)
     flgs1 = 1_ip
     flgs3 = 0_ip
  else if (lawpl_sld(pmate)==1) then                  ! Plasticity model
     flgs1 = 1_ip
     flgs3 = 1_ip
  else if (lawst_sld(pmate)==MAT_MICRO_NO_COUPLING .or. &
           lawst_sld(pmate)==MAT_MICRO_ONE_WAY .or. &
           lawst_sld(pmate)==MAT_MICRO_FULL) then     ! MicroPP microscale model
     flgs1 = 1_ip
     flgs3 = 1_ip
  end if

  ! call built-in material constitutive laws

  ! Use tangent moduli (only for implicit computations)
  if ((kfl_tange_sld == SLD_TANGENT_NUMERICAL .or. kfl_tange_sld == 2_ip) .and. flagt==1_ip ) then
     call sld_dertan(pgaus,pmate,temp0,temp1,&
          gpgi0,gpgdi,gpigd,gpcau,gpdet,gprat,gppio,gpstr,&
          gpene,flagt,gptmo,flags,ggsgo,ielem,gptlo,nfibe,elcod,pnode,lnods,gpsha,gpdds,gpmof)

  else

     if (lawst_sld(pmate)==100) then
                                      ! isolinear elastic
        call sld_stress_model_100(pgaus,pmate,gpgdi,gpstr,gpcau,gpdet,gptlo,gpigd,ielem,elcod,flagt,gpdds,&
             gpigd_eps,gpgdi_eps,gpdet_eps,&
             gpmof)
     else if (lawst_sld(pmate)==101) then                                    ! Neo-Hookean belytschko (textbook)
        call sld_stress_model_101(pgaus,pmate,gpcau,gpgdi,gpene,gpstr,gpdet,flagt,gpdds,&
             gpigd_eps,gpgdi_eps,gpdet_eps,&
             gpmof)
     else if (lawst_sld(pmate)==102) then                                    ! Neo-Hookean Ansys
        call sld_stress_model_102(pgaus,pmate,gpgdi,gpstr,gpdet,flagt,gpdds,&
             gpigd_eps,gpgdi_eps,gpdet_eps,&
             gpmof)
     else if (lawst_sld(pmate)==103) then                                    ! Mooney-Rivlin (Belytschko notation)
        call sld_stress_model_103(pgaus,pmate,gpgdi,gpstr,gpdet,flagt,gpdds,&
             gpigd_eps,gpgdi_eps,gpdet_eps,&
             gpmof)
     else if (lawst_sld(pmate)==105) then                                    ! constant spring
        call sld_stress_model_105(pgaus,pmate,gpstr,flagt,gpdds,&
             gpigd_eps,gpgdi_eps,gpdet_eps,&
             gpmof)
     else if (lawst_sld(pmate)==MAT_MICRO_NO_COUPLING .or. &
              lawst_sld(pmate)==MAT_MICRO_ONE_WAY .or. &
              lawst_sld(pmate)==MAT_MICRO_FULL) then     ! MicroPP microscale model
        call fe2_get_stress_and_ctan(ielem, pgaus, gpstr, gpdds)

     else if (lawst_sld(pmate)==133) then                                    ! Yin & Lin with Peteron's EC coupling integrated
        if (flagt == 1_ip)   &
             call runend('SLD_BUILTIN_MATERIALS: STRESS TANGENT NOT PROGRAMMED YET FOR THIS MODEL')
        call sld_stress_model_133(pgaus,pmate,gpcau,gpstr,gpdet,gptlo,ielem,elcod,&
             gpigd_eps,gpgdi_eps,gpdet_eps,&
             gpmof)
     else if (lawst_sld(pmate)==1331) then                                    ! Yin & Lin only passive properties
        call sld_stress_model_133a(pgaus,pmate,gpgdi,gpigd,gpcau,gpstr,gpdet,gptlo,ielem,elcod,&
             gpigd_eps,gpgdi_eps,gpdet_eps,&
             gpmof)
     else if (lawst_sld(pmate)==134) then
        if (lawho_sld(pmate)==1) then
!          call sld_stress_model_134(&
!          call sld_stress_model_134exp(&
          call sld_stress_model_134implicit(&
               pgaus,pmate,gpvol,gpcau,gpgdi,gpigd,gpstr,gpdet,gptlo,nfibe,ielem,elcod,&
               pnode,lnods,gpsha,flagt,gpdds,&
               gpigd_eps,gpgdi_eps,gpdet_eps,&
               gpmof)
        else if (lawho_sld(pmate)==2) then                                     ! polynomial isochoric term
           call sld_stress_model_134pol(pgaus,pmate,gpvol,gpcau,gpgdi,gpigd,gpstr,gpdet,gptlo,nfibe,ielem,elcod,&
                pnode,lnods,gpsha,flagt,gpdds,&
                gpigd_eps,gpgdi_eps,gpdet_eps,&
                gpmof)
        else if (lawho_sld(pmate)==3) then                                      ! deviatoric formulation
           call sld_stress_model_134dev(pgaus,pmate,gpvol,gpcau,gpgdi,gpigd,gpstr,gpdet,gptlo,nfibe,ielem,elcod,&
                pnode,lnods,gpsha,flagt,gpdds,&
                gpigd_eps,gpgdi_eps,gpdet_eps,&
                gpmof)
        end if
     else if (lawst_sld(pmate)==135) then
        if (flagt == 1_ip)   &
             call runend('SLD_BUILTIN_MATERIALS: STRESS TANGENT NOT PROGRAMMED YET FOR THIS MODEL')
        call sld_stress_model_135(pgaus,pmate,gpcau,gpgdi,gpigd,gpstr,gpdet,pnode,elcod,lnods,ielem,gpsha,nfibe,&
             gpigd_eps,gpgdi_eps,gpdet_eps,&
             gpmof)
     else if (lawst_sld(pmate)==136) then
        call sld_stress_model_136(pgaus,pmate,gpcau,gpgdi,gpigd,gpstr,gpdet,gptlo,nfibe,ielem,elcod,&
             pnode,lnods,gpsha,flagt,gpdds,&
             gpigd_eps,gpgdi_eps,gpdet_eps,&
             gpmof)
     else if (lawst_sld(pmate)==137) then
        ! Holza & Ogden for arteries
        if (lawho_sld(pmate)==4) then
           call sld_stress_model_137hgoc(pgaus,pmate,gpvol,gpcau,gpgdi,gpigd,gpstr,gpdet,gptlo,nfibe,ielem,elcod,&
                pnode,lnods,gpsha,flagt,gpdds,&
                gpigd_eps,gpgdi_eps,gpdet_eps,&
                gpmof)
        else if (lawho_sld(pmate)==5) then
           call sld_stress_model_137ma(pgaus,pmate,gpvol,gpcau,gpgdi,gpigd,gpstr,gpdet,gptlo,nfibe,ielem,elcod,&
                pnode,lnods,gpsha,flagt,gpdds,&
                gpigd_eps,gpgdi_eps,gpdet_eps,&
                gpmof)
        end if
     else if (lawst_sld(pmate)==151) then
        ! sm151
        call sld_stress_model_151(pgaus,pmate,gpgdi,gpigd,gpdet,gpstr,ielem,flagt,gpdds)
     else if (lawst_sld(pmate)==152) then
        ! sm152
        call sld_stress_model_152(pgaus,pmate,gpgdi,gpigd,gpdet,gpstr,ielem,flagt,gpdds)
     else if (lawst_sld(pmate)==153) then
        ! sm153
        call sld_stress_model_153(pgaus,pmate,gpgdi,gpstr,ielem,flagt,gpdds,&
             gpigd_eps,gpgdi_eps,gpdet_eps)
     else if (lawst_sld(pmate)==154) then
        ! sm154
        call sld_stress_model_154(pgaus,pmate,gpgdi,gpigd,gpdet,gpstr,ielem,flagt,gpdds)
     else if (lawpl_sld(pmate)==1 .and. lawst_sld(pmate)==0) then
        ! PLASTICITY MODELS
        call sld_plastic_model_biso(pgaus,gpgdi,gpigd,gpdet,gpstr,ielem,flagt,gpdds)
     end if

  end if

  if ( flgs1 == 1_ip ) then
     !
     ! Calculate P (1st Piola-Kirchhoff stress tensor)
     !
     ! P_jK = F_jI * S_IK
     !
     do igaus=1,pgaus
        do kdime=1,ndime
           do jdime=1,ndime
              gppio(jdime,kdime,igaus)=0.0_rp
              do idime=1,ndime
                 gppio(jdime,kdime,igaus)= &
                      gppio(jdime,kdime,igaus) + &
                      gpgdi(jdime,idime,igaus)*gpstr(idime,kdime,igaus,1)
              enddo
           enddo
        enddo
     end do
  endif

  if (kfl_ninex_sld == 1) then

     !
     ! Compute E and E^eps: esto lo voy a usar para cuando el material no calcule el green lagrange
     !
     ! call sld_greela(pgaus,gpgdi,gpgdi_eps,gpgre,trace)
     !
     ! Compute gpdds (dS/dE) numerically
     !
     do igaus= 1,pgaus
        do mm= 1,ndime
           do nn= 1,ndime
              do jdime=1,ndime
                 do ldime=1,ndime
                    dstre = gpstr(jdime,ldime,igaus,2) - gpstr(jdime,ldime,igaus,1)
                    !                    dgree = abs(gpgre(mm,nn,igaus,2)   - gpgre(mm,nn,igaus,1))
                    dgree = epsex_sld

                    !                    gpdds(mm,jdime,nn,ldime,igaus) = 0.0_rp
                    !                    if (dgree > 0.0_rp) then
                    !                       gpdds(mm,jdime,nn,ldime,igaus) = dstre / dgree
                    !                    end if
                    gpdds(jdime,ldime,mm,nn,igaus) = 0.0_rp
                    if (dgree > 0.0_rp) then
                       gpdds(jdime,ldime,mm,nn,igaus) = dstre / dgree
                    end if
                 end do
              end do
           end do
        end do
     end do

     gptmo= 0.0_rp

     !
     ! Calculate P^eps
     !
     ! P^eps_jK = F^eps_jI * S^eps_IK
     !
     do igaus=1,pgaus
        do kdime=1,ndime
           do jdime=1,ndime
              gppio_eps(jdime,kdime,igaus)=0.0_rp
              do idime=1,ndime
                 gppio_eps(jdime,kdime,igaus)= &
                      gppio_eps(jdime,kdime,igaus) + &
                      gpgdi_eps(jdime,idime,igaus)*gpstr(idime,kdime,igaus,2)
              enddo
           enddo
        enddo
     end do
  end if

  gptmo= 0.0_rp
  if ( flgs3*flagt == 1_ip ) then
     !!  if ((kfl_ninex_sld == 0) .and. ( flgs3*flagt == 1_ip )) then
     !!  if (kfl_ninex_sld == 1) then
     !Calculate dPdF (tangent of 1st Piola-Kirchhoff stress tensor wrt. deformation gradient)
     ! dPdF_iJkL = delta_ik * S_JL + F_iM * FkN * dSdE_MJNL
     !
     do igaus=1,pgaus
        do idime=1,ndime; do jdime=1,ndime; do kdime=1,ndime; do ldime=1,ndime
           gptmo(idime,jdime,kdime,ldime,igaus) = 0.0_rp
           do mm=1,ndime; do nn=1,ndime
              gptmo(idime,jdime,kdime,ldime,igaus)= &
                   gptmo(idime,jdime,kdime,ldime,igaus) + &
                   gpdds(mm,jdime,nn,ldime,igaus)*gpgdi(idime,mm,igaus)*gpgdi(kdime,nn,igaus)
           enddo; enddo
           gptmo(idime,jdime,kdime,ldime,igaus)=&
                gptmo(idime,jdime,kdime,ldime,igaus) + &
                tkron(idime,kdime)*gpstr(jdime,ldime,igaus,1)
        enddo; enddo; enddo; enddo
     enddo

     if (debuggingMode) then
        write(*,*) ''
        write(*,*) 'gpdds = '
        do idime=1,ndime
           write(*,*)''
           do kdime=1,ndime
              write(*,'(1x,3(1x,3(1x,e10.3)))')((gpdds(idime,jdime,kdime,ldime,1),ldime=1,ndime),jdime=1,ndime)
           enddo
        enddo
        write(*,*) ''
        write(*,*) 'gptmo = '
        do idime=1,ndime
           write(*,*)''
           do kdime=1,ndime
              write(*,'(1x,3(1x,3(1x,e10.3)))')((gptmo(idime,jdime,kdime,ldime,1),ldime=1,ndime),jdime=1,ndime)
           enddo
        enddo
     endif

  endif

end subroutine sld_builtin_materials


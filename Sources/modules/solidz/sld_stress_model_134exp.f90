!$-----------------------------------------------------------------------
!> @addtogroup SolidzMaterials
!> @{
!> @file    sld_stress_model_134.f90
!> @author  Eva Casoni 01/07/2017 
!> @date    01/08/2017
!> @brief   Law of Holzapfel et Ogden modified for compressible anisotropic (Nolan, 2014)
!> @details Law of HGO-MA
!> @} 
!-----------------------------------------------------------------------
subroutine sld_stress_model_134exp(&
     pgaus,pmate,gpvol,gpcau,gpgdi,gpigd,gpstr,gpdet,gptlo,nfibe,ielem,elcod,&
     pnode,lnods,gpsha,flagt,gpdds,&
     gpigd_eps,gpgdi_eps,gpdet_eps,&
     gpmof)

  !-----------------------------------------------------------------------
  !****f* Solidz/sld_stress_model_134ma
  ! NAME
  !    sld__stress_model_134ma
  ! DESCRIPTION
  !    Law of Holzapfel et Ogden (2009)
  !
  !    GPGDI ... Deformation tensor ...................... F = grad(u) + I
  !    GPIGD ... Inverse of Deformation tensor ............F^(-1)
  !    GPCAU ... Right Cauchy-Green deformation tensor ... C = F^t x F
  !    GPCAL ... Left Cauchy-Green deformation tensor ...  b = F x F^T
  !    GPSTR ... 2nd P-K Stress tensor ....................S
  !    CASTR ... Cauchy Stress tensor .................... sigma
  !    GPPIO ... 1st Piola-Kirchhoff stress tensor ....... P = F.S
  !    GPENE ... Stored energy function .................. W
  !    GPLEP ... Log strain in the {f s n} system
  !    GPDDS ... Tangent moduli in terms of 2nd P-K stress ....... dS/dE
  !
  !    Special postproc values for this material:
  !    OUTPUT_&_POST-PROCESS
  !    POSTPROCESS LOCEPSILON,   => Log strain in the fibers CS - ln(lambda)
  !    END_OUTPUT_&_POST_PROCESS
  !
  ! USES
  ! USED BY
  !    sld_elmope
  !***
  !-------------------------------------------- ---------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode,lmate,npoin,coord,vmass
  use def_solidz
  use def_master, only       :  &
       cutim,ittim,postp,ittyp,ITASK_ENDRUN,fiber,gpfib,vconc,coupling,kfl_paral,fisoc, kfl_eccty, &
       cell_ca0_ecc
  use def_master, only       :  kfl_cellmod,CELL_TT2006_EXMEDI,CELL_OHARA_EXMEDI,CELL_FITZHUGH_EXMEDI,CELL_NOMOD_EXMEDI


  implicit none
  integer(ip), intent(in)    :: pgaus,pmate,pnode,lnods(pnode),flagt
  real(rp),    intent(in)    :: gpcau(ndime,ndime,pgaus),gpdet(pgaus),gpgdi(ndime,ndime,pgaus),gpmof(pgaus)
  real(rp),    intent(in)    :: gpigd(ndime,ndime,pgaus),gpsha(pnode,pgaus)
  real(rp),    intent(out)   :: gpstr(ndime,ndime,pgaus,2)
  real(rp),    intent(out)   :: gpdds(ndime,ndime,ndime,ndime,pgaus)
  real(rp),    intent(in)        :: gpgdi_eps(ndime,ndime,pgaus)           ! Displacement Gradient F for the perturbed state
  real(rp),    intent(in)        :: gpigd_eps(ndime,ndime,pgaus)           ! Displacement Gradient Inverse F^{-1} for the perturbed state
  real(rp),    intent(in)        :: gpdet_eps(pgaus)                       ! 
  real(rp)                   :: detf0,nfibe_length,tactf,nfibe(ndime,pgaus),nfibe_2(ndime,pgaus)

  integer(ip)                :: igaus,idime,jdime,kdime,ldime,ielem,kangl,inode,ipoin,ikact
  real(rp)                   :: gpcal(ndime,ndime),castr(ndime,ndime)
  real(rp)                   :: castr_active(ndime,ndime),gpstr_active(ndime,ndime)

  real(rp)                   :: nfibe0(ndime,pgaus),nfibe0_2(ndime,pgaus), norma0(ndime,pgaus), nshet0(ndime,pgaus) !f0, s0, n0 (ref. config)
  real(rp)                   ::                     norma(ndime,pgaus) ,nshet(ndime,pgaus)  !f,s,n (fF as define in paper) (nfibe is "output")
  real(rp)                   :: nfibt(ndime,pgaus) ,normt(ndime,pgaus) ,nshtt(ndime,pgaus)  !f,s,n (updated-Unit)

  real(rp)                   :: a,b,af,bf,as,bs,afs,bfs,scalf

  real(rp)                   :: i1,i4f,i4s,i4n,i8fs
  real(rp)                   :: term1
  real(rp)                   :: gpcai(ndime,ndime),Kct

  real(rp)                   :: lamda(3),bidon,dummr,ovlam(3)

  real(rp)                   :: elfib(3,mnode),elfis(3,mnode),elfin(3,mnode),elcac(mnode),gpcac,ca0

  real(rp),    intent(in)    :: gptlo(pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus)
  real(rp),    intent(in)    :: elcod(ndime,mnode)

  real(rp)                   :: tkron(3,3), statvar(6,2)

  real(rp)                   :: fporf(ndime,ndime), fporf_0(ndime,ndime), &
       spors(ndime,ndime), spors_0(ndime,ndime), sporf(ndime,ndime), fpors(ndime,ndime)

  real(rp)                   :: dev_fporf(ndime,ndime), dev_spors(ndime,ndime),dev_fpors(ndime,ndime)

  real(rp) :: fpors_0(ndime,ndime), spk_vol(ndime,ndime), spk_iso(ndime,ndime), spk_4f(ndime,ndime), spk_4s(ndime,ndime), &
       spk_8fs(ndime,ndime), ixinvc(ndime,ndime,ndime,ndime), invcxi(ndime,ndime,ndime,ndime), invcxinvc(ndime,ndime,ndime,ndime), &
       i_invc(ndime,ndime,ndime,ndime), i1xi1(ndime,ndime,ndime,ndime), fxfxinvc(ndime,ndime,ndime,ndime), invcxfxf(ndime,ndime,ndime,ndime), &
       i4fxi4f(ndime,ndime,ndime,ndime), sxsxinvc(ndime,ndime,ndime,ndime), invcxsxs(ndime,ndime,ndime,ndime), i4sxi4s(ndime,ndime,ndime,ndime), &
       fsxfsxinvc(ndime,ndime,ndime,ndime), invcxfsxfs(ndime,ndime,ndime,ndime), i8fsxi8fs(ndime,ndime,ndime,ndime), dev_ipori(ndime,ndime), &
       h4f, h4s, J, J23, J43, tan_opt_aniso_4f(ndime,ndime,ndime,ndime), tan_opt_aniso_4s(ndime,ndime,ndime,ndime), &
       tan_opt_aniso_8fs(ndime,ndime,ndime,ndime), &
       tan_opt_iso(ndime,ndime,ndime,ndime), tan_opt_vol(ndime,ndime,ndime,ndime), term4f, term4s, term8fs, &
       term_vol, iso_i1, iso_i4f, iso_i4n, iso_i4s, iso_i8fs, &
       fxfxfxf(ndime,ndime,ndime,ndime), beta, Cal50, cca2p, n, tan_opt_act(ndime,ndime,ndime,ndime), tmaxi, term_act(ndime,ndime,ndime,ndime)


  real(rp) :: temp(3,3,3,3)



  ! Retrieve the parameters
  scalf = 1.0_rp !default value

  Kct = parco_sld(2,pmate) !0.001 N/cm2 , values of Usyk, 2000
  a   = parco_sld(3,pmate)
  b   = parco_sld(4,pmate)
  af  = parco_sld(5,pmate)
  bf  = parco_sld(6,pmate)
  as  = parco_sld(7,pmate)
  bs  = parco_sld(8,pmate)
  afs = parco_sld(9,pmate)
  bfs = parco_sld(10,pmate)
  scalf = parco_sld(11,pmate)

  ! Initialise some needed matrices and vectors
  gpcai = 0.0_rp
  gpstr = 0.0_rp
  gpdds = 0.0_rp

  nfibt = 0.0_rp
  nshtt = 0.0_rp
  normt = 0.0_rp

  nfibe = 0.0_rp
  nshet = 0.0_rp
  norma = 0.0_rp

  nfibe0 = 0.0_rp
  norma0 = 0.0_rp
  nshet0 = 0.0_rp

  elfib = 0.0_rp

  tkron = 0.0_rp
  do idime = 1,ndime
     tkron(idime,idime) = 1.0_rp
  end do

  !
  ! Gather
  !
  if( modfi_sld(pmate) .ne. 0 ) then
     elfis = 0.0_rp           
     elfin = 0.0_rp
     if (kfl_fiber_sld > 0) then
        if (modfi_sld(pmate) < 0) then
           do inode = 1,pnode
              ipoin = lnods(inode)
              elfib(1:ndime,inode) = fiber(1:ndime,ipoin)
           end do
        else if (modfi_sld(pmate) > 0) then
           gfibe_sld = 0.0_rp
           gfibe_sld(modfi_sld(pmate)) = 1.0_rp
           do inode = 1,pnode
              elfib(1:ndime,inode) = gfibe_sld(1:ndime)
           end do
        end if
        if (kfl_fiber_sld == 2 .or. kfl_fiber_sld == 3) then
           do inode = 1,pnode
              ipoin = lnods(inode)
              elfis(1,inode) = fibts_sld(1,ipoin)
              elfis(2,inode) = fibts_sld(2,ipoin)
              elfis(3,inode) = fibts_sld(3,ipoin)    

              elfin(1,inode) = fibtn_sld(1,ipoin)
              elfin(2,inode) = fibtn_sld(2,ipoin)
              elfin(3,inode) = fibtn_sld(3,ipoin)  
           end do
        end if
     end if
  end if

  if(( coupling('SOLIDZ','EXMEDI') >= 1 ) .or. ( coupling('EXMEDI','SOLIDZ') >= 1 )) then
     !
     ! Coupling stragegy
     !

     if ((kfl_cellmod(pmate) == CELL_TT2006_EXMEDI).or. &
          (kfl_cellmod(pmate) == CELL_OHARA_EXMEDI)) then
        ikact = 1
        ! coupling depends on Ca2+
        do inode = 1,pnode
           ipoin = lnods(inode)
           elcac(inode) = vconc(1,ipoin,1) *1000.0_rp  !unit conversion between EP and ECC model
        end do
     else if (kfl_cellmod(pmate) == CELL_FITZHUGH_EXMEDI)  then  !! FHN
        ikact = 0
        ! coupling depends on the activation index kacti
        do inode = 1,pnode
           ipoin = lnods(inode)
           ikact = ikact + kacti_sld(ipoin)
        end do
     else if (kfl_cellmod(pmate) == CELL_NOMOD_EXMEDI)  then  !
        ikact = 0
     else
        call runend ("SLD_MYOFIL: CELL MODEL NOT INCLUDED IN THE COUPLING WITH SOLIDZ MODULE")
     end if
  end if

  ! Gauss points loop
  do igaus = 1,pgaus

     !Initialize active, passive and volumetric stresses to zero
     gpcai = 0.0_rp
     castr_active = 0.0_rp

     !Calcul of b = F F^T
     gpcal = 0.0_rp
     do idime = 1,ndime
        do jdime = 1,ndime
           do kdime = 1,ndime
              gpcal(idime,jdime) = gpcal(idime,jdime)+gpgdi(idime,kdime,igaus)*gpgdi(jdime,kdime,igaus)
           end do
        end do
     end do

     if( modfi_sld(pmate) .ne. 0 ) then
        !
        ! Fibers interpolated on Gauss point
        !
        nfibe0(1:3,igaus) = 0.0_rp
        do inode = 1,pnode
           nfibe0(1:3,igaus) = nfibe0(1:3,igaus) + gpsha(inode,igaus)*elfib(1:3,inode)
        end do

        !Normalized f_0
        bidon = sqrt(nfibe0(1,igaus)**2.0_rp+nfibe0(2,igaus)**2.0_rp+nfibe0(3,igaus)**2.0_rp)
        if (bidon == 0.0_rp) bidon= 1.0_rp
        do idime = 1,ndime
           nfibe0(idime,igaus) = nfibe0(idime,igaus)/bidon
        end do
     end if

     if( kfl_fiber_sld == 2 .or. kfl_fiber_sld == 3) then

        !
        ! Orthotropic material
        !
        norma0(1:3,igaus) = 0.0_rp
        nshet0(1:3,igaus) = 0.0_rp
        do inode = 1,pnode
           norma0(1:3,igaus) = norma0(1:3,igaus)+gpsha(inode,igaus)*elfin(1:3,inode)
           nshet0(1:3,igaus) = nshet0(1:3,igaus)+gpsha(inode,igaus)*elfis(1:3,inode)
        end do

        !Normalized n_0
        bidon = sqrt(norma0(1,igaus)**2.0_rp+norma0(2,igaus)**2.0_rp+norma0(3,igaus)**2.0_rp)
        if (bidon == 0.0_rp) bidon= 1.0_rp
        do idime=1,ndime
           norma0(idime,igaus)=norma0(idime,igaus)/bidon
        end do

        !Normalized s_0
        bidon = sqrt(nshet0(1,igaus)**2.0_rp+nshet0(2,igaus)**2.0_rp+nshet0(3,igaus)**2.0_rp)
        if (bidon == 0.0_rp) bidon = 1.0_rp
        do idime = 1,ndime
           nshet0(idime,igaus) = nshet0(idime,igaus)/bidon
        end do
     end if

     ! * * * *
     !INVARIANTS
     ! * * * *

     !I_1=Trace(C)
     i1 = 0.0_rp
     do idime = 1,ndime
        i1 = i1 + gpcau(idime,idime,igaus)
     end do

     ! Compute the invariants
     i4f = 0.0_rp
     i4s = 0.0_rp
     i4n = 0.0_rp
     i8fs = 0.0_rp

     do idime = 1,ndime
        do jdime = 1,ndime
           i4f  =  i4f+nfibe0(idime,igaus)*gpcau(idime,jdime,igaus)*nfibe0(jdime,igaus) ! I_4f  = f0 * C_ij * f0 (eq. 5.1)
           i4s  =  i4s+nshet0(idime,igaus)*gpcau(idime,jdime,igaus)*nshet0(jdime,igaus) ! I_4s  = s0 * C_ij * s0
           i4n  =  i4n+norma0(idime,igaus)*gpcau(idime,jdime,igaus)*norma0(jdime,igaus) ! I_4n  = n0 * C_ij * n0
           i8fs = i8fs+0.5_rp*(nfibe0(idime,igaus)*gpcau(idime,jdime,igaus)*nshet0(jdime,igaus)+ &
                nshet0(idime,igaus)*gpcau(idime,jdime,igaus)*nfibe0(jdime,igaus)) ! I_8fs = f0 * C_ij * s0 (eq. 5.3)

           !Transform f=F'*f0 ; s=F'*s0 ; n=F'*n0 (PUSH FORWARD)
           nfibe(idime,igaus) = nfibe(idime,igaus) + gpgdi(idime,jdime,igaus) * nfibe0(jdime,igaus)
           nshet(idime,igaus) = nshet(idime,igaus) + gpgdi(idime,jdime,igaus) * nshet0(jdime,igaus)
           norma(idime,igaus) = norma(idime,igaus) + gpgdi(idime,jdime,igaus) * norma0(jdime,igaus)
        end do
     end do

     do idime = 1,ndime
        gpfib(idime,igaus,ielem) = nfibe(idime,igaus)
     end do

     !
     ! f x f , s x s, f x s, s x f (outer products [3x3])
     !
     do idime = 1,ndime
        do jdime = 1,ndime
           fporf_0(idime,jdime) = nfibe0(idime,igaus)*nfibe0(jdime,igaus)
           spors_0(idime,jdime) = nshet0(idime,igaus)*nshet0(jdime,igaus)
           fpors_0(idime,jdime) = 0.5_rp*(nfibe0(idime,igaus)*nshet0(jdime,igaus)+nshet0(idime,igaus)*nfibe0(jdime,igaus))
        end do
     end do

     ! * * * * * *
     ! Second Piola-Kirchhoff stress tensor
     ! * * * * * *

     ! Heaviside term for tension/compression asymmetry
     h4f = 0.0_rp
     h4s = 0.0_rp
     if (i4f .gt. 1.0_rp) h4f = 1.0_rp
     if (i4s .gt. 1.0_rp) h4s = 1.0_rp

     ! volumetric terms
     J   = gpdet(igaus)
     J23 = gpdet(igaus)**(-(2.0_rp/3.0_rp))
     J43 = gpdet(igaus)**(-(4.0_rp/3.0_rp))
     term_vol = Kct*(J-1.0_rp)*J

     ! Compute the isochoric invariants
     iso_i1 = J23*i1
     iso_i4f = J23*i4f
     iso_i4s = J23*i4s
     iso_i4n = J23*i4n
     iso_i8fs = J23*i8fs

     ! exponential isotropic and isochoric term
     term1 = exp(b*(iso_i1-3.0_rp))

     if (kfl_fiber_sld > 0) then

        ! these checks re necessary, because round-off errors from i4f/s can create spurious small forces in the passive 
        if (abs((i4f-1.0_rp)) < 1.0e-12_rp) i4f = 1.0_rp
        if (abs((i4s-1.0_rp)) < 1.0e-12_rp) i4s = 1.0_rp

        ! exponential anisotropic and isochoric term
        term4f = exp(bf*(iso_i4f-1.0_rp)**2.0_rp)
        term4s = exp(bs*(iso_i4s-1.0_rp)**2.0_rp)
        term8fs = exp(bfs*(iso_i8fs**2.0_rp))

     end if

     ! Invert Right Cauchy-Green strain tensor C (gpcai)
     call invmtx(gpcau(:,:,igaus),gpcai,bidon,ndime)

     !
     ! Deviatoric part: Dev(Aij) = Aij - (1/3)Akk tkron_ij
     !
     dev_ipori = tkron-(1.0_rp/3.0_rp)*i1*gpcai
     dev_fporf = fporf_0-(1.0_rp/3.0_rp)*i4f*gpcai
     dev_spors = spors_0-(1.0_rp/3.0_rp)*i4s*gpcai
     dev_fpors = fpors_0-(1.0_rp/3.0_rp)*i8fs*gpcai

     ! Second Piola-Kirchhoff stress tensor
     spk_vol = 2.0_rp*term_vol*gpcai
     spk_iso = J23*a*term1*dev_ipori
     spk_4f = h4f*(2.0_rp*J23*af*(iso_i4f-1.0_rp)*term4f*dev_fporf)
     spk_4s = h4s*(2.0_rp*J23*as*(iso_i4s-1.0_rp)*term4s*dev_spors)
     spk_8fs = 2.0_rp*J23*afs*iso_i8fs*term8fs*dev_fpors

     gpstr(1:3,1:3,igaus,1) = spk_vol(1:3,1:3)+spk_iso(1:3,1:3)+spk_4f(1:3,1:3)+ & 
          spk_4s(1:3,1:3)+spk_8fs(1:3,1:3)

     ! Calculate the stretch (lambda). In the fiber direction: lambda^2=f_0*C*f_0 = I4f
     lamda(1) = sqrt(i4f)
     lamda(2) = sqrt(i4s)
     lamda(3) = sqrt(i4n)

     !
     ! Updated fiber direction: lambda_f * f_true = F * f_0 = f
     ! => here f is the UNIT fiber direction in the current configuration
     ! used for active stress and postproc
     !
     nfibe_length= 0.0_rp
     do idime = 1,ndime
        nfibt(idime,igaus) = nfibe(idime,igaus)
        nshtt(idime,igaus) = nshet(idime,igaus)
        normt(idime,igaus) = norma(idime,igaus)
        nfibe_length = nfibe_length + nfibe(idime,igaus)*nfibe(idime,igaus) 
     end do

     ! Electromechanical coupling
     ! * * * * * *
     ! Active stress T(lambda, [Ca++])m x m  (Hunter p 688)
     ! * * * * * *

     if(( coupling('SOLIDZ','EXMEDI') >= 1 ) .or. ( coupling('EXMEDI','SOLIDZ') >= 1 )) then

        gpcac= 0.0_rp
        do inode = 1,pnode
           gpcac= gpcac + gpsha(inode,igaus) * elcac(inode)   
        end do

        if((kfl_eccty(pmate) == 1 .or. kfl_eccty(pmate) == 3 .or. kfl_eccty(pmate) == 4) .and. (ikact>0)) then
           ! hunter and land models

           !if  (kfl_eccty(pmate) == 3) then
           !   statvar= 0.0_rp
           !   do inode = 1,pnode
           !      ipoin = lnods(inode)
           !      statvar(1:6,1) =  statvar(1:6,1) + stateland(1:6,ipoin,1) * gpsha(inode,igaus)
           !   end do
           !end if

           !call sld_eccoup(&
           !     pgaus,pmate,igaus,gpdet,gptlo,lamda,nfibt,nshtt,normt,statvar,&
           !     gpigd,castr_active,gpstr_active,gpcac,tactf,ielem)             


           !if  (kfl_eccty(pmate) == 3) then
           !   do inode = 1,pnode
           !      ipoin = lnods(inode)
           !      stateland(1:6,ipoin,2) = &
           !           stateland(1:6,ipoin,2) + statvar(1:6,1) * gpvol(igaus) * gpsha(inode,igaus) / vmass(ipoin)
           !   end do
           !end if


        else if((kfl_eccty(pmate) == 2) .and. (ikact>0)) then
           ! rice model

           call sld_myofil(pgaus,pmate,igaus,gpdet,gptlo,lamda,nfibt,nshtt,normt,nfibe_length,&
                gpigd,castr_active,gpstr_active,gpcac,tactf,ielem)
        end if

        ! Convert the active Cauchy stress to Second Piola-Kirchhoff stress and add it to the passive component
        !castr_active = castr_active/gpdet(igaus)
        gpstr_active = 0.0_rp
        do idime = 1,ndime
           do jdime = 1,ndime
              do kdime = 1,ndime
                 do ldime = 1,ndime
                    gpstr_active(idime,jdime) = gpstr_active(idime,jdime)&
                         + gpdet(igaus)*gpigd(idime,kdime,igaus)*castr_active(kdime,ldime)*gpigd(jdime,ldime,igaus)
                 end do
              end do
           end do
        end do

        gpstr(1:3,1:3,igaus,1) = gpstr(1:3,1:3,igaus,1)+gpstr_active(1:3,1:3)

     end if

     ! * * * * * * * *
     ! Tangent moduli
     ! * * * * * * * *

     if (flagt == 1_ip) then

        !ixinvc = 0.0_rp
        !invcxi = 0.0_rp
        !invcxinvc = 0.0_rp
        !i_invc = 0.0_rp
        !i1xi1 = 0.0_rp
        !fxfxinvc = 0.0_rp
        !invcxfxf = 0.0_rp
        !i4fxi4f = 0.0_rp
        !sxsxinvc = 0.0_rp
        !invcxsxs = 0.0_rp
        !i4sxi4s = 0.0_rp
        !fsxfsxinvc = 0.0_rp
        !invcxfsxfs = 0.0_rp
        !i8fsxi8fs = 0.0_rp
        !fxfxfxf = 0.0_rp

        !ECR: de momento falta el termino de i8fs (el termino cruzado)

        do idime = 1,ndime
           do jdime = 1,ndime
              do kdime = 1,ndime
                 do ldime = 1,ndime
                    ixinvc(idime,jdime,kdime,ldime) = tkron(idime,jdime)*gpcai(kdime,ldime)
                    invcxi(idime,jdime,kdime,ldime) = gpcai(idime,jdime)*tkron(kdime,ldime)
                    invcxinvc(idime,jdime,kdime,ldime) = gpcai(idime,jdime)*gpcai(kdime,ldime)
                    i_invc(idime,jdime,kdime,ldime) = (gpcai(idime,kdime)*gpcai(jdime,ldime)+gpcai(idime,ldime)*gpcai(jdime,kdime))
                    i1xi1(idime,jdime,kdime,ldime) = dev_ipori(idime,jdime)*dev_ipori(kdime,ldime)
                    fxfxinvc(idime,jdime,kdime,ldime) = fporf_0(idime,jdime)*gpcai(kdime,ldime)
                    invcxfxf(idime,jdime,kdime,ldime) = gpcai(idime,jdime)*fporf_0(kdime,ldime)
                    i4fxi4f(idime,jdime,kdime,ldime) = dev_fporf(idime,jdime)*dev_fporf(kdime,ldime)
                    sxsxinvc(idime,jdime,kdime,ldime) = spors_0(idime,jdime)*gpcai(kdime,ldime)
                    invcxsxs(idime,jdime,kdime,ldime) = gpcai(idime,jdime)*spors_0(kdime,ldime)
                    i4sxi4s(idime,jdime,kdime,ldime) = dev_spors(idime,jdime)*dev_spors(kdime,ldime)
                    fsxfsxinvc(idime,jdime,kdime,ldime) = fpors_0(idime,jdime)*gpcai(kdime,ldime)
                    invcxfsxfs(idime,jdime,kdime,ldime) = gpcai(idime,jdime)*fpors_0(kdime,ldime)
                    i8fsxi8fs(idime,jdime,kdime,ldime) = dev_fpors(idime,jdime)*dev_fpors(kdime,ldime)
                    fxfxfxf(idime,jdime,kdime,ldime) = nfibe0(idime,igaus)*nfibe0(jdime,igaus)*nfibe0(kdime,igaus)*nfibe0(ldime,igaus)
                    term_act(idime,jdime,kdime,ldime) = gpstr_active(idime,jdime)*gpcai(kdime,ldime)
                 end do
              end do
           end do
        end do

        tan_opt_iso = &
             2.0_rp*J23*a*term1*(-(1.0_rp/3.0_rp)*ixinvc-(1.0_rp/3.0_rp)*invcxi &
             +(1.0_rp/9.0_rp)*i1*invcxinvc+(1.0_rp/6.0_rp)*i1*i_invc)+2.0_rp*J43*a*b*term1*(i1xi1)

        tan_opt_aniso_4f = &
             h4f*(4.0_rp*J23*af*(iso_i4f-1.0_rp)*term4f*(-(1.0_rp/3.0_rp)*fxfxinvc-(1.0_rp/3.0_rp)*invcxfxf &
             +(1.0_rp/9.0_rp)*i4f*invcxinvc+(1.0_rp/6.0_rp)*i4f*i_invc)+ &
             4.0_rp*J43*af*term4f*i4fxi4f+8.0_rp*J43*af*bf*((iso_i4f-1.0_rp)**2.0_rp)*term4f*i4fxi4f)

        tan_opt_aniso_4s = &
             h4s*(4.0_rp*J23*as*(iso_i4s-1.0_rp)*term4s*(-(1.0_rp/3.0_rp)*sxsxinvc-(1.0_rp/3.0_rp)*invcxsxs &
             +(1.0_rp/9.0_rp)*i4s*invcxinvc+(1.0_rp/6.0_rp)*i4s*i_invc)+ &
             4.0_rp*J43*as*term4s*i4sxi4s+8.0_rp*J43*as*bs*((iso_i4s-1.0_rp)**2.0_rp)*term4s*i4sxi4s)

        tan_opt_aniso_8fs = &
             4.0_rp*J23*afs*iso_i8fs*term8fs*(-(1.0_rp/3.0_rp)*fsxfsxinvc-(1.0_rp/3.0_rp)*invcxfsxfs &
             +(1.0_rp/9.0_rp)*i8fs*invcxinvc+(1.0_rp/6.0_rp)*i8fs*i_invc)+ &
             4.0_rp*J43*afs*term8fs*i8fsxi8fs+8.0_rp*J43*afs*bfs*(iso_i8fs**2.0_rp)*term8fs*i8fsxi8fs

        tan_opt_vol = 4.0_rp*Kct*J*((2.0_rp*J-1.0_rp)*invcxinvc-0.5_rp*(J-1.0_rp)*i_invc)

        gpdds(1:3,1:3,1:3,1:3,igaus) = tan_opt_iso(1:3,1:3,1:3,1:3)+tan_opt_aniso_4f(1:3,1:3,1:3,1:3)+tan_opt_aniso_4s(1:3,1:3,1:3,1:3)+ &
             tan_opt_aniso_8fs(1:3,1:3,1:3,1:3)+tan_opt_vol(1:3,1:3,1:3,1:3)

        ! Active component of the tangent operator
        if(( coupling('SOLIDZ','EXMEDI') >= 1 ) .or. ( coupling('EXMEDI','SOLIDZ') >= 1 )) then
           if (kfl_eccty(pmate) == 1) then
              !Ca0 is not necesarily zero for the EP model.
              ca0=sum(cell_ca0_ecc(:,pmate))*1000.0_rp/3.0_rp
              beta = 1.45_rp 
              n = hillc_sld(pmate)
              cca2p = gpcac-ca0
              Cal50 = cal50_sld(pmate)
              tmaxi= 1.0e6_rp * cocof_sld(pmate)
              tan_opt_act = gpdet(igaus)*(term_act+(cca2p**n/(cca2p**n + Cal50**n))*tmaxi*(beta/lamda(1))*fxfxfxf)
              gpdds(1:3,1:3,1:3,1:3,igaus) = gpdds(1:3,1:3,1:3,1:3,igaus)+tan_opt_act(1:3,1:3,1:3,1:3)
           end if
        end if
     end if

  end do !gauss points

end subroutine sld_stress_model_134exp

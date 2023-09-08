!$-----------------------------------------------------------------------
!> @addtogroup SolidzMaterials
!> @ingroup    Solidz
!> @{
!> @file    sld_stress_model_134.f90
!> @author  Pierre Lafortune
!> @date    16/11/1966
!> @brief   Law of Holzapfel et Ogden (2009)
!> @details Law of Holzapfel et Ogden (2009)
!> @} 
!-----------------------------------------------------------------------
subroutine sld_stress_model_134(&
     pgaus,pmate,gpvol,gpcau,gpgdi,gpigd,gpstr,gpdet,gptlo,nfibe,ielem,elcod,&
     pnode,lnods,gpsha,flagt,gpdds,&
     gpigd_eps,gpgdi_eps,gpdet_eps,&
     gpmof)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_stress_model_134
  ! NAME
  !    sld__stress_model_134
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
       cutim,ittim,postp,ittyp,ITASK_ENDRUN,fiber,gpfib,vconc,coupling,kfl_paral,kfl_eccty,stretlam
  use def_master, only       :  kfl_cellmod,CELL_TT2006_EXMEDI,CELL_OHARA_EXMEDI,CELL_FITZHUGH_EXMEDI,CELL_NOMOD_EXMEDI


  implicit none
  integer(ip), intent(in)    :: pgaus,pmate,pnode,lnods(pnode),flagt
  real(rp),    intent(in)    :: gpcau(ndime,ndime,pgaus),gpdet(pgaus),gpgdi(ndime,ndime,pgaus),gpmof(pgaus)
  real(rp),    intent(in)    :: gpigd(ndime,ndime,pgaus),gpsha(pnode,pgaus)
  real(rp),    intent(out)   :: gpstr(ndime,ndime,pgaus,2),nfibe(ndime,pgaus)
  real(rp),    intent(out)   :: gpdds(ndime,ndime,ndime,ndime,pgaus)
  real(rp),    intent(in)        :: gpgdi_eps(ndime,ndime,pgaus)           ! Displacement Gradient F for the perturbed state
  real(rp),    intent(in)        :: gpigd_eps(ndime,ndime,pgaus)           ! Displacement Gradient Inverse F^{-1} for the perturbed state
  real(rp),    intent(in)        :: gpdet_eps(pgaus)                       ! 
  real(rp)                   :: detf0,nfibe_length,tactf

  integer(ip)                :: igaus,idime,jdime,kdime,ldime,ielem,kangl,inode,ipoin,ikact
  real(rp)                   :: gpcal(ndime,ndime),castr(ndime,ndime),castv(ndime,ndime),sleci(3,3,3)

  real(rp)                   :: nfibe0(ndime,pgaus),norma0(ndime,pgaus),nshet0(ndime,pgaus) !f0, s0, n0 (ref. config)
  real(rp)                   ::                     norma(ndime,pgaus) ,nshet(ndime,pgaus)  !f,s,n (fF as define in paper) (nfibe is "output")
  real(rp)                   :: nfibt(ndime,pgaus) ,normt(ndime,pgaus) ,nshtt(ndime,pgaus)  !f,s,n (updated-Unit)

  real(rp)                   :: a,b,af,afm,bf,as,bs,afs,bfs,afc

  real(rp)                   :: i1,i4f,i4s,i4n,i8fs,i4ff,i4sf
  real(rp)                   :: fporf(ndime,ndime),fporf_t(ndime,ndime),spors(ndime,ndime)
  real(rp)                   :: fpors(ndime,ndime),sporf(ndime,ndime)
  real(rp)                   :: term1,term2,term3,term4,termK
  real(rp)                   :: gpcai(ndime,ndime),bidon,Kct,dwvdc(ndime,ndime),trace

  real(rp)                   :: lamda(3),Cal50,n,norme1,norme2
  real(rp)                   :: scasta,scastr,scastv,sdwvdc,bidon2,dummr,ovlam(3)

  real(rp)                   :: castr_active(ndime,ndime),gpstr_active(ndime,ndime),gpgre(ndime,ndime)

  real(rp)                   :: xcord,ycord,zcord,vuni1(3),vuni2(3),&
       elfib(3,mnode),elfis(3,mnode),elfin(3,mnode),elcac(mnode),gpcac

  real(rp)                   :: lenght_pnt,angle1,omega,a_epi,b_epi,yp_epi,zp_epi,lenght_end,l_star
  real(rp)                   :: rotma(3,3),lenght_epi,yp_end,zp_end,a_end,b_end,nore1,scalf

  real(rp),    intent(in)    :: gptlo(pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus)
  real(rp),    intent(in)    :: elcod(ndime,mnode)

  real(rp)                   :: tkron(3,3)
  real(rp)                   :: statvar(6,2)

  !real(rp)                   :: Fpert(ndime,ndime,pgaus), Cpert(ndime,ndime), bpert(ndime,ndime), epsil





  scalf = 1.0_rp !default value 

  Kct = parco_sld( 2,pmate) !0.001 N/cm2 , values of Usyk, 2000
  a   = parco_sld( 3,pmate)
  b   = parco_sld( 4,pmate)
  af  = parco_sld( 5,pmate)
  bf  = parco_sld( 6,pmate)
  as  = parco_sld( 7,pmate)
  bs  = parco_sld( 8,pmate)
  afs = parco_sld( 9,pmate)
  bfs = parco_sld(10,pmate)
  scalf = parco_sld(11,pmate)


  !  if (scalf == 0.0_rp) then
  !       write(*,*) '--| ALYA  sld_stress_model_134. Must define a scaling parameter (11th position) '
  !       stop
  !  end if 


  !Now in sld_velsnd.f90
  ! Compute sound velocity only at time-step 0 (initialization)
  !if (ittim == 0_ip) then
  !   velas_sld(1,1) = sqrt(parco_sld(1,1)/densi_sld(1,1))
  !end if

  ! if (ielem==1) write(*,*) 'param1 ',Kct,a,b,af,bf,pmate
  ! if (ielem==84) write(*,*) 'param84 ',Kct,a,b,af,bf,pmate

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
  do idime=1,ndime
     tkron(idime,idime) = 1.0_rp
  end do

  !Fpert = 0.0_rp
  !epsil = 0.0000000001_rp


  ! levi-civita symbol
  sleci= 0.0_rp
  sleci(1,2,3)= 1.0_rp
  sleci(2,3,1)= 1.0_rp
  sleci(3,1,2)= 1.0_rp
  sleci(3,2,1)= 1.0_rp
  sleci(1,3,2)= 1.0_rp
  sleci(2,1,3)= 1.0_rp

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
           gfibe_sld= 0.0_rp
           gfibe_sld(modfi_sld(pmate))= 1.0_rp
           do inode = 1,pnode
              elfib(1:ndime,inode)= gfibe_sld(1:ndime)
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
  elcac = 0.0_rp

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
        ikact= 0
        ! coupling depends on the activation index kacti
        do inode = 1,pnode
           ipoin = lnods(inode)
           ikact= ikact + kacti_sld(ipoin)
        end do
     else if (kfl_cellmod(pmate) == CELL_NOMOD_EXMEDI)  then  !
        ikact= 0
     else 
        call runend ("SLD_MYOFIL: CELL MODEL NOT INCLUDED IN THE COUPLING WITH SOLIDZ MODULE")
     end if
  end if

  !write(*,*) "stress", kacti_sld(1:5), elcac

  do igaus = 1,pgaus

     !Initialize active, passive and volumetric stresses to zero
     castr = 0.0_rp
     castv = 0.0_rp
     castr_active = 0.0_rp
     gpstr_active = 0.0_rp
     gpcai = 0.0_rp
     statvar = 0.0_rp

     !Calcul of b = F F^T
     gpcal = 0.0_rp
     do idime = 1,ndime
        do jdime = 1,ndime
           do kdime = 1,ndime
              gpcal(idime,jdime) = &
                   gpcal(idime,jdime) + gpgdi(idime,kdime,igaus) * gpgdi(jdime,kdime,igaus)
           end do
        end do
     end do

     if( modfi_sld(pmate) .ne. 0 ) then
        !
        ! Fibers interpolated on Gauss point
        !
        nfibe0(1,igaus) = 0.0_rp
        nfibe0(2,igaus) = 0.0_rp
        nfibe0(3,igaus) = 0.0_rp
        do inode = 1,pnode
           nfibe0(1,igaus) = nfibe0(1,igaus) + gpsha(inode,igaus) * elfib(1,inode)
           nfibe0(2,igaus) = nfibe0(2,igaus) + gpsha(inode,igaus) * elfib(2,inode)
           nfibe0(3,igaus) = nfibe0(3,igaus) + gpsha(inode,igaus) * elfib(3,inode)
        end do

        !Normalized f_0
        bidon=sqrt(nfibe0(1,igaus)**2.0_rp+nfibe0(2,igaus)**2.0_rp+nfibe0(3,igaus)**2.0_rp)
        if (bidon == 0.0_rp) bidon= 1.0_rp
        do idime=1,ndime
           nfibe0(idime,igaus)=nfibe0(idime,igaus)/bidon
        end do

     end if

     if( kfl_fiber_sld == 2 .or. kfl_fiber_sld == 3) then

        !
        ! Orthotropic material
        !
        norma0(1,igaus) = 0.0_rp
        norma0(2,igaus) = 0.0_rp
        norma0(3,igaus) = 0.0_rp
        nshet0(1,igaus) = 0.0_rp
        nshet0(2,igaus) = 0.0_rp
        nshet0(3,igaus) = 0.0_rp
        do inode = 1,pnode
           norma0(1,igaus) = norma0(1,igaus) + gpsha(inode,igaus) * elfin(1,inode)
           norma0(2,igaus) = norma0(2,igaus) + gpsha(inode,igaus) * elfin(2,inode)
           norma0(3,igaus) = norma0(3,igaus) + gpsha(inode,igaus) * elfin(3,inode)
           nshet0(1,igaus) = nshet0(1,igaus) + gpsha(inode,igaus) * elfis(1,inode)
           nshet0(2,igaus) = nshet0(2,igaus) + gpsha(inode,igaus) * elfis(2,inode)
           nshet0(3,igaus) = nshet0(3,igaus) + gpsha(inode,igaus) * elfis(3,inode)
        end do

     end if

     ! * * * *
     !INVARIANTS
     ! * * * *

     !I_1=Trace(C)
     i1=0.0_rp
     do idime = 1,ndime
        i1 = i1 + gpcau(idime,idime,pgaus)
     end do

     i4f = 0.0_rp
     i4s = 0.0_rp
     i4n = 0.0_rp
     i8fs = 0.0_rp


     do idime = 1,ndime
        do jdime = 1,ndime
           i4f  =  i4f + nfibe0(idime,igaus) * gpcau(idime,jdime,igaus) * nfibe0(jdime,igaus) ! I_4f  = f0 * C_ij * f0 (eq. 5.1)
           i4s  =  i4s + nshet0(idime,igaus) * gpcau(idime,jdime,igaus) * nshet0(jdime,igaus) ! I_4s  = s0 * C_ij * s0
           i4n  =  i4n + norma0(idime,igaus) * gpcau(idime,jdime,igaus) * norma0(jdime,igaus) ! I_4n  = n0 * C_ij * n0
           i8fs = i8fs + nfibe0(idime,igaus) * gpcau(idime,jdime,igaus) * nshet0(jdime,igaus) ! I_8fs = f0 * C_ij * s0 (eq. 5.3)
           ! * * * * *
           !Transform f=F'*f0 ; s=F'*s0 ; n=F'*n0 (PUSH FORWARD)

           nfibe(idime,igaus) = nfibe(idime,igaus) + gpgdi(idime,jdime,igaus) * nfibe0(jdime,igaus)
           nshet(idime,igaus) = nshet(idime,igaus) + gpgdi(idime,jdime,igaus) * nshet0(jdime,igaus)
           norma(idime,igaus) = norma(idime,igaus) + gpgdi(idime,jdime,igaus) * norma0(jdime,igaus)
        end do
     end do

     ! Store nfibe for postproc

     do idime = 1,ndime
        gpfib(idime,igaus,ielem) = nfibe(idime,igaus)
     end do

     !
     ! f x f , s x s, f x s, s x f (outer products [3x3])
     !
     fporf = 0.0_rp
     spors = 0.0_rp
     fpors = 0.0_rp
     sporf = 0.0_rp
     do idime = 1,ndime
        do jdime = 1,ndime
           fporf(idime,jdime) = nfibe(idime,igaus) * nfibe(jdime,igaus)
           spors(idime,jdime) = nshet(idime,igaus) * nshet(jdime,igaus)
           fpors(idime,jdime) = nfibe(idime,igaus) * nshet(jdime,igaus)
           sporf(idime,jdime) = nshet(idime,igaus) * nfibe(jdime,igaus)
        end do
     end do



     ! * * * * * *
     ! Passive terms (calculate directy the cauchy stress (castr))
     ! * * * * * *

     ! Non-collageous material
     term1 = (a/gpdet(igaus))*exp(b*(i1-3.0_rp))-(a/gpdet(igaus))

     term2= 0.0_rp
     term3= 0.0_rp
     term4= 0.0_rp


     if (kfl_fiber_sld > 0) then

        !muscle fibers in tension
        if( i4f >= 1.0_rp ) then
           afm=af
        else  !muscle fibers in compression
           afm=6.5_rp*af
        end if

        ! these checks re necessary, because round-off errors from i4f/s can create spurious small forces in the passive 
        !if (abs((i4f-1.0_rp)) < 1.0e-12) i4f = 1.0_rp
        !if (abs((i4s-1.0_rp)) < 1.0e-12) i4s = 1.0_rp

        term2 = (2.0_rp/gpdet(igaus))*afm*(i4f-1.0_rp)*exp(bf *(i4f-1.0_rp)**2.0_rp)

        !collagen fibers in tension
        if (i4s >= 1.0_rp) then
           term3 = (2.0_rp/gpdet(igaus))*as*(i4s-1.0_rp)*exp(bs*(i4s-1.0_rp)**2.0_rp)
        else
           term3 = (2.0_rp/gpdet(igaus))*as*(i4s-1.0_rp)*exp(bs*(i4s-1.0_rp)**2.0_rp)
        end if

        !shear
        term4 = (afs/gpdet(igaus))*i8fs*exp(bfs*i8fs**2.0_rp)
     end if

     do idime = 1,ndime
        do jdime = 1,ndime
           castr(idime,jdime) = castr(idime,jdime)   &
                + term1 * gpcal(idime,jdime)        &
                + term2 * fporf(idime,jdime)        &
                + term3 * spors(idime,jdime)        &
                + term4 * (fpors(idime,jdime)+sporf(idime,jdime))
        end do
     end do


 if (ielem==0 .and. igaus==1) then
     print*, ' '
     print*, 'castr'
     print*, castr(1,1), castr(1,2), castr(1,3)
     print*, castr(2,1), castr(2,2), castr(2,3)
     print*, castr(3,1), castr(3,2), castr(3,3)

     print*, ' '
     print*, 'TERM2'
     print*, (2.0_rp/gpdet(igaus))*afm*(i4f-1.0_rp)
     print*, exp(bf *(i4f-1.0_rp)**2.0_rp)
     print*, term2

end if

     !NEW - Scale the stress with a user defined factor 
     !  do idime = 1,ndime
     !     do jdime = 1,ndime
     !        castr(idime,jdime) = scalf*castr(idime,jdime)
     !     end do
     !  end do


     ! * * * * * *
     ! Volumetric term
     ! * * * * * *

     !term1 = gpdet(igaus)*Kct*(1.0_rp-(1.0_rp/gpdet(igaus)))    
     term1 = 1/gpdet(igaus)* Kct * gpmof(igaus)*log(gpdet(igaus))!*gpcai(igaus)  
     !term1 = 4.0_rp * Kct * ( gpdet(igaus)**3.0_rp - gpdet(igaus) )

     do idime = 1,ndime
        castv(idime,idime) = castv(idime,idime)+term1
     end do


     !NEW - Scale the stress with a user defined factor 
     !     do idime = 1,ndime
     !        do jdime = 1,ndime
     !           castv(idime,jdime) = scalf*castv(idime,jdime)
     !        end do
     !     end do

     !   !(Kroon, 2007)
     !
     !    term1=4.0_rp*Kct*(gpdet(igaus)**3.0_rp - gpdet(igaus))
     !
     !    do idime=1,ndime
     !         castv(idime,idime) = castv(idime,idime)+term1
     !    end do


     !Calculate the stretch (lamda). In the fiber direction: lambda^2=f_0*C*f_0 = I_4f

     lamda(1) = sqrt(i4f)
     lamda(2) = sqrt(i4s)
     lamda(3) = sqrt(i4n)

     ovlam(1) = 0.0_rp
     ovlam(2) = 0.0_rp
     ovlam(3) = 0.0_rp
     if( lamda(1) /= 0.0_rp ) then
        !         gplep(1,1,igaus) = log(lamda(1))
        ovlam(1) = 1.0_rp / lamda(1)
     end if
     if( lamda(2) /= 0.0_rp ) then
        !         gplep(2,2,igaus) = log(lamda(2))
        ovlam(2) = 1.0_rp / lamda(2)
     end if
     if( lamda(3) /= 0.0_rp ) then
        !         gplep(3,3,igaus) = log(lamda(3))
        ovlam(3) = 1.0_rp / lamda(3)
     end if

     !
     ! Assemble fiber stretch (lamda(1)), parallel final assembly done in 
     ! 
     do inode= 1,pnode
        ipoin = lnods(inode)
        stretlam(ipoin)= stretlam(ipoin) + lamda(1) * gpvol(igaus) * gpsha(inode,igaus) / vmass(ipoin)
     end do

     !
     ! Updated fiber direction: lambda_f * f_true = F * f_0 = f
     ! => here f is the UNIT fiber direction in the current configuration
     ! used for active stress and postproc
     !
     nfibe_length= 0.0_rp
     do idime = 1,ndime
        nfibt(idime,igaus) = nfibe(idime,igaus) * ovlam(1)
        nshtt(idime,igaus) = nshet(idime,igaus) * ovlam(2)
        normt(idime,igaus) = norma(idime,igaus) * ovlam(3)
        nfibe_length=   nfibe_length + nfibe(idime,igaus)*nfibe(idime,igaus) 
     end do
     nfibe_length= sqrt(nfibe_length) 

     ! Electromechanical coupling
     ! * * * * * *
     ! Active stress T(lambda, [Ca++])m x m  (Hunter p 688)
     ! * * * * * *

     if(( coupling('SOLIDZ','EXMEDI') >= 1 ) .or. ( coupling('EXMEDI','SOLIDZ') >= 1 )) then

        gpcac= 0.0_rp
        do inode = 1,pnode
           gpcac= gpcac + gpsha(inode,igaus) * elcac(inode)   
        end do

        if((kfl_eccty(pmate) == 1 .or. kfl_eccty(pmate) == 3 .or. kfl_eccty(pmate)==4 ) .and. (ikact>0)) then
           ! hunter and land unidirectional/bidirectional models

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
     end if


     !if (ielem == 1253) write(695,*) nfibe_length,tactf
     !if (ielem ==    1) write(696,*) nfibe_length,tactf

     !Transform Cauchy stress to S
     do idime = 1,ndime
        do jdime = 1,ndime
           do kdime = 1,ndime
              do ldime = 1,ndime
                 gpstr(idime,jdime,igaus,1) = gpstr(idime,jdime,igaus,1)+ &
                      gpdet(igaus) *&
                      gpigd(idime,kdime,igaus) * &
                      (castr(kdime,ldime) + castv(kdime,ldime) + castr_active(kdime,ldime))*&
                      gpigd(jdime,ldime,igaus)
              end do
           end do
        end do
     end do


     !!Usik 2000
     !!
     !call invmtx(gpcau(1,1,igaus),gpcai(1,1),bidon,ndime)
     !term1=Kct*gpdet(igaus)*log(gpdet(igaus))
     !do idime = 1,ndime
     !    do jdime = 1,ndime
     !       gpstr(idime,jdime,igaus,1) = gpstr(idime,jdime,igaus,1) + &
     !                                  term1*gpcai(idime,jdime)
     !    end do
     !end do


     !do jdime=1,ndime
     !   do idime=1,ndime
     !      call invmtx(gpcau(idime,jdime,igaus),gpcai(idime,jdime,igaus),bidon,ndime)
     !   end do
     !end do
     !     call invmtx(gpcau(:,:,igaus),gpcai,bidon,ndime)



     if (flagt == 1_ip) then
        ! Tangent moduli d2W/dE_ij*dE_kl
        call invmtx(gpcau(:,:,igaus),gpcai,bidon,ndime)
        do idime = 1,ndime
           do jdime = 1,ndime
              do kdime = 1,ndime
                 do ldime=1,ndime
                    term1 = 2.0_rp*a*b*exp(b*(i1-3.0_rp))*tkron(idime,jdime)*tkron(kdime,ldime)
                    term2 = (8.0_rp*af*bf*exp(bf*(i4f-1.0_rp)**2.0_rp)*((i4f-1.0_rp)**2.0_rp)&
                         + 4.0_rp*af*exp(bf*(i4f-1.0_rp)**2.0_rp))*fporf(idime,jdime)*fporf(kdime,ldime)
                    term3 = (8.0_rp*as*bs*exp(bs*(i4s-1.0_rp)**2.0_rp)*((i4s-1.0_rp)**2.0_rp)&
                         + 4.0_rp*as*exp(bs*(i4s-1.0_rp)**2.0_rp))*spors(idime,jdime)*spors(kdime,ldime)
                    term4 = 4.0_rp*(afs+2.0_rp*afs*bfs*(i8fs**2.0_rp))*exp(bfs*(i8fs**2.0_rp))&
                         *(sporf(idime,jdime)+fpors(idime,jdime))*(sporf(kdime,ldime)+fpors(kdime,ldime)) 
                    !  termK=2.0_rp*Kct*((gpdet(igaus)**2.0_rp)-gpdet(igaus))*(-0.5_rp*(tkron(idime,kdime)*tkron(jdime,ldime)+tkron(idime,ldime)*tkron(jdime,kdime)))&
                    !        +2.0_rp*Kct*(gpdet(igaus)**2.0_rp-gpdet(igaus)*0.5_rp)*gpcai(idime,jdime)*gpcai(kdime,ldime)
                    termK = Kct*gpmof(igaus)*(gpcai(idime,jdime)*gpcai(kdime,ldime)&
                         -log(gpdet(igaus))*(gpcai(idime,kdime)*gpcai(jdime,ldime)+gpcai(idime,ldime)*gpcai(jdime,kdime)))
                    gpdds(idime,jdime,kdime,ldime,igaus) = &
                         gpdds(idime,jdime,kdime,ldime,igaus)+term1+term2+term3+term4+termK
                 end do
              end do
           end do
        end do

     end if

  end do !gauss points


  !100 format (5(E16.8,','))
  !101 format (6(F16.8,','))


end subroutine sld_stress_model_134

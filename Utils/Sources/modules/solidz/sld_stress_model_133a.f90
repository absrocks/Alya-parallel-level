subroutine sld_stress_model_133a(pgaus,pmate,gpgdi,gpigd,gpcau,gpstr,gpdet,gptlo,ielem,elcod,&
             gpigd_eps,gpgdi_eps,gpdet_eps,&
             gpmof)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_stress_model_133
  ! NAME
  !    sld__stress_model_133
  ! DESCRIPTION
  !    Law use by Lin & Yan, passive behavior. Do NOT include the EC coupling. 
  ! 
  !    IN:
  !    GPCAU ... Right Cauchy-Green deformation tensor ... C = F^t x F
  !    GPGDI ... Deformation tensor ...................... F = grad(u) + I    
  !
  !    OUT:
  !    GPSTR ... 2nd P-K Stress tensor ........................... S
  !    GPENE ... Stored energy function .................. W
  !
  !    LOCAL: 
  !    GPCAS....EQUAL TO C (GPCAU) WHEN NO DECOMPISITION IS USED (if multi=1) ....C
  !         ....EQUAL TO C_{iso} WHEN DECOMP. IS USED ....C_{iso}=C*J^(-2/3)
  !    
  !    GPCAV....EQUAL TO C (GPCAU) WHEN NO DECOMPISITION IS USED (if multi=1) ....C
  !         ....EQUAL TO C_{vol} WHEN DECOMP. IS USED ....C_{vol}=I*J^(2/3)
  !
  !
  !
  ! Use the multiplicative decomposition for the compressibility (if multi=1): (Chaves,vol.2, p. 62)
  !   F = F(iso)F(vol)
  !   C = C(iso)C(vol)
  !   J(iso),J(vol)
  !   Notation:
  !   C(iso) = gpcas 
  !   C(vol) = gpcav
  !
  ! USES
  ! USED BY
  !    sld_sld_stress_model_133
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode
  use def_solidz, only       :  lawco_sld,parco_sld,forcm_sld,densi_sld,velas_sld
  use def_master, only       :  cutim,ittim,postp,ittyp,ITASK_ENDRUN,kfl_modul,coupling

  implicit none
  integer(ip), intent(in)    :: pgaus,pmate,ielem
  real(rp),    intent(in)    :: gpcau(ndime,ndime,pgaus),gpdet(pgaus),gpgdi(ndime,ndime,pgaus),gpigd(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: gpstr(ndime,ndime,pgaus,2)
  integer(ip)                :: igaus,idime,jdime,kdime,multi
  real(rp)                   :: nfibe0(ndime,pgaus),norma0(ndime,pgaus),nshet0(ndime,pgaus),gpmof(pgaus),&
       nfibe(ndime,pgaus),norma(ndime,pgaus),nshet(ndime,pgaus),&
       nfibt(ndime,pgaus),normt(ndime,pgaus),nshtt(ndime,pgaus)                                
  real(rp)                   :: p1,p2,p3,p4,a0,a1,a2,a3,a4,a5
  real(rp)                   :: dwadc(ndime,ndime),dwpdc(ndime,ndime),&
       i1,i4f,i4s,i4n,nporn(ndime,ndime),term1,term2,qterm,dqpdc(ndime,ndime),&
       gpcai(ndime,ndime),bidon,bidon2,K,dwvdc(ndime,ndime),gpcav(ndime,ndime),&
       gpcas(ndime,ndime,pgaus),forcb,sdwvdc,sdwpdc,sdwadc,trace,lamda(3)
  real(rp),    intent(in)    :: gptlo(pgaus) 
  real(rp),    intent(in)    :: elcod(ndime,mnode)  
  real(rp),    intent(in)        :: gpgdi_eps(ndime,ndime,pgaus)           ! Displacement Gradient F for the perturbed state
  real(rp),    intent(in)        :: gpigd_eps(ndime,ndime,pgaus)           ! Displacement Gradient Inverse F^{-1} for the perturbed state
  real(rp),    intent(in)        :: gpdet_eps(pgaus)                       ! 

  real(rp)                   :: castr_active(ndime,ndime),gpstr_active(ndime,ndime)

  K = parco_sld(2,1) !0.001 N/cm² , values of Usyk, 2000

  p1=parco_sld(3,1) 
  p2=parco_sld(4,1)     
  p3=parco_sld(5,1)     
  p4=parco_sld(6,1)     


  gpstr=0.0_rp
  gpcas=0.0_rp
  multi = 0

  nfibe=0.0_rp 
  nshet=0.0_rp  
  norma=0.0_rp 

  nfibt=0.0_rp
  nshtt=0.0_rp 
  normt=0.0_rp 

  ! Compute sound velocity only at time-step 0 (initialization)
  if (ittim == 0_ip) then
     velas_sld(1,1) = sqrt(parco_sld(1,1)/densi_sld(1,pmate))
  end if

  do igaus=1,pgaus

     castr_active=0.0_rp
     gpstr_active=0.0_rp

     !orientation fiber option (for the prolate spheroid or other regular geometry)
     call sld_angleh(3,pgaus,igaus,ielem,elcod,nfibe0,norma0,nshet0)

     if (multi==1) then
        !C = C(iso)C(vol) = gpcas gpcav
        do idime=1,ndime
           do jdime=1,ndime
              gpcas(idime,jdime,igaus)=gpcau(idime,jdime,igaus)*gpdet(igaus)**(-2.0_rp/3.0_rp)
           end do
        end do

        gpcav=0.0_rp
        do idime=1,ndime
           gpcav(idime,idime)= gpdet(igaus)**(2.0_rp/3.0_rp)
        end do

     else if (multi==0) then

        !gpcas = gpcav = gpcau
        do idime=1,ndime
           do jdime=1,ndime
              gpcas(idime,jdime,igaus)=gpcau(idime,jdime,igaus)
           end do
        end do

        gpcav=0.0_rp
        do idime=1,ndime
           do jdime=1,ndime
              gpcav(idime,jdime)=gpcau(idime,jdime,igaus)
           end do
        end do

     end if


     !I_1=Trace(C) (isochoric)
     i1=0.0_rp
     do idime=1,ndime
        i1 = i1 + gpcas(idime,idime,igaus)
     end do

     !I_4f = f0 * C_ij * f0 (isochoric)
     i4f = 0.0_rp
     do idime=1,ndime
        do jdime=1,ndime
           i4f = i4f + nfibe0(idime,igaus)*gpcas(idime,jdime,igaus)*nfibe0(jdime,igaus)
        end do
     end do

     !I_4s = s0 * C_ij * s0
     i4s = 0.0_rp
     do idime=1,ndime
        do jdime=1,ndime
           i4s = i4s + nshet0(idime,igaus)*gpcas(idime,jdime,igaus)*nshet0(jdime,igaus)
        end do
     end do

     !I_4n = n0 * C_ij * n0
     i4n = 0.0_rp
     do idime=1,ndime
        do jdime=1,ndime
           i4n = i4n + norma0(idime,igaus)*gpcas(idime,jdime,igaus)*norma0(jdime,igaus)
        end do
     end do


     !N x N (outer product)
     nporn=0.0_rp
     do idime=1,ndime
        do jdime=1,ndime
           nporn(idime,jdime) = nfibe0(idime,igaus)*nfibe0(jdime,igaus)
        end do
     end do

     !Transform f=F·f0 ; s=F·s0 ; n=F·n0

     do idime=1,ndime
        do jdime=1,ndime
           nfibe(idime,igaus) = nfibe(idime,igaus) + gpgdi(idime,jdime,igaus)*nfibe0(jdime,igaus)
           nshet(idime,igaus) = nshet(idime,igaus) + gpgdi(idime,jdime,igaus)*nshet0(jdime,igaus)
           norma(idime,igaus) = norma(idime,igaus) + gpgdi(idime,jdime,igaus)*norma0(jdime,igaus)
        end do
     end do

     ! * * * * * *
     ! Passive terms d(W_p)/dC (Lin & Yin (1998) - eq(1))
     ! * * * * * *

     !dQ/dC
     dqpdc=0.0_rp

     term1 = p3*i1 - 3.0_rp*p3 + 2.0_rp*p4*i4f - 2.0_rp*p4
     do idime=1,ndime
        do jdime=1,ndime
           dqpdc(idime,jdime) = term1 * nporn(idime,jdime)
        end do
     end do

     term2 = 2.0_rp*p2*i1 - 6.0_rp*p2 + p3*i4f - p3
     do idime=1,ndime
        dqpdc(idime,idime)=dqpdc(idime,idime) + term2
     end do

     !dW_p/dC

     qterm = p2*(i1-3.0_rp)**2.0_rp + p3*(i1-3.0_rp)*(i4f-1.0_rp) + p4*(i4f-1.0_rp)**2.0_rp

     dwpdc=0.0_rp
     do idime=1,ndime
        do jdime=1,ndime
           dwpdc(idime,jdime) = p1*exp(qterm)*dqpdc(idime,jdime)
        end do
     end do


     ! * * * * * *
     ! Actif terms d(W_a)/dC (Lin & Yan eq(3))
     ! * * * * * *

     !EC coupling
     if(( coupling('SOLIDZ','EXMEDI') >= 1 ) .or. ( coupling('EXMEDI','SOLIDZ') >= 1 )) then

        !Calculate the stretch (lamda). In the fiber direction: lambda^2=f_0*C*f_0 = I_4f
        lamda(1)=sqrt(i4f)
        lamda(2)=sqrt(i4s)
        lamda(3)=sqrt(i4n)



        !Updated fiber direction: lambda_f * f_true = F * f_0 = f
        !=> here f is realy the fiber direction in the ref. conf.
        do idime=1,ndime
           nfibt(idime,igaus) = nfibe(idime,igaus)/lamda(1)
           nshtt(idime,igaus) = nshet(idime,igaus)/lamda(2)
           normt(idime,igaus) = norma(idime,igaus)/lamda(3)
        end do

        call sld_eccoup(&
             pgaus,pmate,igaus,gpdet,gptlo,lamda,nfibt,nshtt,normt,&
             gpigd,castr_active,gpstr_active,ielem)


     end if

     ! * * * * * *    
     ! Volumetric term d(W_v)/dC (Dorri et al. eq(8))
     ! * * * * * * 

     !(Usyk 2000)

     !C^{-1}
     gpcai=0.0_rp        
     call invmtx(gpcau(1,1,igaus),gpcai(1,1),bidon,ndime)  

     !K - See header of the routine. 
     term1=K*gpdet(igaus)*log(gpdet(igaus))*0.5_rp

     dwvdc=0.0_rp
     do idime=1,ndime
        do jdime=1,ndime  
           dwvdc(idime,jdime) = term1*gpcai(idime,jdime)    
        end do
     end do

     ! * * * * * *            
     ! Stress S = 2*(dW_a/dC) + 2*(dW_p/dC) + 2*(dW_v/dC)
     ! * * * * * *    

     do idime=1,ndime
        do jdime=1,ndime                                      

           gpstr(idime,jdime,igaus,1)=gpstr(idime,jdime,igaus,1)+&
                gpstr_active(idime,jdime)+&  
                2.0_rp*dwpdc(idime,jdime)+&  
                2.0_rp*dwvdc(idime,jdime)         
        end do
     end do


  end do !gauss points     


100 format (8(F16.8,','))
end subroutine sld_stress_model_133a

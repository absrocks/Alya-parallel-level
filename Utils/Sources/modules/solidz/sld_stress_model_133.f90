subroutine sld_stress_model_133(pgaus,pmate,gpcau,gpstr,gpdet,gptlo,ielem,elcod,&
             gpigd_eps,gpgdi_eps,gpdet_eps,&
             gpmof)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_stress_model_133
  ! NAME
  !    sld__stress_model_133
  ! DESCRIPTION
  !    Law use by Lin & Yan.
  ! 
  !    IN:
  !    GPCAU ... Right Cauchy-Green deformation tensor ... C = F^t x F
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
  !    sld_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode
  use def_solidz, only       :  lawco_sld,parco_sld,forcm_sld,densi_sld,velas_sld
  use def_master, only       :  cutim,ittim,postp,ittyp,ITASK_ENDRUN,coupling

  implicit none
  integer(ip), intent(in)    :: pgaus,pmate,ielem
  real(rp),    intent(in)    :: gpcau(ndime,ndime,pgaus),gpdet(pgaus),gpmof(pgaus),elcod(ndime,mnode)
  real(rp),    intent(out)   :: gpstr(ndime,ndime,pgaus,2)
  integer(ip)                :: igaus,idime,jdime,kdime,multi
  real(rp)                   :: nfibe(ndime,pgaus),alpha,alphr,trace,trac2,acnte,bcnte,ccnte
  real(rp)                   :: p1,p2,p3,p4,a0,a1,a2,a3,a4,a5
  real(rp)                   :: green(ndime,ndime),dwadc(ndime,ndime),dwpdc(ndime,ndime),&
       i1,i4,nporn(ndime,ndime),term1,term2,qterm,dqpdc(ndime,ndime),&
       gpcai(ndime,ndime),bidon,bidon2,K,dwvdc(ndime,ndime),gpcav(ndime,ndime),&
       gpcas(ndime,ndime,pgaus),forcb,sdwvdc,sdwpdc,sdwadc
  real(rp),    intent(in)    :: gptlo(pgaus) 
  real(rp),    intent(in)        :: gpgdi_eps(ndime,ndime,pgaus)           ! Displacement Gradient F for the perturbed state
  real(rp),    intent(in)        :: gpigd_eps(ndime,ndime,pgaus)           ! Displacement Gradient Inverse F^{-1} for the perturbed state
  real(rp),    intent(in)        :: gpdet_eps(pgaus)                       ! 

  K = parco_sld(2,1) !0.001 N/cm² , values of Usyk, 2000

  p1=parco_sld(3,1) 
  p2=parco_sld(4,1)     
  p3=parco_sld(5,1)     
  p4=parco_sld(6,1)     

  a0=parco_sld(7,1)
  a1=parco_sld(8,1)
  a2=parco_sld(9,1)  
  a3=parco_sld(10,1)           
  a4=parco_sld(11,1)           
  a5=parco_sld(12,1)           

  ! Compute sound velocity only at time-step 0 (initialization)
  !if (ittim == 0_ip) then
  !   velas_sld(1,1) = sqrt(parco_sld(1,1)/densi_sld(1,1))
  !end if

  K= 0.001_rp !N/cm² , values of Usyk, 2000

  gpstr=0.0_rp
  gpcas=0.0_rp
  multi = 0


  do igaus=1,pgaus

     !EC coupling 
     if(( coupling('SOLIDZ','EXMEDI') >= 1 ) .or. ( coupling('EXMEDI','SOLIDZ') >= 1 )) then

        ! PETERSON'S METHOD (find active force)
!!!        call sld_activf(gptlo,ielem,pgaus)    

        ! Modify the active parameters with the forces coming from the EC coupling 
!!!!        forcb = ((bridg_sld(2,igaus,ielem)+bridg_sld(4,igaus,ielem)))/forcm_sld  

        a0 = a0*forcb
        a1 = a1*forcb
        a2 = a2*forcb
        a3 = a3*forcb
        a4 = a4*forcb
        a5 = a5*forcb          
     else 
        a0 = 0.0_rp
        a1 = 0.0_rp
        a2 = 0.0_rp
        a3 = 0.0_rp
        a4 = 0.0_rp
        a5 = 0.0_rp                          
     end if


     !provisional N - must be unitar
     nfibe(1,igaus) = 1.0_rp
     nfibe(2,igaus) = 0.0_rp
     nfibe(3,igaus) = 0.0_rp 

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

     !I_4=N_i*C_ij*N_j (isochoric)
     i4 = 0.0_rp
     do idime=1,ndime
        do jdime=1,ndime
           i4 = i4 + nfibe(idime,igaus)*gpcas(idime,jdime,igaus)*nfibe(jdime,igaus)
        end do
     end do

     !N x N (outer product)  
     nporn=0.0_rp   
     do idime=1,ndime
        do jdime=1,ndime
           nporn(idime,jdime) = nfibe(idime,igaus)*nfibe(jdime,igaus)
        end do
     end do

     ! * * * * * *    
     ! Passive terms d(W_p)/dC (Lin & Yin (1998) - eq(1))
     ! * * * * * *  

     !dQ/dC 
     dqpdc=0.0_rp     

     term1 = p3*i1 - 3.0_rp*p3 + 2.0_rp*p4*i4 - 2.0_rp*p4  
     do idime=1,ndime
        do jdime=1,ndime                     
           dqpdc(idime,jdime) = term1 * nporn(idime,jdime)       
        end do
     end do

     term2 = 2.0_rp*p2*i1 - 6.0_rp*p2 + p3*i4 - p3     
     do idime=1,ndime
        dqpdc(idime,idime)=dqpdc(idime,idime) + term2     
     end do

     !dW_p/dC

     qterm = p2*(i1-3.0_rp)**2.0_rp + p3*(i1-3.0_rp)*(i4-1.0_rp) + p4*(i4-1.0_rp)**2.0_rp    

     dwpdc=0.0_rp  
     do idime=1,ndime
        do jdime=1,ndime      
           dwpdc(idime,jdime) = p1*exp(qterm)*dqpdc(idime,jdime)
        end do
     end do


     ! * * * * * *    
     ! Actif terms d(W_a)/dC (Lin & Yan eq(3))
     ! * * * * * * 

     !dW_a/dC
     dwadc=0.0_rp    

     term1 = a1*i1 - 3.0_rp*a1 + 2.0_rp*a3*i4 - 2.0_rp*a3 + a5   
     do idime=1,ndime
        do jdime=1,ndime                     
           dwadc(idime,jdime) = term1 * nporn(idime,jdime)       
        end do
     end do

     term2 = a1*i4 - a1 + 2.0_rp*a2*i1 - 6.0_rp*a2 + a4     
     do idime=1,ndime
        dwadc(idime,idime)=dwadc(idime,idime) + term2     
     end do


     ! * * * * * *    
     ! Volumetric term d(W_v)/dC (Dorri et al. eq(8))
     ! * * * * * * 

     !(Ruter et al.)

     !C^{-1}
     !gpcai=0.0_rp        
     !call invmtx(gpcav(1,1),gpcai(1,1),bidon,ndime)  

     !K=5000.0_rp!70000.0_rp!100000.0_rp !2000_rp!2000.0_rp
     !term1=K*(gpdet( igaus)**2.0_rp-gpdet(igaus))
     !!term1=K*log(gpdet(igaus))

     !dwvdc=0.0_rp
     !do idime=1,ndime
     !   do jdime=1,ndime  
     !      dwvdc(idime,jdime) = term1*gpcai(idime,jdime)    
     !   end do
     !end do


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




     !OJO - cutre: Ecrire la force dans un fichier txt - Pour 2 sets seulement
     !
     !
     !
     sdwadc=dwadc(1,1)+dwadc(2,2)+dwadc(3,3)+dwadc(1,2)+dwadc(1,3)+dwadc(2,3)+dwadc(2,1)+dwadc(3,1)+dwadc(3,2)       
     sdwpdc=dwpdc(1,1)+dwpdc(2,2)+dwpdc(3,3)+dwpdc(1,2)+dwpdc(1,3)+dwpdc(2,3)+dwpdc(2,1)+dwpdc(3,1)+dwpdc(3,2)    
     sdwvdc=dwvdc(1,1)+dwvdc(2,2)+dwvdc(3,3)+dwvdc(1,2)+dwvdc(1,3)+dwvdc(2,3)+dwvdc(2,1)+dwvdc(3,1)+dwvdc(3,2)    
     if (ielem==1099) then
        if (igaus==1) then


           bidon2=real(mod(ittim,postp(1) % npp_stepi(1)),rp)  !write every "postp(1) % npp_stepi(1)" itt
           if (bidon2 == 0.0_rp .or. ittyp == ITASK_ENDRUN) then 
              write(6333,100) cutim,forcb,sdwadc,sdwpdc,sdwvdc 
           end if
        end if
     end if





     ! * * * * * *            
     ! Stress S = 2*(dW_a/dC) + 2*(dW_p/dC) + 2*(dW_v/dC)
     ! * * * * * *    


     do idime=1,ndime
        do jdime=1,ndime                                      
           !gpstr(idime,jdime,igaus,1)=gpstr(idime,jdime,igaus,1) +&
           !     2.0_rp*dwvdc(idime,jdime)+&
           !     2.0_rp*dwpdc(idime,jdime)      

           !gpstr(idime,jdime,igaus,1)=gpstr(idime,jdime,igaus,1) +&
           !                         2.0_rp*dwpdc(idime,jdime) +&
           !                         2.0_rp*dwadc(idime,jdime) +&                                  
           !                         2.0_rp*dwvdc(idime,jdime)

           gpstr(idime,jdime,igaus,1)=gpstr(idime,jdime,igaus,1)+&
                2.0_rp*dwadc(idime,jdime)+&  
                2.0_rp*dwpdc(idime,jdime)+&  
                2.0_rp*dwvdc(idime,jdime)         
        end do
     end do


  end do !gauss points     

  !write(*,*) 'gpstr ',gpstr
100 format (8(F16.8,','))
end subroutine sld_stress_model_133

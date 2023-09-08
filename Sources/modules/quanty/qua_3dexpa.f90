subroutine arma_pp_3D(nato,nespec,atnumber) 

  use def_master
  use def_domain
  use def_quanty

  implicit none
  integer(ip), intent(in) :: nato,nespec,atnumber
  integer(ip)             :: lmax,loc,nrad,nspin,idime,nk,nderiv,nconta,nopa(10)
  real(rp)                :: x0,y0,z0,epsil,ainte
  !real(rp), intent(in)    :: radio(nrad),ppseu(lmax+1,nrad),ppphi(lmax+1,nrad)
  complex(rp)             :: ylm,armonico,funteor,funteor1d
  integer(ip)             :: kk,jj,ll,ii,kpos,NP,ncm,nesta,ncp0,ncl0,ncm0,nl,nm,mm,ncore,narma_rho
  real(rp)                :: RAD_NP,APEND,YL,XB,FUNYLM,PI,vale,zz ,noocp,ANTE
  character(20)           :: filefin,num,fin,FILEVPS

  real(rp), allocatable   :: v_aux(:,:)
  complex(rp),allocatable :: phionda(:),phi_ps(:,:,:),auxi_lm(:)

  !                                        _l=lmax  _m=l   _
  !   problema rd ==>  Vps(x,y,z) = V0(r) +>        >      Ylm(x,y,z)*(V0(r)-Vl(r))*Ylm(x,y,z)
  !                                        -l=1     -m=-l 
  !    
  ! En 3d necesito la funcion que dado x y z me de Ylm(tita,phi) en el espacio.


  x0 = atomo_qua(1)%coorx(nato)
  y0 = atomo_qua(1)%coory(nato)
  z0 = atomo_qua(1)%coorz(nato)

  lmax = especie_qua(nespec)%nlmax
  loc  = especie_qua(nespec)%nlocal
  nrad = especie_qua(nespec)%nrad
  vale =especie_qua(nespec)%valencia
  nspin= atomo_qua(1)%spin(nato)
  zz = atnumber 

  PI=4.0_rp*atan(1.0_rp)
  NP=npoin

  ! especie_qua(kato)%ppphi(ii,kk)
  narma_rho = 1


  ALLOCATE(phionda(NP),V_AUX(LMAX+1,NP),phi_ps(LMAX+1,2*(LMAX)+1,np),auxi_lm(np))


  DO JJ=1,LMAX+1

     DO KK=1,NP

        rad_np=sqrt( (coord(1,kk)-x0)**2 + (coord(2,kk)-y0)**2 + (coord(3,kk)-z0)**2) 

        KPOS=1
        DO WHILE( KPOS<= nrad .and. especie_qua(nespec)%radio(KPOS)<RAD_NP )
           KPOS=KPOS+1
        ENDDO
        IF(KPOS.EQ.1) THEN

           V_AUX(JJ,KK) = especie_qua(nespec)%ppseu(JJ,KPOS)

        ELSEIF(KPOS > 1 .and. KPOS <= nrad) THEN

           APEND = (especie_qua(nespec)%ppseu(JJ,KPOS)-especie_qua(nespec)%ppseu(JJ,KPOS-1))/(especie_qua(nespec)%radio(KPOS)-especie_qua(nespec)%radio(KPOS-1))
           XB    = especie_qua(nespec)%ppseu(JJ,KPOS) - APEND*especie_qua(nespec)%radio(KPOS)

           V_AUX(JJ,KK) =  rad_np * APEND  + XB

        ELSEIF(KPOS > nrad) THEN

           V_AUX(JJ,KK)=especie_qua(nespec)%ppseu(JJ,nrad)    

        ENDIF

     ENDDO

  ENDDO

  !  AHORA CALCULO EL PSEUDOPOTENCIAL DE VERDAD 

  write(num,'(i1)') NATO
  FILEVPS ='PSPOT'//NUM(1:1)//'.3D'  

  !OPEN(UNIT=111, FILE=FILEVPS,STATUS='UNKNOWN')


  DO KK=1,NP
     v_pot_ps(KK)=v_pot_ps(KK) + V_AUX(LOC,KK)
  ENDDO



  if(kfl_spher==1) then
     !write(111,*)  ' x  V'
     !DO KK=1,NP
     !   WRITE(111,'(2E15.5)') coord(1,KK),v_pot_ps(KK)
     !ENDDO

  else
     !write(111,*)  ' x     y      z      V'
     !DO KK=1,NP
     !   WRITE(111,'(4E15.5)') coord(1,KK),coord(2,KK),coord(3,KK),v_pot_ps(KK)
     !ENDDO
  endif

  !CLOSE(111)


  ! USO loca para crear la RHOON inicial
  ! creo que debo sumar valencia con los estados / 2
  call qua_idespe(atnumber,nspin,ncp0,ncl0,ncm0) 

  ncore=(atnumber-int(vale,ip))/nspin

  nesta = 0
  do nk=1+ncore,ncp0 
     do nl=0,nk-1
        if(nl.le.ncl0) then
           do mm=1,2*nl+1
              nesta=nesta + 1  
           enddo
        endif
     enddo
  enddo

  if(narma_rho.eq.0) then ! arma el rho con LOCAS

     noocp=0
     do nk=ncp0,1,-1 

        do nl=0,nk-1

           do mm=1,2*nl+1

              nm = mm - (nl+1) 

              if(nesta.le.0) then
                 noocp = 100
              else
                 noocp = noocp + especie_qua(nespec)%nocupa(nesta)
              endif

              if(noocp.le.vale*1.05) then
                 do kk=1,npoin

                    if(kfl_spher==1) then
                       ylm = 1.0/sqrt(4*PI)
                       phionda(kk)=funteor1d(NK,NL,NM,ZZ,coord(1,kk),X0) * YLM

                    else

                       ylm=armonico(coord(1,kk)-x0,coord(2,kk)-y0,coord(3,kk)-z0,NL,NM)

                       phionda(kk)=funteor(NK,NL,NM,ZZ,coord(1,kk),coord(2,kk), &
                            coord(3,kk),X0,Y0,Z0) * YLM

                    endif
                    rhoon(kk,1)=rhoon(kk,1) + especie_qua(nespec)%nocupa(nesta)*conjg(phionda(kk))*phionda(kk) 

                 enddo

              endif

              nesta = nesta - 1  
           enddo
        enddo
     enddo

  else ! lo arma con psusofunciones

     ! primero las expando!

     !especie_qua(nespec)%nocupa(nesta)

     nconta=0

     DO JJ=1,LMAX+1
        LL=JJ-1

        DO II=1,2*LL+1

           nconta=nconta+1     

           NCM = -LL + II - 1 


           DO KK=1,NP

              rad_np=sqrt( (coord(1,kk)-x0)**2 + (coord(2,kk)-y0)**2 + (coord(3,kk)-z0)**2) 

              KPOS=1
              DO WHILE( KPOS<= nrad .and. especie_qua(nespec)%radio(KPOS)<RAD_NP )
                 KPOS=KPOS+1
              ENDDO

              ylm=armonico(coord(1,kk)-x0,coord(2,kk)-y0,coord(3,kk)-z0,NL,NM)

              IF(KPOS.EQ.1) THEN

                 phi_ps(jj,ii,kk) = especie_qua(nespec)%ppphi(JJ,KPOS) * ylm

              ELSEIF(KPOS > 1 .and. KPOS < nrad) THEN

                 APEND = (especie_qua(nespec)%ppphi(JJ,KPOS)-especie_qua(nespec)%ppphi(JJ,KPOS-1))/(especie_qua(nespec)%radio(KPOS)-especie_qua(nespec)%radio(KPOS-1))
                 XB    = especie_qua(nespec)%ppphi(JJ,KPOS) - APEND*especie_qua(nespec)%radio(KPOS)

                 phi_ps(jj,ii,kk) =  (rad_np * APEND  + XB) * ylm

              ELSEIF(KPOS > nrad) THEN

                 phi_ps(jj,ii,kk) = especie_qua(nespec)%ppphi(JJ,nrad) *ylm   

              ENDIF

              rhoon(kk,1)=rhoon(kk,1) + especie_qua(nespec)%nocupa(nconta)*conjg(phi_ps(jj,ii,kk))*phi_ps(jj,ii,kk) 

           ENDDO

           ! verifico las seudo-funciones de onda
           do kk=1,np
              auxi_lm(kk)=CONJG(phi_ps(jj,ii,kk))*phi_ps(jj,ii,kk)
           enddo
           NK=0
           NDERIV=0
           call  qua_3dinte(NK,NDERIV,auxi_lm,AINTE) 

           ! esto me da el denom para KB

           DO KK=1,NP

              auxi_lm(kk)=CONJG(phi_ps(JJ,II,KK))*phi_ps(JJ,II,KK)*(V_AUX(JJ,KK)-V_AUX(Loc,KK))

           enddo

           NK=0
           call  qua_3dinte(NK,NDERIV,auxi_lm,AINTE) 
           NL_denom(jj,ii)=NL_denom(jj,ii)+(AINTE)             

           !  construyo los vectores NO locales del PP con KB
           if(NL_denom(jj,ii).ne.0.0) call kb(jj,ii,V_AUX,lmax+1,phi_ps,loc) 

        ENDDO
     ENDDO

  endif


  DEALLOCATE(V_AUX,phionda,phi_ps,auxi_lm)

end subroutine arma_pp_3D



subroutine kb(ll,mm,V_AUX,lmm,phi,loc) 

  use def_parame
  use def_master
  use def_domain
  use def_quanty
  implicit none

  integer(ip),intent(in) :: ll,mm,lmm,loc
  real(rp),intent(in)    :: V_aux(lmm,*)
  complex(rp),intent(in) :: phi(lmm,2*(lmm-1)+1,*)
  integer(ip)            :: inode,jnode   

  real(rp)    :: elmat(mnode,mnode)
  integer(ip) :: ielem,igaus                           ! Indices and dimensions
  integer(ip) :: pelty,pmate,pnode
  integer(ip) :: pgaus,plapl,porde,ptopo,ipoin

  real(rp)    :: elpro(mnode)                          ! Gather 
  real(rp)    :: elcod(ndime,mnode)

  real(rp)    :: gpvol(mgaus)                          ! |J|*w
  real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dNk/dxj
  real(rp)    :: gphes(ntens,mnode,mgaus)              ! dNk/dxidxj

  real(rp)    :: gpcon(mgaus)                          ! k
  real(rp)    :: gposc(mgaus)
  real(rp)    :: gpcod(ndime,mgaus)
  complex(rp) :: Har_pro,fact2
  complex(rp) :: gprhs(mgaus),elrhs(mnode)         ! f (all terms)

  !
  ! Loop over elements
  !
  do ielem=1,nelem
     !
     ! Element dimensions
     !
     pelty = ltype(ielem)
     pnode = nnode(pelty)
     pgaus = ngaus(pelty)
     plapl = llapl(pelty)
     porde = lorde(pelty)
     ptopo = ltopo(pelty)
     !
     ! Gather operations
     !
     call qua_elmgat(&
          pnode,lnods(1,ielem),elpro,elcod)
     !
     ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
     !
     call elmcar(&
          pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
          elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
          gphes,ielem)

     Har_pro=0.0

     do inode=1,pnode
        ipoin=lnods(inode,ielem)
        elrhs(inode)=0.0
        Har_pro = Har_pro + (V_AUX(ll,ipoin)-V_AUX(loc,ipoin))*PHI(ll,mm,ipoin)
     enddo
     Har_pro = Har_pro/pnode

     do igaus=1,pgaus
        gprhs(igaus)=Har_pro
     end do

     do igaus=1,pgaus
        fact2=gprhs(igaus)*gpvol(igaus)
        do inode=1,pnode
           elrhs(inode)=elrhs(inode)+fact2**elmar(pelty)%shape(inode,igaus)
        end do
     end do

     do inode=1,pnode           
        ipoin=lnods(inode,ielem)
        NonLoc(ll,mm,ipoin)=NonLoc(ll,mm,ipoin)+elrhs(inode)
     end do

  enddo


  RETURN
END subroutine kb


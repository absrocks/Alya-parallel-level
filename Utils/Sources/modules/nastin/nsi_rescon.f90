subroutine nsi_rescon()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_rescon()
  ! NAME 
  !    nsi_rescont
  ! DESCRIPTION
  !    This routine computed the continuous residuals
  ! USES
  ! USED BY
  !    nsi_cvgunk
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_nastin
  implicit none

  if( solve(ivari_nsi) % kfl_algso /= 0 ) then
     !
     ! Sparse assembly
     !
     call nsi_bcsrax(1_ip,ndime+1,amatr,rhsid)  ! Mom.+Cont.

  else
     !
     ! Profile assembly: not supported
     !
     resin_nsi(1) = 0.0_rp  ! Momentum residual
     resin_nsi(2) = 0.0_rp  ! Continuity residual     
     resss_nsi(1) = 0.0_rp  ! Steady momentum residual
     resss_nsi(2) = 0.0_rp  ! Steady continuity residual     
  end if

end subroutine nsi_rescon

subroutine nsi_bcsrax(itask,nbvar,amatr,rhsid) 
  !----------------------------------------------------------------------
  !****f* Nastin/nsi_bcsrax
  ! NAME 
  !     nsi_bcsrax
  ! DESCRIPTION
  !     The RHS (RHSID) be exchanged at the present residual is computed
  !     before calling the solver.
  ! INPUT
  !    NBVAR ....... Number of varr_dombles
  !    AMATR ....... Matrix
  !    RHSID ....... Right-hand side
  ! OUTPUT
  !    RESIN_NSI ... Continuous residual
  ! USES
  ! USED BY
  !***
  !----------------------------------------------------------------------
  use def_kintyp, only    :  ip,rp
  use def_master, only    :  kfl_paral,NPOIN_REAL_12DI,&
       &                     NPOIN_REAL_1DIM,parr1,IMASTER,&
       &                     veloc,press,nparr,parre,mem_modul,modul,&
       &                     unkno,densi,INOTMASTER,solve
  use def_domain, only    :  npoin,ndime,&
       &                     c_dom,r_dom,c_sym,r_sym
  use def_nastin, only    :  resin_nsi,&
       &                     kfl_trres_nsi,kfl_timei_nsi,&
       &                     kfl_regim_nsi,ndbgs_nsi,&
       &                     ncomp_nsi,resss_nsi
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: itask,nbvar
  real(rp),    intent(in) :: amatr(nbvar,nbvar,*)
  real(rp),    intent(in) :: rhsid(nbvar,npoin)
  integer(ip)             :: ii,jj,kk,ll,col,idime,ipoin,idofv,icomp
  integer(4)              :: istat
  real(rp),    pointer    :: yy(:,:),zz(:),bb2(:,:),bb1(:)
  real(rp)                :: raux
  real(rp),    pointer    :: unkn2(:,:)

  if(kfl_regim_nsi==2) then
     unkn2 => densi           ! 2nd variable is density
  else
     unkn2 => press           ! 2nd variable is pressure
  end if
  icomp = min(ncomp_nsi,3_ip) ! Last time step component

  select case(itask)

  case(1_ip)
     !
     ! Momentum and continuity equations
     !
     call nsi_rotunk(1_ip,unkno)  ! global to local

     if( IMASTER ) then
        allocate(yy(1,1),stat=istat)
     else 
        allocate(yy(ndime,npoin),stat=istat)
        call memchk(0_ip,istat,mem_modul(1:2,modul),'YY','nsi_rescon',yy)

     end if

     if( .not. IMASTER ) then
        do ii=1,npoin
           do kk=1,ndime
              yy(kk,ii) = 0.0_rp
           end do
           do jj=r_dom(ii),r_dom(ii+1)-1
              col=c_dom(jj)
              idofv=(col-1)*nbvar
              do ll=1,ndime  
                 idofv=idofv+1
                 raux=unkno(idofv)
                 do kk=1,ndime
                    yy(kk,ii)=yy(kk,ii)+amatr(ll,kk,jj)*raux
                 end do
              end do
              raux=unkno(col*nbvar)
              do kk=1,ndime
                 yy(kk,ii)=yy(kk,ii)+amatr(nbvar,kk,jj)*raux  
              end do
           end do
        end do
     end if
     call rhsmod(ndime,yy)
     call residu(2_ip,ndime+1_ip,ndime,rhsid,yy,1_ip,1_ip,ndime,1.0_rp,resin_nsi(1))

     if( .not. IMASTER ) then
        do ii=1,npoin
           do kk=1,ndime
              yy(kk,ii) = 0.0_rp
           end do
           do jj=r_dom(ii),r_dom(ii+1)-1
              col=c_dom(jj)
              idofv=(col-1)*nbvar
              do ll=1,ndime  
                 idofv=idofv+1
                 raux=unkno(idofv)
                 do kk=1,ndime
                    yy(kk,ii)=yy(kk,ii)+amatr(ll,kk,jj)*veloc(ll,col,icomp)
                 end do
              end do
              do kk=1,ndime
                 yy(kk,ii)=yy(kk,ii)+amatr(nbvar,kk,jj)*press(col,icomp)
              end do
           end do
        end do
     end if
     call rhsmod(ndime,yy)
     call residu(2_ip,ndime+1_ip,ndime,rhsid,yy,1_ip,1_ip,ndime,1.0_rp,resss_nsi(1))

     if( IMASTER ) then
        deallocate(yy,stat=istat)
     else
        call memchk(2_ip,istat,mem_modul(1:2,modul),'YY','nsi_rescon',yy)
        deallocate(yy,stat=istat)
        if(istat/=0) call memerr(2_ip,'YY','nsi_rescon',0_ip)
     end if
     !
     ! Continuity equation
     !
     if( IMASTER ) then
        allocate(zz(1),stat=istat)
     else
        allocate(zz(npoin),stat=istat)
        call memchk(0_ip,istat,mem_modul(1:2,modul),'ZZ','nsi_rescon',zz)
     end if

     if( INOTMASTER ) then
        do ii=1,npoin
           zz(ii) = 0.0_rp
           do jj=r_dom(ii),r_dom(ii+1)-1
              col=c_dom(jj)
              idofv=(col-1)*nbvar
              do ll=1,ndime
                 idofv=idofv+1
                 zz(ii)=zz(ii)+amatr(ll,nbvar,jj)*unkno(idofv)
              end do
              idofv=col*nbvar
              zz(ii)=zz(ii)+amatr(nbvar,nbvar,jj)*unkno(idofv)
           end do
        end do
     end if
     call rhsmod(1_ip,zz)
     call residu(2_ip,ndime+1,1_ip,rhsid,zz,ndime+1,1_ip,1_ip,1.0_rp,resin_nsi(2))

     if( INOTMASTER ) then
        do ii=1,npoin
           zz(ii) = 0.0_rp
           do jj=r_dom(ii),r_dom(ii+1)-1
              col=c_dom(jj)
              idofv=(col-1)*nbvar
              do ll=1,ndime
                 idofv=idofv+1
                 zz(ii)=zz(ii)+amatr(ll,nbvar,jj)*veloc(ll,col,icomp)
              end do
              zz(ii)=zz(ii)+amatr(nbvar,nbvar,jj)*press(col,icomp)
           end do
        end do
     end if
     call rhsmod(1_ip,zz)
     call residu(2_ip,ndime+1,1_ip,rhsid,zz,ndime+1,1_ip,1_ip,1.0_rp,resss_nsi(2))

     if( IMASTER ) then
        deallocate(zz,stat=istat)
     else
        call memchk(2_ip,istat,mem_modul(1:2,modul),'ZZ','nsi_rescon',zz)
        deallocate(zz,stat=istat)
        if(istat/=0) call memerr(2_ip,'ZZ','nsi_rescon',0_ip)
     end if
 
     call nsi_rotunk(2_ip,unkno)  ! local to global

  case(2_ip)
     !
     ! Momentum equation
     !
     if( IMASTER ) then
        allocate(yy(1,1), stat=istat)
        allocate(bb2(1,1),stat=istat)
     else
        allocate(yy(ndime,npoin),stat=istat)
        call memchk(0_ip,istat,mem_modul(1:2,modul),'YY', 'nsi_rescon',yy)
     end if

     if ( INOTMASTER ) then
        call nsi_rotunk(1_ip,unkno)   ! global to local
        call bcsrax(1_ip,npoin,ndime,amatr,c_dom,r_dom,unkno,yy)
        call nsi_rotunk(2_ip,unkno)   ! local to global
        allocate(bb2(ndime,npoin),stat=istat)
        call memchk(0_ip,istat,mem_modul(1:2,modul),'BB2','nsi_rescon',bb2)
        do ipoin=1,npoin
           do idime=1,ndime
              bb2(idime,ipoin)=rhsid(idime,ipoin)
           end do
        end do    
     end if
     call rhsmod(ndime,bb2)
     call residu(2_ip,ndime,ndime,bb2,yy,1_ip,1_ip,ndime,1.0_rp,resin_nsi(1))

     if ( INOTMASTER ) then
        call nsi_rotunk(3_ip,veloc(1,1,3))  ! global to local
        call bcsrax(1_ip,npoin,ndime,amatr,c_dom,r_dom,veloc(1,1,icomp),yy)
        call nsi_rotunk(4_ip,veloc(1,1,3))  ! local to global
     end if
     call residu(2_ip,ndime,ndime,bb2,yy,1_ip,1_ip,ndime,1.0_rp,resss_nsi(1))
 
     if( IMASTER ) then
        deallocate(bb2,stat=istat)
        deallocate(yy,stat=istat)
     else
        call memchk(2_ip,istat,mem_modul(1:2,modul),'BB2','nsi_rescon',bb2)
        deallocate(bb2,stat=istat)
        if(istat/=0) call memerr(2_ip,'BB2','nsi_rescon',0_ip)
        call memchk(2_ip,istat,mem_modul(1:2,modul),'YY','nsi_rescon',yy)
        deallocate(yy,stat=istat)
        if(istat/=0) call memerr(2_ip,'YY','nsi_rescon',0_ip)
     end if
 
  case(3_ip)
     !
     ! Continuity equation
     !
     if( IMASTER ) then
        allocate(zz(1), stat=istat)
        allocate(bb1(1),stat=istat)
     else
        allocate(zz(npoin),stat=istat)
        call memchk(0_ip,istat,mem_modul(1:2,modul),'ZZ','nsi_rescon',zz)
     end if

     if( INOTMASTER ) then
        if(solve(2)%kfl_symme==1) then
           call bsymax(1_ip,npoin,1_ip,amatr,c_sym,r_sym,press,zz)
        else
           call bcsrax(1_ip,npoin,1_ip,amatr,c_dom,r_dom,press,zz)
        end if
        allocate(bb1(npoin),stat=istat)
        call memchk(0_ip,istat,mem_modul(1:2,modul),'BB1','nsi_rescon',bb1)
        do ipoin=1,npoin
           bb1(ipoin)=rhsid(1,ipoin)
        end do  
     end if
     call rhsmod(1_ip,bb1)
     call residu(2_ip,1_ip,1_ip,bb1,zz,1_ip,1_ip,1_ip,1.0_rp,resin_nsi(2))

     if( INOTMASTER ) then
        if(solve(2)%kfl_symme==1) then
           call bsymax(1_ip,npoin,1_ip,amatr,c_sym,r_sym,press(1,icomp),zz)
        else
           call bcsrax(1_ip,npoin,1_ip,amatr,c_dom,r_dom,press(1,icomp),zz)
        end if
     end if
     call residu(2_ip,1_ip,1_ip,bb1,zz,1_ip,1_ip,1_ip,1.0_rp,resss_nsi(2))

     if( IMASTER ) then
        deallocate(bb1,stat=istat)
        deallocate(zz,stat=istat)
     else
        call memchk(2_ip,istat,mem_modul(1:2,modul),'BB1','nsi_rescon',bb1)
        deallocate(bb1,stat=istat)
        if(istat/=0) call memerr(2_ip,'BB1','nsi_rescon',0_ip)
        call memchk(2_ip,istat,mem_modul(1:2,modul),'ZZ','nsi_rescon',zz)
        deallocate(zz,stat=istat)
        if(istat/=0) call memerr(2_ip,'ZZ','nsi_rescon',0_ip)
     end if

  case(4_ip)
     !
     ! Momentum+continuity equations
     !
     if( IMASTER ) then
        allocate(yy(1,1),stat=istat)
     else 
        allocate(yy(nbvar,npoin),stat=istat)
        call memchk(0_ip,istat,mem_modul(1:2,modul),'YY','nsi_rescon',yy)
        do ii=1,npoin
           do kk=1,nbvar
              yy(kk,ii) = 0.0_rp
           end do
           do jj=r_dom(ii),r_dom(ii+1)-1
              col=c_dom(jj)
              do kk=1,ndime
                 yy(kk,ii)=yy(kk,ii)+amatr(nbvar,kk,jj)*press(col,1)             
                 do ll=1,ndime
                    yy(kk,ii)=yy(kk,ii)+amatr(ll,kk,jj)*veloc(ll,col,1)
                 end do
              end do
              kk=ndime+1
              yy(kk,ii)=yy(kk,ii)+amatr(nbvar,nbvar,jj)*press(col,1)              
              do ll=1,ndime
                 yy(kk,ii)=yy(kk,ii)+amatr(ll,nbvar,jj)*veloc(ll,col,1)
              end do
           end do
        end do
     end if
     call bcsrax(1_ip,npoin,ndime+1,amatr,c_dom,r_dom,unkno,yy)
     call rhsmod(nbvar,yy)
     call residu(2_ip,nbvar,nbvar,rhsid,yy,1_ip,1_ip,nbvar,1.0_rp,resin_nsi(2))
     resin_nsi(1)=resin_nsi(2)
     if( IMASTER ) then
        deallocate(yy,stat=istat)
     else
        call memchk(2_ip,istat,mem_modul(1:2,modul),'YY','nsi_rescon',yy)
        deallocate(yy,stat=istat)
        if(istat/=0) call memerr(2_ip,'YY','nsi_rescon',0_ip)
     end if

  end select

end subroutine nsi_bcsrax


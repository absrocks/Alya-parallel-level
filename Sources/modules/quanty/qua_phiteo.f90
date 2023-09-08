subroutine all_electron(nato,nspin,atnumber,x0,y0,z0,atommax) 
  use def_parame
  use def_master
  use def_domain
  use def_quanty
  use def_solver
  use mod_postpr
  implicit none
  complex(rp)  :: ylm,armonico,funteor,funteor1d
  integer(ip), intent(in) :: nato,atnumber,nspin,atommax
  real(rp),    intent(in) :: x0,y0,z0
  integer(ip)             :: memo,kk,nk,jj,ncp0,ncl0,ncm0,nl,nm,mm,nespe,nesta,nderiv
  real(rp)                :: zz,epsil,intphi
  character(25) :: filefin,num,fin,FILEVPS
  complex(rp), allocatable :: phionda(:)
character(5) :: wopos(2)


  zz=atnumber
  !PI=4.0_rp*atan(1.0_rp)

  call qua_idespe(atnumber,nspin,ncp0,ncl0,ncm0)   

  allocate(phionda(npoin)) 

  !write(num,'(i1)') NATO
  !FILEVPS ='VPOT'//NUM(1:1)//'.3D'  

  !if(noutput/=0.and.atommax==nato) then
  !  OPEN(UNIT=111, FILE=FILEVPS,STATUS='UNKNOWN')
  !  write(111,*)  ' x     y      z      Vpot'
  !endif

  if(kfl_spher==1) then

     do nk=1,npoin
        epsil = 0.0
        if(abs(x0-coord(1,nk)).lt.1e-9) epsil = 1e-9
        v_pot_ps(nk)=v_pot_ps(nk) - zz/sqrt( (coord(1,nk)-(x0+epsil))**2 )
        if(noutput/=0.and.atommax==nato) then
           !WRITE(111,'(2E15.5)') coord(1,nk),v_pot_ps(nk)
        endif
     enddo

  else

     do nk=1,npoin
        epsil=0.0_rp
 
        if((coord(1,nk)-x0)==0.0_rp .and. (coord(2,nk)-y0)==0.0_rp .and. (coord(3,nk)-z0)==0.0_rp) epsil = 1.0e-2_rp

        v_pot_ps(nk) = v_pot_ps(nk) - zz/sqrt( (coord(1,nk)-(x0+epsil))**2 + (coord(2,nk)-y0)**2 + (coord(3,nk)-z0)**2 )
        if(noutput/=0.and.atommax==nato) then
           !WRITE(111,'(4E15.5)') coord(1,nk),coord(2,nk),coord(3,nk),v_pot_ps(nk)
        endif
     enddo

  endif
 
  if(noutput/=0.and.atommax==nato) then
     !CLOSE(111)
  endif
!wopos(1)='POTEN'
!wopos(2)='SCALA'
!call postpr(v_pot_ps,wopos,ittim,cutim)
!call runend('OPOPSWS')

  memo = 0 

  nespe= atomo_qua(1)%Tiespe(NATO)

  do nk=1,ncp0 

     do nl=0,nk-1

        if(nl.le.ncl0) then

           do mm=1,2*nl+1

              nm = mm - (nl+1) 
              memo = memo + 1  
              !write(num,'(4i1)') NATO,NK,NL,MM

              !fin='ARMONICO'
              !filefin=fin(1:8)//num(1:4)//'.3D'

              if(noutput/=0) then
                 !open(unit=111,file= filefin,status = 'unknown')
              endif

              if(kfl_spher==0) then
                 if(noutput/=0) then
                    !write(111,*) ' x     y     z     ylm '
                 endif
              endif

              do kk=1,npoin

                 if(kfl_spher==1) then
                    ylm = 1.0_rp/sqrt(4.0_rp*PI)
                    phionda(kk)=funteor1d(NK,NL,NM,ZZ,coord(1,kk),X0) * YLM

                    if(noutput/=0) then
                       !write(111,'(2E15.5)') coord(1,kk)-x0,ABS(YLM)
                    endif

                 else

                    ylm=armonico(coord(1,kk)-x0,coord(2,kk)-y0,coord(3,kk)-z0,NL,NM)

                    phionda(kk)=funteor(NK,NL,NM,ZZ,coord(1,kk),coord(2,kk), &
                         coord(3,kk),X0,Y0,Z0) * YLM

                    if(noutput/=0) then
                       !write(111,'(4E15.5)') coord(1,kk)-x0,coord(2,kk)-y0,coord(3,kk)-z0,ABS(YLM)
                    endif
                 endif

                 rhoon(kk,1) = rhoon(kk,1) + especie_qua(nespe)%nocupa(memo)*conjg(phionda(kk))*phionda(kk) 

              enddo

              if(noutput/=0) then
                 !close(111)
              endif


              fin='PHI'
              filefin = "data\\"//fin(1:3)//num(1:4)//'.3D'
              if(noutput/=0) then
                 !OPEN(UNIT=111, FILE=filefin,STATUS='UNKNOWN')
              endif

              if(kfl_spher==0) then

                 if(noutput/=0) then
                    !WRITE(111,*) 'x      y       z      phi'
                    DO KK=1,npoin
                       if(abs(PHIONDA(KK)) < 1E-20) phionda(kk) = (0.0_rp,0.0_rp)
                       !WRITE(111,'(4E15.8)') coord(1,KK),coord(2,KK),coord(3,KK),ABS(PHIONDA(KK)) 
                    ENDDO
                    nderiv=0
                    !print*,'all electron: integral=',nk,nl,mm
                    !call qua_3dinte(1_ip,NDERIV,PHIONDA,intphi)
                    !print*,'all electron: integral finished'
                    !write(6,*) 'integral de Phionda  ',intphi

                 endif

              else
                 if(noutput/=0) then
                    !WRITE(111,*) 'x phi'
                    DO KK=1,npoin
                       if(abs(PHIONDA(KK)) < 1E-20) phionda(kk)=(0.0,0.0)
                       !WRITE(111,'(2E15.8)') coord(1,KK),ABS(PHIONDA(KK)) 
                    ENDDO
                 endif
              endif

              if(noutput/=0) then
                 !close(111)
              endif
           enddo
        endif
     enddo
  enddo

  !if(noutput/=0) then
  !   fin='RHO'
  !   filefin = "data\\"//fin(1:3)//'.3D'
  !   OPEN(UNIT=111, FILE=filefin,STATUS='UNKNOWN')
  !   write(111,*)  ' x     y      z      Rho'
     !NK=0
     !NDERIV=0
     !call qua_3dinte_R(NK,NDERIV,RHOON(:,1),intphi)
     !print*,'Integral de rho=',intphi
  !   write(eigen_sol(1)%lun_solei ,*) 'integral de Rho  ',NK,intphi
  !   DO KK=1,npoin
  !      WRITE(111,'(4E15.8)') coord(1,KK),coord(2,KK),coord(3,KK),rhoon(KK,1) 
  !   ENDDO
  !   close(111)
  !endif

  deallocate(phionda)

  print*,'all electron finished)=',kfl_paral

end subroutine all_electron

subroutine saca_phi(NCP0,NCL0,NCM0,X0,Y0,Z0,PHI,NK)
  use def_master
  use def_domain
  implicit none
  integer(ip), intent(in) :: NCP0,NCL0,NCM0,nk
  real(rp), intent(in)    :: X0,Y0,Z0
  complex(rp),intent(in)  :: phi(*)
  integer(ip)                 :: kk
  character(20) filefin
  character(6)  num,fin

  fin='PHI'

  !write(num,'(4i1)') NK,NCP0,NCL0,NCM0

  !filefin = fin(1:3)//num(1:4)//'.3D'


  !OPEN(UNIT=111, FILE=filefin,STATUS='UNKNOWN')

  !WRITE(111,*) 'x      y       z      phi'
  !DO KK=1,npoin

     !WRITE(111,'(4E15.5)') coord(1,KK),coord(2,KK),coord(3,KK),ABS(PHI(KK)) ! SQRT(CONJG(PHION(KK))*PHION(KK))

  !ENDDO


  !close(111)

end subroutine saca_phi


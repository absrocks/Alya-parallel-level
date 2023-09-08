subroutine qua_poibou()
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_poibou
  ! NAME 
  !    qua_poibou
  ! DESCRIPTION
  !    This routine calculate the boundaries conditions for the Poisson ecuation  
  !    depending of the rho and cuadrupolar terms
  ! USES
  !    qua_integracuad
  ! USED BY
  !    qua_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_quanty, only : bvessH_qua,kfl_fixno_qua 
  use def_master
  use def_domain
  use def_quanty
  use mod_postpr
  implicit none
  real(rp)     :: cpu_refe1,cpu_refe2,rmedio,ante,vlim
  integer(ip)  :: inode,kk,l,m,ni 
  complex(rp)  :: armonico,funq00,funq10,funq1m1,funq11,funq20,funq2m2
  complex(rp)  :: funq2m1,funq21,funq22
  character(5) :: wopos(2)

  L=0
  M=0
  call integracuad(l,M,funq00)

  L=1
  M=0
  call integracuad(l,M,funq10)

  L=1
  M=-1
  call integracuad(l,M,funq1m1)

  L=1
  M=1
  call integracuad(l,M,funq11)

  L=2
  M=0
  call integracuad(l,M,funq20)

  L=2
  M=-2
  call integracuad(l,M,funq2m2)

  L=2
  M=-1
  call integracuad(l,M,funq2m1)

  L=2
  M=1
  call integracuad(l,M,funq21)

  L=2
  M=2
  call integracuad(l,M,funq22)

  if( INOTSLAVE ) then
     print*,'Integrals (poibou)=',funq00,funq10,funq1m1,funq11
     print*,'Integrals (poibou)=',funq20,funq2m2,funq2m1,funq21
     print*,'Integrals (poibou)=',funq22
  end if

  if( INOTMASTER ) then

     do ni=1,npoin

        if(  kfl_fixno_qua(1,ni) == 1 .or. &
             kfl_fixno_qua(1,ni) == 4 .or. &
             kfl_fixno_qua(1,ni) == 5 ) then

           rmedio = sqrt(  coord(1,ni)**2 +  coord(2,ni)**2 +  coord(3,ni)**2)
           ante = 16.0_rp**atan(1.0_rp)/rmedio 

           vlim = ante *0.28209_rp*funq00

           bvessH_qua(ni,1)= vlim + ante/(3.0_rp*rmedio) * (                       &
                &                  funq10*armonico(coord(1,ni),coord(2,ni),coord(3,ni),1,0)   +    &
                &                  funq1m1*armonico(coord(1,ni),coord(2,ni),coord(3,ni),1,-1) +    &
                &                  funq11*armonico(coord(1,ni),coord(2,ni),coord(3,ni),1,1)   )   

           bvessH_qua(ni,1)= bvessH_qua(ni,1) + ante/(3.0_rp*rmedio)/(5.0_rp*rmedio) * (   &
                &                  funq20*armonico(coord(1,ni),coord(2,ni),coord(3,ni),2,0)   +    &
                &                  funq2m2*armonico(coord(1,ni),coord(2,ni),coord(3,ni),2,-2) +    &
                &                  funq2m1*armonico(coord(1,ni),coord(2,ni),coord(3,ni),2,-1) +    &
                &                  funq21*armonico(coord(1,ni),coord(2,ni),coord(3,ni),2,1) +      &
                &                  funq22*armonico(coord(1,ni),coord(2,ni),coord(3,ni),2,2)   )   

        endif
     enddo
  end if

  !wopos(1)='BVESS'
  !wopos(2)='SCALA'
  !call postpr(bvessH_qua(:,1),wopos,ittim,cutim)

end subroutine qua_poibou

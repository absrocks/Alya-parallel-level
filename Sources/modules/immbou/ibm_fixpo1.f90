subroutine ibm_fixpo1(dt,iimbo,tol,maxit)
  !-----------------------------------------------------------------------
  !****f* ibm_fixpo1/ibm_fixpo1
  ! NAME
  !    ibm_fixpo1
  ! DESCRIPTION
  !    This routines calculate the angular quantities for a rigid solid
  !    with fixed point iterative method. Full coupled algorithm.   
  ! USED BY
  !    eulerib
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_immbou
  implicit none
  real(rp),    intent(in) :: dt
  integer(ip), intent(in) :: iimbo,maxit
  real(rp),    intent(in) :: tol
  integer(ip)             :: idime,jdime,itera
  real(rp)                :: deter,numer,denom,auxir
  real(rp)                :: Taver(3)
  real(rp),    pointer    :: T(:,:),z(:,:),w(:,:),s(:,:),IT(:),R(:,:),q(:,:) ! Angular motion
  real(rp)                :: Io(3,3),In1(3,3),Iinv(3,3),wi(3),qtem(4),vtem(3),v2tem(3),Ttem(3,3),error

  T  => imbou(iimbo)%torqu
  IT => imbou(iimbo)%momin
  z  => imbou(iimbo)%accea
  w  => imbou(iimbo)%veloa
  s  => imbou(iimbo)%posia
  R  => imbou(iimbo)%rotac
  q  => imbou(iimbo)%quate

  Taver(3) = 0.0_rp
  if (kfl_nforc_ibm == 1) then
     do idime = 1,3
        Taver(idime) = T(idime,1)
     end do
  elseif (kfl_nforc_ibm == 2) then
     do idime = 1,3
        Taver(idime) = 0.5_rp*(T(idime,1) + T(idime,2))
     end do
  end if
  ! Initialitation
  Io(1,1)=IT(1); Io(1,2)=IT(4); Io(1,3)=IT(5)
  Io(2,1)=IT(4); Io(2,2)=IT(2); Io(2,3)=IT(6)
  Io(3,1)=IT(5); Io(3,2)=IT(6); Io(3,3)=IT(3)

  wi(1)=w(1,1); wi(2)=w(2,1) ;wi(3)=w(3,1)

  vtem(1)=0.0_rp; vtem(2)=0.0_rp; vtem(3)=0.0_rp

  error=1.0_rp
  itera=0_ip

  ! Iterate to find the angular quantities
  do while (error>tol .and. itera<maxit)  
  
     ! Compute the actual quaternion 
     ! 1. Compute the product the actual quaternion with the actual angular velocity
     qtem(1)= -wi(1)*q(2,1) - wi(2)*q(3,1) - wi(3)*q(4,1)
     qtem(2)=  q(1,1)*wi(1) + wi(2)*q(4,1) - wi(3)*q(3,1)
     qtem(3)=  q(1,1)*wi(2) + wi(3)*q(2,1) - wi(1)*q(4,1)
     qtem(4)=  q(1,1)*wi(3) + wi(1)*q(3,1) - wi(2)*q(2,1)
     ! 2. Add this result, multiplied by 0.5, with the quarterion from the previous time step
     numer=0.0_rp
     do idime=1,ndime+1
        q(idime,1)= q(idime,2) + 0.5_rp * dt * qtem(idime)              
        numer=numer + q(idime,1)*q(idime,1)
     end do     
     

     ! 3. Normalize actual quaternion
     if (numer > zeror) then
        do idime=1,ndime+1
           q(idime,1)= q(idime,1) / sqrt(numer)
        end do
     end if



     ! Obtain rotation matrix from actual quaternion
     R(1,1)= 1_rp - 2_rp*q(3,1)*q(3,1) - 2_rp*q(4,1)*q(4,1)
     R(2,2)= 1_rp - 2_rp*q(2,1)*q(2,1) - 2_rp*q(4,1)*q(4,1)
     R(3,3)= 1_rp - 2_rp*q(2,1)*q(2,1) - 2_rp*q(3,1)*q(3,1)
     R(1,2)=        2_rp*q(2,1)*q(3,1) - 2_rp*q(1,1)*q(4,1)
     R(2,1)=        2_rp*q(2,1)*q(3,1) + 2_rp*q(1,1)*q(4,1)
     R(1,3)=        2_rp*q(2,1)*q(4,1) + 2_rp*q(1,1)*q(3,1)
     R(3,1)=        2_rp*q(2,1)*q(4,1) - 2_rp*q(1,1)*q(3,1)
     R(2,3)=        2_rp*q(3,1)*q(4,1) - 2_rp*q(1,1)*q(2,1)
     R(3,2)=        2_rp*q(3,1)*q(4,1) + 2_rp*q(1,1)*q(2,1)

     ! Compute the actual inertia tensor
     call mbmab0(Ttem,R,Io,3,3,3)
     call mbmabt(In1,Ttem,R,3,3,3)

     ! Compute the actual inverse inertia tensor
     call invmtx(In1,Iinv,deter,3)

     ! Compute the actual angular acceleration
     call mbvab0(vtem,In1,wi,3,3)
     call vecpro(wi,vTem,v2tem,3)

     do idime = 1,ndime
        vtem(idime) = Taver(idime) - v2Tem(idime)                      
     end do

     ! Compute the actual angular quantities           
     do idime = 1,ndime 
        if ( kfl_ralei_ibm == 0 ) then          
           z(idime,1) = xrota_ibm(idime) * ( Iinv(idime,1)*vTem(1) + Iinv(idime,2)*vTem(2) + Iinv(idime,3)*vTem(3) )
        else
           
           if ( ittim <= nstra_ibm ) then 
              auxir = 1.0_rp
           else if ( ittim >= nenra_ibm ) then 
              auxir = 0.0_rp
           else
              auxir = dble ( nenra_ibm - ittim ) / dble ( nenra_ibm - nstra_ibm )
           end if

           z(idime,1) = xrota_ibm(idime) * ( Iinv(idime,1)*vTem(1) + Iinv(idime,2)*vTem(2) + Iinv(idime,3)*vTem(3) &
                - auxir * ralei_ibm * w(idime,2) )

        end if
        if( ittim == 1 ) z(idime,2) = z(idime,1)
        w(idime,1) = w(idime,2) + xrota_ibm(idime) * ( dt*z(idime,2) + gamma_ibm*dt*(z(idime,1) - z(idime,2))   )   
     end do

     ! Compute the angular velocity error
     numer=0.0_rp
     denom=0.0_rp
     do idime = 1,ndime              
        numer = numer + (w(idime,1) - wi(idime))*(w(idime,1) - wi(idime))
        denom = denom + w(idime,1)*w(idime,1)
     end do

     if (denom < zeror) then
        error = 0.0_rp
     else
        error = sqrt(numer)/sqrt(denom)
     end if

     ! Actualice the actual angular velocity iteration
     wi(1)=w(1,1); wi(2)=w(2,1) ;wi(3)=w(3,1) 

     vtem(1)=0_rp; vtem(2)=0_rp; vtem(3)=0_rp     
     itera=itera+1_ip

  end do
  !s(idime,1) = s(idime,2) + xrota_ibm(idime) * ( dt*w(idime,2) + 0.5_rp*dt*dt*z(idime,2) + beta_ibm*dt*dt*( z(idime,1) - z(idime,2) ) )


end subroutine ibm_fixpo1


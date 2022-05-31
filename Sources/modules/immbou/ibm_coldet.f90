subroutine ibm_coldet(ittib)
  !-----------------------------------------------------------------------
  !****f* ibm_coldet/ibm_coldet
  ! NAME
  !    ibm_coldet
  ! DESCRIPTION
  !    This routines detect a collision between particles A PRIORI
  ! OUTPUT
  !    twall_ibm(iwaib) % cotim =  Collision time particle-wall
  !    imbou(iimbo) % cotim = Collision time particle-particle
  !    dtime_ibm: minimum of all collision times
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only    :  ip,rp
  use def_master
  use def_domain
  use def_immbou
  use mod_kdtree
  implicit none

  integer(ip), intent(in) :: ittib
  integer(ip)             :: iimbo,jimbo,idime,iwaib
  integer(ip)             :: icond,ifind
  real(rp)                :: dt,tstep,minti
  real(rp)                :: dista,coori(3),coorj(3),norma(3)
  real(rp)                :: boboi(3,2),boboj(3,2)

  real(rp)                :: cotim(nimbo + nwaib)
  !
  ! Minimum time step (proportional to the global time step )
  ! 
  minti = 1.0_rp/1000.0_rp
  
  if( kfl_colli_ibm /= 0 ) then
     do iimbo = 1,nimbo 
        cotim(iimbo) = -1.0e10_rp
     end do
     do iwaib = 1,nwaib  
        cotim(nimbo + iwaib) = -1.0e10_rp
     end do
     !--------------------------------------------------------------------------
     !
     ! Estimate the time of collision between particles for an external time step
     !
     !---------------------------------------------------------------------------    
     dt = dtime_ibm
     do iimbo = 1,nimbo
        if ( imbou(iimbo) % cotim == -1.0e10_rp ) then
           call ibm_movbox(iimbo,2_ip,dtime_ibm,boboi)
           !
           ! IIMBO vs. JIMBO PARTICLES
           !
           do jimbo = iimbo+1,nimbo
              if ( imbou(jimbo) % cotim == -1.0e10_rp ) then
                 call ibm_movbox(jimbo,2_ip,dtime_ibm,boboj)
                 !
                 ! Determine if the bounding boxes are intersected
                 !
                 icond=0_ip
                 do idime=1,ndime
                    if ( boboi(idime,1) < boboj(idime,2) .and. boboi(idime,2) > boboj(idime,1) ) then
                       icond=icond+1_ip
                    end if
                 end do

                 if ( icond == ndime ) then


                    !
                    ! Determine the minimum distance betwwen the particles and its direction
                    !
                    call ibm_dpapar(iimbo,jimbo,coori,coorj,norma,dista)
                    call ibm_tpapar(iimbo,jimbo,coori,coorj,norma,dista,minti,tstep)       
                    cotim(jimbo) = min( abs(cotim(jimbo)) , cutim_ibm + tstep )                    
                    if ( tstep < dt  ) then
                       dt        = tstep
                    end if                    
                 end if
              end if
           end do
           !
           ! IIMBO PARTICLE vs. IWAIB WALL
           !
           do iwaib = 1,nwaib
              if ( twall_ibm(iwaib) % cotim == -1.0e10_rp ) then                 
                 icond=0_ip
                 ! Determine if the bounding boxes are intersected
                 !
                 do idime=1,ndime
                    if ( boboi(idime,1) < twall_ibm(iwaib) % bobox(idime,2) .and. boboi(idime,2) >  twall_ibm(iwaib) % bobox(idime,1) ) then
                       icond=icond+1_ip
                    end if
                 end do
                 if (icond==ndime) then  
                    call ibm_boxwal(boboi,iwaib,ifind)                    
                    if ( ifind == 1 ) then            
                       !
                       ! Determine the minimum distance betwwen the particles and its direction
                       !                                  
                       call ibm_dpawal(iimbo,iwaib,coori,coorj,norma,dista) 
                       call ibm_tpawal(iimbo,coori,norma,dista,minti,tstep) 
        
                       cotim(iimbo)         = min( abs(cotim(iimbo))         , cutim_ibm + tstep )
                       cotim(nimbo + iwaib) = min( abs(cotim(nimbo + iwaib)) , cutim_ibm + tstep )
                       if ( tstep < dt  ) then
                          dt = tstep
                       end if
                    end if
                 end if
              end if
           end do
        end if
        if( cotim(iimbo) == -1.0e10_rp) cotim(iimbo) = 1.0e10_rp
     end do
     dtime_ibm = dt

     !----------------------------------------------------
     !
     ! Reestimate the time of collision if it is neccesary  
     !
     !----------------------------------------------------   
     do iimbo = 1,nimbo  

        if ( imbou(iimbo) % cotim < cutim_ibm + dtime_ibm - zeror) then
           !
           ! Calcute the linear acceleration and the angular and linear displacement in the next time step. We suppose a linear acceleration
           ! We use the Taylor Series for the approximation
           !
           call ibm_movbox(iimbo,2_ip,dtime_ibm,boboi)
           !
           ! IIMBO vs. JIMBO PARTICLES
           !
           do jimbo = iimbo+1,nimbo
              if ( imbou(jimbo) % cotim < cutim_ibm + dtime_ibm - zeror) then                 
                 call ibm_movbox(jimbo,2_ip,dtime_ibm,boboj)
                 !
                 ! Determine if the bounding boxes are intersected
                 !
                 icond=0_ip
                  do idime=1,ndime
                    if ( boboi(idime,1) < boboj(idime,2) .and. boboi(idime,2) > boboj(idime,1) ) then
                       icond=icond+1_ip
                    end if
                 end do
                 if ( icond==ndime ) then
                    !
                    ! Determine the minimum distance betwwen the particles and its direction
                    !
                    call ibm_dpapar(iimbo,jimbo,coori,coorj,norma,dista)
                    call ibm_tpapar(iimbo,jimbo,coori,coorj,norma,dista,minti,tstep)                    
                    imbou(iimbo) % cotim = min( imbou(iimbo) % cotim, cutim_ibm + tstep )
                    imbou(jimbo) % cotim = min( imbou(jimbo) % cotim, cutim_ibm + tstep )        
                    if ( tstep < dt  ) then
                       dt = tstep
                    end if                               
                 end if
              end if
           end do
           !
           ! IIMBO PARTICLE vs. IWAIB WALL
           !
           do iwaib = 1,nwaib
              if ( twall_ibm(iwaib) % cotim < cutim_ibm + dtime_ibm - zeror) then
                 !
                 ! Determine if the bounding boxes are intersected
                 !
                 icond=0_ip
                 do idime=1,ndime
                    if ( boboi(idime,1) < twall_ibm(iwaib) % bobox(idime,2) .and. boboi(idime,2) >  twall_ibm(iwaib) % bobox(idime,1) ) then
                       icond=icond+1_ip
                    end if
                 end do
                 if (icond==ndime) then     
                    call ibm_boxwal(boboi,iwaib,ifind)                    
                    if ( ifind == 1 ) then
                       !
                       ! Determine the minimum distance betwwen the particles and its direction
                       !           
                       call ibm_dpawal(iimbo,iwaib,coori,coorj,norma,dista)
                       call ibm_tpawal(iimbo,coori,norma,dista,minti,tstep)                               
                       imbou(iimbo) % cotim     = min( imbou(iimbo) % cotim    , cutim_ibm + tstep )
                       twall_ibm(iwaib) % cotim = min( twall_ibm(iwaib) % cotim, cutim_ibm + tstep )                          
                       if ( tstep < dt  ) then
                          dt = tstep
                       end if
                    end if
                 end if
              end if
           end do
        end if
     end do
     dtime_ibm = dt
     !
     ! Parall: Compute minimum time step
     !
     call pararr('MIN',0_ip,1_ip,dtime_ibm)
     do iimbo = 1,nimbo
        cotim(iimbo) = imbou(iimbo) % cotim
     end do
     call pararr('MIN',0_ip,nimbo,cotim)
     do iimbo = 1,nimbo
        imbou(iimbo) % cotim = cotim(iimbo) 
     end do
  end if
end subroutine ibm_coldet


subroutine ibm_tpapar(iimbo,jimbo,coori,coorj,direc,dista,minti,tstep)
  !-----------------------------------------------------------------------
  ! NAME
  !    ibm_contib
  ! DESCRIPTION
  !    This routines calculate the approximate time step of collision beetween two particles
  ! USED BY
  !    
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  imbou,cutim,dtime,zeror
  use def_domain, only     :  ndime,nnode
  use def_immbou, only     :  dtime_ibm,cutim_ibm
  implicit none
  integer(ip), intent(in)  :: iimbo,jimbo
  real(rp),    intent(in)  :: coori(3),coorj(3),direc(3),dista,minti
  real(rp),    intent(out) :: tstep
  integer(ip)              :: idime,itera,salir
  real(rp)                 :: toler,fdt,dfdt,newdt,olddt,dt
  real(rp)                 :: root1,root2,fact1,fact2
  real(rp)                 :: rmaxi,rmaxj,ri(ndime),rj(ndime)
  real(rp), pointer        :: ai(:,:),vi(:,:),xi(:,:),zi(:,:),wi(:,:),si(:,:)
  real(rp), pointer        :: aj(:,:),vj(:,:),xj(:,:),zj(:,:),wj(:,:),sj(:,:)
  !
  ! Determine the distance from the closest point to the mass center for IMMBO 
  !
  xi    => imbou(iimbo) % posil
  vi    => imbou(iimbo) % velol
  ai    => imbou(iimbo) % accel
  si    => imbou(iimbo) % posia
  wi    => imbou(iimbo) % veloa
  zi    => imbou(iimbo) % accea
  rmaxi = 0.0_rp
  do idime = 1,ndime
     rmaxi  = rmaxi + direc(idime)*(coori(idime) - imbou(iimbo)%posil(idime,2))
  end do
  rmaxi = abs(rmaxi)
  !
  ! Determine the distance from the closest point to the mass center for JMMBO 
  !  
  xj    => imbou(jimbo) % posil
  vj    => imbou(jimbo) % velol
  aj    => imbou(jimbo) % accel
  sj    => imbou(jimbo) % posia
  wj    => imbou(jimbo) % veloa
  zj    => imbou(jimbo) % accea
  rmaxj = 0.0_rp
  do idime = 1,ndime
     rmaxj  = rmaxj + direc(idime)*(coorj(idime) - imbou(jimbo)%posil(idime,2))
  end do
  rmaxj = abs(rmaxj)

  dt    = 0.0
  salir = 0_ip
  !
  ! Incerement
  !      
  do while (salir==0)
     !
     ! Newthon-Raphson
     !      
     toler = 1.0_rp
     itera = 0_ip
     newdt = dt
     do while ( abs(toler) > 1.0e-10 .and. itera < 100_ip )
        itera = itera + 1_ip
        olddt = newdt
        !
        ! Determine the value of the function
        !     
        fdt = 0.0_rp
        root1 = 0.0_rp
        root2 = 0.0_rp
        do idime = 1,ndime
           fdt = fdt &
                - (vi(idime,2)*olddt + 0.5_rp*ai(idime,3)*olddt*olddt + &
                (1.0_rp/6.0_rp)*(1.0_rp/dtime)*(ai(idime,1) - ai(idime,2))*olddt*olddt*olddt)*direc(idime) &
                + (vj(idime,2)*olddt + 0.5_rp*aj(idime,3)*olddt*olddt + &
                (1.0_rp/6.0_rp)*(1.0_rp/dtime)*(aj(idime,1) - aj(idime,2))*olddt*olddt*olddt)*direc(idime)
           root1 = root1 &
                + (wi(idime,2)*olddt + 0.5_rp*zi(idime,3)*olddt*olddt + &
                (1.0_rp/6.0_rp)*(1.0_rp/dtime)*(zi(idime,1) - zi(idime,2))*olddt*olddt*olddt)**2.0_rp
           root2 = root2 &
                + (wj(idime,2)*olddt + 0.5_rp*zj(idime,3)*olddt*olddt + &
                (1.0_rp/6.0_rp)*(1.0_rp/dtime)*(zj(idime,1) - zj(idime,2))*olddt*olddt*olddt)**2.0_rp
        end do
        fdt = fdt + min( sqrt(root1)*imbou(iimbo) %maxdi,imbou(iimbo) % maxdi - rmaxi)
        fdt = fdt + min( sqrt(root2)*imbou(jimbo) %maxdi,imbou(jimbo) % maxdi - rmaxj) &
             - dista + min(imbou(iimbo) %maxdi,imbou(jimbo) %maxdi)*2.5e-3_rp
        !
        ! Determine the value of the derivative
        !     
        dfdt = 0.0_rp
        fact1 = 0.0_rp
        fact2 = 0.0_rp
        do idime = 1,ndime
           dfdt = dfdt &
                - (vi(idime,2) + ai(idime,3)*olddt + 0.5_rp*(1.0_rp/dtime)*(ai(idime,1) - ai(idime,2))*olddt*olddt)*direc(idime) &
                + (vj(idime,2) + aj(idime,3)*olddt + 0.5_rp*(1.0_rp/dtime)*(aj(idime,1) - aj(idime,2))*olddt*olddt)*direc(idime)        
           fact1 = fact1 &
                + 2.0_rp*(wi(idime,2)*olddt + 0.5_rp*zi(idime,3)*olddt*olddt + &
                (1.0_rp/6.0_rp)*(1.0_rp/dtime)*(zi(idime,1) - zi(idime,2))*olddt*olddt*olddt) &
                *(wi(idime,2) + zi(idime,3)*olddt + 0.5_rp*(1.0_rp/dtime)*(zi(idime,1) - zi(idime,2))*olddt*olddt)
           fact2 = fact2 &
                + 2.0_rp*(wj(idime,2)*olddt + 0.5_rp*zj(idime,3)*olddt*olddt + &
                (1.0_rp/6.0_rp)*(1.0_rp/dtime)*(zj(idime,1) - zj(idime,2))*olddt*olddt*olddt) &
                *(wj(idime,2) + zj(idime,3)*olddt + 0.5_rp*(1.0_rp/dtime)*(zj(idime,1) - zj(idime,2))*olddt*olddt)
        end do
        if ( abs(root1) > zeror ) then
           if ( sqrt(root1)*imbou(iimbo) %maxdi < imbou(iimbo) %maxdi - rmaxi) then
              dfdt = dfdt + 0.5_rp*(1.0_rp/sqrt(root1))*fact1*imbou(iimbo) %maxdi
           end if
        end if
        if ( abs(root2) > zeror ) then
           if ( sqrt(root2)*imbou(jimbo) %maxdi < imbou(jimbo) %maxdi - rmaxj) then
              dfdt = dfdt + 0.5_rp*(1.0_rp/sqrt(root2))*fact2*imbou(jimbo) %maxdi
           end if
        end if
        !
        ! Determine the value of the actual iteration 
        !          
        if ( abs(dfdt) > zeror ) then
           newdt =  olddt - (fdt/dfdt)
        else
           newdt = cutim - cutim_ibm
           itera = 100_ip
        end if
        !
        ! Normalized error
        !
        if ( (newdt*newdt) > zeror ) then
           toler = ((newdt - olddt)*(newdt - olddt)) / (newdt*newdt)         
        else
           toler = ((newdt - olddt)*(newdt - olddt)) / zeror
        end if
     end do

     if (newdt > dt .or. dt >= dtime_ibm) then
        salir = 1_ip 
     else
        dt = dt + dtime*100.0_rp*minti
     end if
  end do


  if ( newdt < 0.0_rp ) then
     newdt = 1.0e10_rp
  else if ( newdt < dtime*minti ) then
     newdt = dtime*minti
  end if
  tstep = newdt
end subroutine ibm_tpapar




subroutine ibm_tpawal(iimbo,coori,direc,dista,minti,tstep)
  !-----------------------------------------------------------------------
  ! NAME
  !    ibm_contib
  ! DESCRIPTION
  !    This routines calculate the approximate time step of collision beetween a particle
  !    and a wall
  ! USED BY
  !    
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  imbou,cutim,dtime,zeror
  use def_domain, only     :  ndime,nnode
  use def_immbou, only     :  dtime_ibm,cutim_ibm
  implicit none
  integer(ip), intent(in)  :: iimbo
  real(rp),    intent(in)  :: coori(3),direc(3),dista,minti
  real(rp),    intent(out) :: tstep
  integer(ip)              :: idime,itera,salir
  real(rp)                 :: toler,fdt,dfdt,newdt,olddt,dt
  real(rp)                 :: root1,fact1
  real(rp)                 :: rmaxi,ri(ndime)
  real(rp), pointer        :: ai(:,:),vi(:,:),xi(:,:),zi(:,:),wi(:,:),si(:,:)
  !
  ! Determine the distance from the closest point to the mass center for IMMBO 
  !
  xi    => imbou(iimbo) % posil
  vi    => imbou(iimbo) % velol
  ai    => imbou(iimbo) % accel
  si    => imbou(iimbo) % posia
  wi    => imbou(iimbo) % veloa
  zi    => imbou(iimbo) % accea
  rmaxi = 0.0_rp
  do idime = 1,ndime
     rmaxi  = rmaxi + direc(idime)*(coori(idime) - imbou(iimbo)%posil(idime,2))
  end do
  rmaxi = abs(rmaxi)  

  dt    = 0.0
  salir = 0_ip
  !
  ! Incerement
  !      
  do while (salir==0)
     !
     ! Newthon-Raphson
     !      
     toler = 1.0_rp
     itera = 0_ip
     newdt = dt
     do while ( abs(toler) > 1.0e-10 .and. itera < 100_ip)
        itera = itera + 1_ip
        olddt = newdt
        !
        ! Determine the value of the function
        !
        fdt = 0.0_rp
        root1 = 0.0_rp
        do idime = 1,ndime
           fdt = fdt &
                - (vi(idime,2)*olddt + 0.5_rp*ai(idime,3)*olddt*olddt + &
                (1.0_rp/6.0_rp)*(1.0_rp/dtime)*(ai(idime,1) - ai(idime,2))*olddt*olddt*olddt)*direc(idime)
           root1 = root1 &
                + (wi(idime,2)*olddt + 0.5_rp*zi(idime,3)*olddt*olddt + &
                (1.0_rp/6.0_rp)*(1.0_rp/dtime)*(zi(idime,1) - zi(idime,2))*olddt*olddt*olddt)**2.0_rp
        end do
        fdt = fdt + min( sqrt(root1)*imbou(iimbo) %maxdi,imbou(iimbo) % maxdi - rmaxi) &
             - dista + imbou(iimbo) %maxdi*2.5e-3_rp
        !
        ! Determine the value of the derivative
        !     
        dfdt = 0.0_rp
        fact1 = 0.0_rp
        do idime = 1,ndime
           dfdt = dfdt &
                - (vi(idime,2) + ai(idime,3)*olddt + 0.5_rp*(1.0_rp/dtime)*(ai(idime,1) - ai(idime,2))*olddt*olddt)*direc(idime)        
           fact1 = fact1 &
                + 2.0_rp*(wi(idime,2)*olddt + 0.5_rp*zi(idime,3)*olddt*olddt + &
                (1.0_rp/6.0_rp)*(1.0_rp/dtime)*(zi(idime,1) - zi(idime,2))*olddt*olddt*olddt) &
                *(wi(idime,2) + zi(idime,3)*olddt + 0.5_rp*(1.0_rp/dtime)*(zi(idime,1) - zi(idime,2))*olddt*olddt)
        end do
        if ( abs(root1) > zeror ) then
           if ( sqrt(root1)*imbou(iimbo) %maxdi < imbou(iimbo) %maxdi - rmaxi) then
              dfdt = dfdt + 0.5_rp*(1.0_rp/sqrt(root1))*fact1*imbou(iimbo) %maxdi
           end if
        end if
        !
        ! Determine the value of the actual iteration 
        !          
        if ( abs(dfdt) > zeror ) then
           newdt =  olddt - (fdt/dfdt)
        else
           newdt = cutim - cutim_ibm
           itera = 100_ip
        end if
        !
        ! Normalized error
        !
        toler = ((newdt - olddt)*(newdt - olddt)) / (newdt*newdt)
     end do

     if (newdt > dt .or. dt >= dtime_ibm) then
        salir = 1_ip 
     else
        dt = dt + dtime*100.0_rp*minti
     end if
  end do

  if ( newdt < 0.0_rp ) then
     newdt = 1.0e10_rp
  else if ( newdt < dtime*minti ) then
     newdt = dtime*minti
  end if
  tstep = newdt
end subroutine ibm_tpawal


!!$subroutine ibm_tpapar(iimbo,jimbo,coori,coorj,direc,dista,minti,tstep)
!!$  !-----------------------------------------------------------------------
!!$  ! NAME
!!$  !    ibm_contib
!!$  ! DESCRIPTION
!!$  !    This routines calculate the approximate time step of collision beetween two particles
!!$  ! USED BY
!!$  !    
!!$  !--------------------------------- --------------------------------------
!!$  use def_kintyp, only     :  ip,rp
!!$  use def_master, only     :  imbou,dtime,zeror
!!$  use def_domain, only     :  ndime
!!$  use def_immbou, only     :  dtime_ibm
!!$  implicit none
!!$  integer(ip), intent(in)  :: iimbo,jimbo  
!!$  real(rp),    intent(in)  :: coori(ndime),coorj(ndime)
!!$  real(rp),    intent(in)  :: direc(ndime)
!!$  real(rp),    intent(in)  :: dista
!!$  real(rp),    intent(in)  :: minti
!!$  real(rp),    intent(out) :: tstep
!!$
!!$  integer(ip)              :: idime,iroot,zmaxi,zmaxj
!!$  real(rp)                 :: magi1,magi2,magj1,magj2
!!$  real(rp)                 :: discr,roots(3),dt,tenor(ndime),tedis
!!$  real(rp)                 :: A,B,C,D ! coefficients of the cubic equation 
!!$  real(rp)                 :: F,G,H,I,J,K,L,M,N,P,R,S,T,U ! constants use to solve cubic equation
!!$  
!!$  ! Actual time step
!!$  dt    =  dtime_ibm
!!$  tstep = 1.0e10_rp
!!$
!!$  magi1 = 0.0_rp; magi2 = 0.0_rp
!!$  magj1 = 0.0_rp; magj2 = 0.0_rp
!!$  do idime = 1,ndime
!!$     magi1 = magi1 + imbou(iimbo)%accel(idime,1)*imbou(iimbo)%accel(idime,1)
!!$     magi2 = magi2 + imbou(iimbo)%accel(idime,3)*imbou(iimbo)%accel(idime,3)
!!$     magj1 = magj1 + imbou(jimbo)%accel(idime,1)*imbou(jimbo)%accel(idime,1)
!!$     magj2 = magj2 + imbou(jimbo)%accel(idime,3)*imbou(jimbo)%accel(idime,3)
!!$  end do
!!$  if (magi1 > magi2) then
!!$     zmaxi = sqrt(magi1)
!!$  else
!!$     zmaxi = sqrt(magi2)
!!$  end if
!!$  if (magj1 > magj2) then
!!$     zmaxj = sqrt(magj1)
!!$  else
!!$     zmaxj = sqrt(magj2)
!!$  end if
!!$ 
!!$  ! Obtain the coefficients of the cubic equation. Use Newmark equation for displacement and the fact that the acceleration is considered linear 
!!$  A = 0.0_rp; B = 0.0_rp; C = 0.0_rp; D = 0.0_rp
!!$  do idime = 1,ndime 
!!$     A = A + ( 1.0_rp/(6.0_rp*dt)) * ( imbou(iimbo)%accel(idime,1) - imbou(iimbo)%accel(idime,2) & 
!!$                                      -imbou(jimbo)%accel(idime,1) + imbou(jimbo)%accel(idime,2) ) * direc(idime)     
!!$     B = B + ( 1.0_rp/2.0_rp     ) * ( imbou(iimbo)%accel(idime,3) - imbou(jimbo)%accel(idime,3) ) * direc(idime)
!!$     C = C +                         ( imbou(iimbo)%velol(idime,2) - imbou(jimbo)%velol(idime,2) ) * direc(idime)
!!$  end do
!!$  B = B + zmaxi*imbou(iimbo)%maxdi + zmaxj*imbou(jimbo)%maxdi
!!$  D = -dista + dista*1.0e-2_rp
!!$
!!$  roots(1) = 1.0e10_rp
!!$  roots(2) = 1.0e10_rp
!!$  roots(3) = 1.0e10_rp
!!$ 
!!$  ! Solve the cubic equation
!!$  if ( abs(A) > zeror ) then
!!$     F = ( 3.0_rp*C/A - (B*B)/(A*A) )/3.0_rp
!!$     G = ( 2.0_rp*(B*B*B)/(A*A*A) - 9.0_rp*B*C/(A*A) + 27.0_rp*D/A  )/27.0_rp
!!$     H = ( (G*G)/4.0_rp + (F*F*F)/27.0_rp )
!!$     if ( abs(F)>zeror .and. abs(G)>zeror .and. H<=-zeror ) then
!!$        I = sqrt( (G*G)/4.0_rp - H )
!!$        J = 1.0_rp/(I*I*I)
!!$        K = acos( - G/(2.0_rp*I) )
!!$        L = J*-1.0_rp
!!$        M = cos(K/3.0_rp)
!!$        N = sqrt(3.0_rp)*sin(K/3.0_rp)
!!$        P = B/(3.0_rp*A)*-1.0_rp        
!!$        roots(1) = 2.0_rp*J*M - B/(3.0_rp*A)        
!!$        roots(2) = L*(M+N) + P
!!$        roots(3) = L*(M-N) + P
!!$     elseif ( H > zeror  ) then
!!$        R = - (G/2.0_rp) + sqrt(H)
!!$        S = 1.0_rp/(R*R*R)
!!$        T = -G/2.0_rp - sqrt(H)
!!$        U = 1.0_rp/(T*T*T)
!!$        roots(1) = (S+U) - B/(3*A)
!!$     elseif ( abs(F)<zeror .and. abs(G)<zeror .and. abs(H)<zeror ) then
!!$        roots(1) = (D/A)**(1.0_rp/3.0_rp)*-1.0_rp
!!$        roots(2) = (D/A)**(1.0_rp/3.0_rp)*-1.0_rp
!!$        roots(3) = (D/A)**(1.0_rp/3.0_rp)*-1.0_rp
!!$     end if
!!$ ! Solve a quadratic equation        
!!$  elseif ( abs(B) > zeror ) then
!!$     discr = C*C - 4.0_rp*B*D
!!$     if ( discr > zeror ) then
!!$        roots(1) = ( -C + sqrt(C*C - 4.0_rp*B*D) ) / ( 2.0_rp*B )
!!$        roots(2) = ( -C - sqrt(C*C - 4.0_rp*B*D) ) / ( 2.0_rp*B )
!!$     end if
!!$ ! Solve a linear equation        
!!$  else             
!!$     if (abs(C) > zeror) then
!!$        roots(1) = - D/C
!!$     end if
!!$  end if
!!$
!!$  print *,"cubes: ", iimbo,jimbo
!!$  print *,"roots: ", roots(:)
!!$  
!!$  tstep = 1.0e10_rp
!!$  do iroot = 1,3
!!$     if (roots(iroot) >= 0.0_rp) then
!!$        tstep = min(roots(iroot), tstep)
!!$     end if
!!$  end do  
!!$  
!!$end subroutine ibm_tpapar
!!$
!!$
!!$subroutine ibm_tpawal(iimbo,coori,direc,dista,minti,tstep)
!!$  !-----------------------------------------------------------------------
!!$  ! NAME
!!$  !    ibm_contib
!!$  ! DESCRIPTION
!!$  !    This routines calculate the approximate time step of collision beetween two particles
!!$  ! USED BY
!!$  !    
!!$  !--------------------------------- --------------------------------------
!!$  use def_kintyp, only     :  ip,rp
!!$  use def_master, only     :  imbou,dtime,zeror
!!$  use def_domain, only     :  ndime
!!$  use def_immbou, only     :  dtime_ibm
!!$  implicit none
!!$  integer(ip), intent(in)  :: iimbo
!!$  real(rp),    intent(in)  :: coori(ndime)
!!$  real(rp),    intent(in)  :: direc(ndime)
!!$  real(rp),    intent(in)  :: dista
!!$  real(rp),    intent(in)  :: minti
!!$  real(rp),    intent(out) :: tstep
!!$
!!$  integer(ip)              :: idime,iroot,zmaxi
!!$  real(rp)                 :: magi1,magi2
!!$  real(rp)                 :: discr,roots(3),dt
!!$  real(rp)                 :: A,B,C,D ! coefficients of the cubic equation 
!!$  real(rp)                 :: F,G,H,I,J,K,L,M,N,P,R,S,T,U ! constants use to solve cubic equation
!!$  
!!$  ! Actual time step
!!$  dt    =  dtime_ibm
!!$  tstep = 1.0e10_rp
!!$
!!$  magi1 = 0.0_rp; magi2 = 0.0_rp
!!$  do idime = 1,ndime
!!$     magi1 = magi1 + imbou(iimbo)%accel(idime,1)*imbou(iimbo)%accel(idime,1)
!!$     magi2 = magi2 + imbou(iimbo)%accel(idime,3)*imbou(iimbo)%accel(idime,3)
!!$  end do
!!$  if (magi1 > magi2) then
!!$     zmaxi = sqrt(magi1)
!!$  else
!!$     zmaxi = sqrt(magi2)
!!$  end if
!!$
!!$  ! Obtain the coefficients of the cubic equation. Use Newmark equation for displacement and the fact that the acceleration is considered linear 
!!$  A = 0.0_rp; B = 0.0_rp; C = 0.0_rp; D = 0.0_rp
!!$  do idime = 1,ndime 
!!$     A = A + ( 1.0_rp/(6.0_rp*dt)) * ( imbou(iimbo)%accel(idime,1) - imbou(iimbo)%accel(idime,2) ) * direc(idime)     
!!$     B = B + ( 1.0_rp/2.0_rp     ) * ( imbou(iimbo)%accel(idime,3)                               ) * direc(idime)
!!$     C = C +                         ( imbou(iimbo)%velol(idime,2)                               ) * direc(idime)
!!$  end do
!!$  B = B + zmaxi*imbou(iimbo)%maxdi
!!$  D = -dista + dista*1.0e-2_rp
!!$
!!$  roots(1) = 1.0e10_rp
!!$  roots(2) = 1.0e10_rp
!!$  roots(3) = 1.0e10_rp
!!$ 
!!$  ! Solve the cubic equation
!!$  if ( abs(A) > zeror ) then
!!$     F = ( 3.0_rp*C/A - (B*B)/(A*A) )/3.0_rp
!!$     G = ( 2.0_rp*(B*B*B)/(A*A*A) - 9.0_rp*B*C/(A*A) + 27.0_rp*D/A  )/27.0_rp
!!$     H = ( (G*G)/4.0_rp + (F*F*F)/27.0_rp )
!!$     if ( abs(F)>zeror .and. abs(G)>zeror .and. H<=-zeror ) then
!!$        I = sqrt( (G*G)/4.0_rp - H )
!!$        J = 1.0_rp/(I*I*I)
!!$        K = acos( - G/(2.0_rp*I) )
!!$        L = J*-1.0_rp
!!$        M = cos(K/3.0_rp)
!!$        N = sqrt(3.0_rp)*sin(K/3.0_rp)
!!$        P = B/(3.0_rp*A)*-1.0_rp        
!!$        roots(1) = 2.0_rp*J*M - B/(3.0_rp*A)        
!!$        roots(2) = L*(M+N) + P
!!$        roots(3) = L*(M-N) + P
!!$     elseif ( H > zeror  ) then
!!$        R = - (G/2.0_rp) + sqrt(H)
!!$        S = 1.0_rp/(R*R*R)
!!$        T = -G/2.0_rp - sqrt(H)
!!$        U = 1.0_rp/(T*T*T)
!!$        roots(1) = (S+U) - B/(3*A)
!!$     elseif ( abs(F)<zeror .and. abs(G)<zeror .and. abs(H)<zeror ) then
!!$        roots(1) = (D/A)**(1.0_rp/3.0_rp)*-1.0_rp
!!$        roots(2) = (D/A)**(1.0_rp/3.0_rp)*-1.0_rp
!!$        roots(3) = (D/A)**(1.0_rp/3.0_rp)*-1.0_rp
!!$     end if
!!$ ! Solve a quadratic equation        
!!$  elseif ( abs(B) > zeror ) then
!!$     discr = C*C - 4.0_rp*B*D
!!$     if ( discr > zeror ) then
!!$        roots(1) = ( -C + sqrt(C*C - 4.0_rp*B*D) ) / ( 2.0_rp*B )
!!$        roots(2) = ( -C - sqrt(C*C - 4.0_rp*B*D) ) / ( 2.0_rp*B )
!!$     end if
!!$ ! Solve a linear equation        
!!$  else             
!!$     if (abs(C) > zeror) then
!!$        roots(1) = - D/C
!!$     end if
!!$  end if
!!$
!!$
!!$  print *,"walls: ", iimbo
!!$  print *,"roots: ", roots(:)
!!$
!!$  tstep = 1.0e10_rp
!!$  do iroot = 1,3
!!$     if (roots(iroot) >= 0.0_rp) then
!!$        tstep = min(roots(iroot), tstep)
!!$     end if
!!$  end do  
!!$end subroutine ibm_tpawal

subroutine ibm_newmak(dt)
  !-----------------------------------------------------------------------
  !****f* ibm_newmak/ibm_newmak
  ! NAME
  !    ibm_newmak
  ! DESCRIPTION
  !    This routines detect a collision between particles A PRIORI
  !                             gamma   beta
  !     Fox Goodwin             1 / 2   1 / 12    conditionnaly stable
  !     Linear acceleration     1 / 2   1 / 6     conditionnaly stable
  !     Average acceleration    1 / 2   1 / 4     inconditionnaly stable
  !     Nosotros                0.75    0.390625  diffusive
  !     Nosotros                1.0     0.5625    super diffusive
  !     External                1 / 2   0         explicit 
  !
  !     Relation between beta and gamma: beta = 0.25*(gamma+1/2)^2
  !
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master
  use def_domain
  use def_immbou
  use mod_kdtree
  implicit none

  real(rp),  intent(in) :: dt
  integer(ip)           :: iimbo,idime,ipoib
  integer(ip)           :: itera,maxit
  real(rp)              :: error,numer,denom,tol
  real(rp)              :: qtem(4),xcoor(3),rotco(3)
  real(rp), pointer     :: ai(:,:),vi(:,:),xi(:,:),zi(:,:),wi(:,:),si(:,:),qi(:,:),Ri(:,:),Rt(:,:),xt(:)

  !----------------------------------------------------------------------
  !
  ! Update the particle data
  ! 
  !----------------------------------------------------------------------
  do iimbo = 1,nimbo  
     xi => imbou(iimbo) % posil
     vi => imbou(iimbo) % velol
     ai => imbou(iimbo) % accel
     si => imbou(iimbo) % posia
     wi => imbou(iimbo) % veloa
     zi => imbou(iimbo) % accea
     qi => imbou(iimbo) % quate
     Ri => imbou(iimbo) % rotac
     do idime = 1,3        

        xi(idime,1) = xi(idime,2) + xline_ibm(idime) * (dt*vi(idime,2) + 0.5_rp*dt*dt*ai(idime,3) + &
             beta_ibm*dt*dt*( ai(idime,1) - ai(idime,2)))
        vi(idime,1) = vi(idime,2) + xline_ibm(idime) * (dt*ai(idime,3) + gamma_ibm*dt*(ai(idime,1) - ai(idime,2)))
        ai(idime,3) = ai(idime,2) + xline_ibm(idime) * ((1.0_rp/dtime) * (ai(idime,1) - ai(idime,2))*(cutim_ibm - cutim + dtime))

        si(idime,1) = si(idime,2) + xrota_ibm(idime) * (dt*wi(idime,2) + 0.5_rp*dt*dt*zi(idime,3) + &
             beta_ibm*dt*dt*( zi(idime,1) - zi(idime,2)))
        wi(idime,1) = wi(idime,2) + xrota_ibm(idime) * (dt*zi(idime,3) + gamma_ibm*dt*(zi(idime,1) - zi(idime,2)))
        zi(idime,3) = zi(idime,2) + xrota_ibm(idime) * &
             ((1.0_rp/dtime) * (imbou(iimbo) % accel(idime,1) - imbou(iimbo) % accel(idime,2))*(cutim_ibm - cutim + dtime))

     end do


     if (ndime == 2) then
        Ri(1,1) =  COS(si(3,1))
        Ri(1,2) = -SIN(si(3,1))
        Ri(1,3) =  0.0_rp
        Ri(2,1) =  SIN(si(3,1))
        Ri(2,2) =  COS(si(3,1))
        Ri(2,3) =  0.0_rp
        Ri(3,1) =  0.0_rp
        Ri(3,2) =  0.0_rp
        Ri(3,3) =  1.0_rp
     else
        !
        ! 2. Add this result, multiplied by 0.5, with the quarterion from the previous time step
        !
        numer=0.0_rp

        qtem(1) = qi(1,2) + 0.5_rp * dt *( -wi(1,1)*qi(2,1) - wi(2,1)*qi(3,1) - wi(3,1)*qi(4,1))
        qtem(2) = qi(2,2) + 0.5_rp * dt *(  qi(1,1)*wi(1,1) + wi(2,1)*qi(4,1) - wi(3,1)*qi(3,1))
        qtem(3) = qi(3,2) + 0.5_rp * dt *(  qi(1,1)*wi(2,1) + wi(3,1)*qi(2,1) - wi(1,1)*qi(4,1))
        qtem(4) = qi(4,2) + 0.5_rp * dt *(  qi(1,1)*wi(3,1) + wi(1,1)*qi(3,1) - wi(2,1)*qi(2,1))

        numer   = qtem(1)*qtem(1) + qtem(2)*qtem(2) + qtem(3)*qtem(3) + qtem(4)*qtem(4)
        !
        ! 3. Normalize actual quaternion
        !
        if (numer/=0.0_rp) then
           numer = sqrt(numer)
           qtem(1) = qtem(1) / numer
           qtem(2) = qtem(2) / numer
           qtem(3) = qtem(3) / numer
           qtem(4) = qtem(4) / numer
        end if
        do idime=1,ndime+1
           qi(idime,1)= qtem(idime)
        end do
        Ri(1,1)= 1_rp - 2_rp*qi(3,1)*qi(3,1) - 2_rp*qi(4,1)*qi(4,1)
        Ri(2,2)= 1_rp - 2_rp*qi(2,1)*qi(2,1) - 2_rp*qi(4,1)*qi(4,1)
        Ri(3,3)= 1_rp - 2_rp*qi(2,1)*qi(2,1) - 2_rp*qi(3,1)*qi(3,1)
        Ri(1,2)=        2_rp*qi(2,1)*qi(3,1) - 2_rp*qi(1,1)*qi(4,1)
        Ri(2,1)=        2_rp*qi(2,1)*qi(3,1) + 2_rp*qi(1,1)*qi(4,1)
        Ri(1,3)=        2_rp*qi(2,1)*qi(4,1) + 2_rp*qi(1,1)*qi(3,1)
        Ri(3,1)=        2_rp*qi(2,1)*qi(4,1) - 2_rp*qi(1,1)*qi(3,1)
        Ri(2,3)=        2_rp*qi(3,1)*qi(4,1) - 2_rp*qi(1,1)*qi(2,1)
        Ri(3,2)=        2_rp*qi(3,1)*qi(4,1) + 2_rp*qi(1,1)*qi(2,1)
     end if
     !
     ! Update coordinates
     !          
     do ipoib = 1,imbou(iimbo) % npoib
        do idime = 1,ndime
           xcoor(idime) = imbou(iimbo) % cooin(idime,ipoib)
        end do
        call mbvab0(rotco,Ri,xcoor,3_ip,3_ip)      
        do idime = 1,ndime
           imbou(iimbo)%cooi2(idime,ipoib) = imbou(iimbo)%cooib(idime,ipoib)
           imbou(iimbo)%cooib(idime,ipoib) = rotco(idime) + xi(idime,1) 
        end do
     end do
     !
     ! Update the bounding box for each particle and for each face in a particle.
     ! Also, update the tree structures used by find points inside the particle
     !           
     if( IMASTER .or. imbou(iimbo) % kfl_typeb >= 1 ) then
        continue
     else
        call kdtree(&
             0_ip,mnoib,imbou(iimbo) % npoib,imbou(iimbo) % nboib,&
             imbou(iimbo) % cooib,imbou(iimbo) %lnoib,imbou(iimbo) %ltyib,&
             imbou(iimbo) % fabox,imbou(iimbo) % bobox,imbou(iimbo) % sabox,&
             imbou(iimbo) % blink,imbou(iimbo) % struc,imbou(iimbo) % ldist,&
             imbou(iimbo) % lnele)

     end if
  end do

end subroutine ibm_newmak

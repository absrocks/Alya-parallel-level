subroutine nsa_bouwal(&
     pnodb,iboun,lboel,gbsha,bovel,bovfi,tract,gbvis,&
     gbden,baloc,velfr,rough)

  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_bouwal
  ! NAME 
  !    nsa_bouwal
  ! DESCRIPTION
  !    This routine computes the surface traction for the NS equations at
  !    a given integration point of a boundary IBOUN received by argument
  !    due to the use of a turbulent wall law. The algorithm is:
  !    - Compute the tangential velocity u at the integration point x.
  !      In fact, u is not necessarily tangential at x. For example:
  !      -> u    -> u      
  !      o---x---o---x---o
  !                      |\ u
  !                      | 
  !    - Given the distance to the wall y, compute U*
  !    - Compute y+=yU*/nu
  !      if(y+>5) then
  !        t=-rho*U*^2*(u_tan-u_fix_tan)/|u_tan-u_fix_tan|
  !      else
  !        u+=y+ => U*^2=u*nu/y so that
  !        t=-mu*u/y
  !      end if
  ! USES
  !    vecnor
  !    frivel
  ! USED BY
  !    nsa_bouope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_kermod, only     :  kfl_delta,kfl_ustar,delta_dom
  use def_domain, only     :  ndime,ywalb
  use def_nastal

  implicit none
  integer(ip), intent(in)  :: pnodb,iboun
  integer(ip), intent(in)  :: lboel(pnodb)
  real(rp),    intent(in)  :: gbvis,gbden
  real(rp),    intent(in)  :: bovel(ndime,pnodb),gbsha(pnodb)
  real(rp),    intent(in)  :: bovfi(ndime,pnodb)
  real(rp),    intent(in)  :: baloc(ndime,ndime),rough
  real(rp),    intent(out) :: velfr
  real(rp),    intent(out) :: tract(ndime)
  integer(ip)              :: ldime,idime,inodb,idofn
  integer(ip)              :: jevab,ievab,jdime,jnodb
  real(rp)                 :: veloc(3),fact1,fact2,velfi(3)
  real(rp)                 :: vikin,yplus                      ! nu, U*, y+, y
  real(rp)                 :: tveno                            ! |u_tan-u_fix_tan|
  real(rp)                 :: tvelo(3)                         !  u_tan
  real(rp)                 :: tvefi(3)                         !  u_fix_tan  - it comes from bvess_nsi put in global system
  real(rp)                 :: tvedi(3)                         !  u_tan - u_fix_tan
  real(rp)                 :: velsh(ndime,pnodb*ndime)
  real(rp)                 :: delta_aux

  if( kfl_delta == 1 ) then
     delta_aux = ywalb(iboun)                                  ! variable wall distance
  else
     delta_aux = delta_dom                                     ! fixed wall distance
  end if

  if( ( delta_aux > zensa ) .or. ( kfl_delta == 1 ) ) then

     do idime = 1,ndime                                        ! Velocity
        veloc(idime) = 0.0_rp
        velfi(idime) = 0.0_rp                                  ! velocity of the wall, due to moving domain
     end do
     do inodb = 1,pnodb
        do idime = 1,ndime            
           veloc(idime) = veloc(idime) + gbsha(inodb) * bovel(idime,inodb)
           velfi(idime) = velfi(idime) + gbsha(inodb) * bovfi(idime,inodb)
        end do
     end do
     !
     ! Tangent velocity (TVELO), tangent component of the prescribed velocity (TVEFI) and TVENO = |TVELO-TVEFI| 
     !
     do idime = 1,ndime     
        tvelo(idime) = veloc(idime)
        tvefi(idime) = velfi(idime)
        do ldime = 1,ndime
           tvelo(idime) = tvelo(idime)   &
                - baloc(idime,ndime) &
                * baloc(ldime,ndime) * veloc(ldime)
           tvefi(idime) = tvefi(idime)   &
                - baloc(idime,ndime) &
                * baloc(ldime,ndime) * velfi(ldime)
        end do
        tvedi(idime) = tvelo(idime) - tvefi(idime)
     end do
     !
     ! Compute U*: VELFR
     !
     call vecnor(tvedi,ndime,tveno,2_ip)                    ! |u_tan-u_fix_tan|
     vikin = gbvis/gbden                                    ! nu
     call frivel(delta_aux,rough,tveno,vikin,velfr)         ! U*
     !
     ! Compute prescribed traction
     !
     yplus = delta_aux * velfr / vikin

     if( yplus < 5.0_rp .and. kfl_ustar == 0 ) then

        fact1 = gbden * vikin / delta_aux                ! t = - mu/y * u

     else if( kfl_ustar == 1 ) then

        velfr = 0.41_rp / log( (delta_aux+rough) / rough )
        fact1 = gbden * tveno * velfr * velfr

     else

        fact1 = gbden * velfr * velfr / tveno            ! t = - rho*U*^2/|u_tan-u_fix_tan| * (u_tan-u_fix_tan)

     end if

     do idime = 1,ndime
        tract(idime) = tract(idime) - fact1 * tvedi(idime) 
     end do

  end if

end subroutine nsa_bouwal

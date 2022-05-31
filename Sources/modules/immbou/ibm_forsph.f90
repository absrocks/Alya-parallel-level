subroutine ibm_forsph()
  !-----------------------------------------------------------------------
  !****f* ibm_forsph/ibm_forsph
  ! NAME
  !    ibm_forsph
  ! DESCRIPTION
  !    Compute the drag force using a drag formula
  !    Fd = pi * mu * d/8 * (Cd*Re) * u
  ! USES
  !    dragfo
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_immbou
  implicit none
  integer(ip)       :: iimbo,idime,inode,ipoin,pelty,dummi,ielem
  integer(ip)       :: kfl_veloc
  real(rp)          :: shaib(mnode),v(3),deriv(3*64)
  real(rp)          :: coloc(3),Cd,CdRe,Du,diame,Re,x(3),dCddRe
  real(rp), pointer :: up(:,:),Fv(:),Fp(:) 
  !
  ! KFL_VELOC: Check if a velocity is available from other module
  !
  if( kfl_modul(ID_NASTIN) /= 0 .or. kfl_modul(ID_NASTAL) /= 0 .or. kfl_vefun > 0 ) then
     kfl_veloc = 1
  else
     kfl_veloc = 0
  end if
  call memgen(1_ip,nimbo,0_ip)
  !
  ! Initialize Fv and Fp
  !
  do iimbo = 1,nimbo
     Fv    => imbou(iimbo) % vforce
     Fp    => imbou(iimbo) % pforce
     do idime = 1,ndime
        Fv(idime) = 0.0_rp
        Fp(idime) = 0.0_rp
     end do
  end do
  !
  ! Loop over particles
  !
  if( INOTMASTER ) then

     do iimbo = 1,nimbo

        Fv    => imbou(iimbo) % vforce
        Fp    => imbou(iimbo) % pforce
        up    => imbou(iimbo) % velol
        diame =  ( 6.0_rp / pi * imbou(iimbo) % volum ) ** ( 1.0_rp / 3.0_rp )
        v(1) =  0.0_rp
        v(2) =  0.0_rp
        v(3) =  0.0_rp
        do idime = 1,ndime
           x(idime) = imbou(iimbo) % posil(idime,1) + imbou(iimbo)%posgr(idime)
        end do

        if( kfl_veloc == 1 ) then
       
           call runend('IBM NOT CODED')
           !call elsest(&
           !     2_ip,1_ip,ielse,mnode,ndime,npoin,nelem,nnode(1:),&
           !     lnods,ltype,ltopo,coord,x,relse,ielem,&
           !     shaib,deriv,coloc,dummi)

           if( ielem > 0 ) then
              !
              ! Fluid velocity at particle c.o.g.
              !    
              gisca(iimbo) = 1
              pelty = ltype(ielem) 
              if( pelty < 0 ) call runend('IBM_FORSPH: ERROR!')
              do inode = 1,nnode(pelty)
                 ipoin = lnods(inode,ielem)
                 do idime = 1,ndime
                    v(idime) = v(idime) + veloc(idime,ipoin,1) * shaib(inode)
                 end do
              end do

           else
              !
              ! This subdomain does not have the particle c.o.g 
              !
              gisca(iimbo) = 0

           end if

        else

           gisca(iimbo) = 1

        end if

        if( gisca(iimbo) == 1 ) then
           !
           ! Velocity difference Du = u -us
           !
           Du = 0.0_rp
           do idime = 1,ndime
              Du = Du + ( up(idime,1) - v(idime) ) ** 2
           end do
           Du = sqrt( Du )
           !
           ! Drag Cd
           ! 
           call dragfo(kfl_drafo_ibm,Du,visme,denme,diame,1.0_rp,CdRe,Cd,Re,dCddRe)
           !
           ! Force
           !
           do idime = 1,ndime
              Fv(idime) = -pi * visme * diame/8.0_rp  * CdRe * ( up(idime,1) - v(idime) )
              Fp(idime) =  0.0_rp
           end do

        end if

     end do

  else
     
     do iimbo = 1,nimbo
        gisca(iimbo) = 0
     end do

  end if
  !
  ! Average viscous force over subdomains
  !
  if( kfl_veloc == 1 ) then

     do iimbo = 1,nimbo
        if( gisca(iimbo) == 0 ) then
           call runend('IBM_FORSPH: COULD NOT FIND VELOCITY TO INTERPOLATE')
        end if
     end do

     if( IPARALL ) then
        call ibm_defini(9_ip)
        call parari('SUM',0_ip,nimbo,gisca)
        do iimbo = 1,nimbo
           Fv => imbou(iimbo) % vforce
           do idime = 1,ndime
              Fv(idime) = Fv(idime) / real(gisca(iimbo))
           end do
        end do
     end if

  end if

  call memgen(3_ip,nimbo,0_ip)

end subroutine ibm_forsph

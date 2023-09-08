subroutine qua_elmpre(itask,&
     pnode,pgaus,gpden,gpdif,gpsha,gpcar,elpro,&
     elcod,gprhs,gpcod,gpcou,gpaxi,lnods,gposc,&
     gpelf,gpesp,gppxc,gphar,gpion)
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_elmpre
  ! NAME
  !   qua_elmpre
  ! DESCRIPTION
  !    Compute some Gauss values
  ! OUTPUT 
  !    GPRHS .......... righ side contribution by element
  ! USES
  ! USED BY
  !    qua_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  mnode,mgaus,ndime,ntens,kfl_naxis,kfl_spher
  use def_quanty

  implicit none
  integer(ip), intent(in)    :: itask,pnode,pgaus,lnods(pnode)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus),gpdif(pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: elpro(mnode)
  real(rp),    intent(in)    :: elcod(ndime,pnode)
  real(rp),    intent(inout) :: gpden(pgaus)
  real(rp),    intent(inout) :: gpcou(pgaus),gposc(pgaus),gpaxi(pgaus),gpelf(pgaus),gpesp(pgaus)
  real(rp),    intent(out)   :: gpcod(ndime,pgaus)
  real(rp),    intent(out)   :: gprhs(pgaus)
  real(rp),    intent(out)   :: gppxc(pgaus),gphar(pgaus),gpion(pgaus)
  integer(ip)                :: idime,inode,igaus,itime,ipoin
  real(rp)                   :: dummr,gpdkt,rho_pro,Har_pro,xc_pro,ion_pro,aa_1,bb_1
  real(rp)                   :: ax_1,cc_1,xcoef_1,ycoef_1,zcoef_1,thrd_1
  !
  ! Coordinates
  !
  do igaus = 1,pgaus
     do idime = 1,ndime
        gpcod(idime,igaus) = 0.0_rp
     end do
     do inode = 1,pnode
        do idime = 1,ndime
           gpcod(idime,igaus) = gpcod(idime,igaus) &
                + gpsha(inode,igaus) * elcod(idime,inode)
        end do
     end do
  end do
  !
  ! termino de fuente para DFT u all electron case
  !
  ! Right-hand side: GPRHS
  !
  if( itask == 1 ) then
     rho_pro = 0.0_rp
     do inode = 1,pnode
        rho_pro = rho_pro + elpro(inode)
     enddo
     rho_pro = rho_pro/real(pnode)

     do igaus = 1,pgaus
        gprhs(igaus) = 16.0_rp * atan(1.0_rp) * rho_pro
     end do
  else
     do igaus = 1,pgaus
        gprhs(igaus) = 0.0_rp
     end do
  endif
  !
  ! coulomb term: GPCOU
  !
  if(kfl_coulo_qua==1 ) then
     if(ndime==3) then
        do igaus=1,pgaus
           gpcou(igaus)= coulo_qua/(sqrt( (gpcod(1,igaus)-x0_qua)**2 + &
                (gpcod(2,igaus)-y0_qua)**2+(gpcod(3,igaus)-z0_qua)**2 ))    
        end do
     elseif(ndime==2) then
        do igaus=1,pgaus
           gpcou(igaus)= coulo_qua/(sqrt( (gpcod(1,igaus)-x0_qua)**2 + &
                (gpcod(2,igaus)-y0_qua)**2))    
        end do
     elseif(ndime==1) then
        do igaus=1,pgaus
           gpcou(igaus)= coulo_qua/sqrt( (gpcod(1,igaus)-x0_qua)**2)    
        end do
     endif
  else
     do igaus = 1,pgaus
        gpcou(igaus) = 0.0_rp
     end do
  end if
  !
  ! axisym term ON
  !
  if( kfl_naxis == 1 ) then

     do igaus = 1,pgaus
        gpaxi(igaus) = 0.5_rp*lcuanorb_qua**2/(gpcod(1,igaus)-x0_qua)**2    
     end do

  else
     do igaus = 1,pgaus
        gpaxi(igaus) = 0.0_rp
     end do
  end if
  !
  ! espheric term ON
  !
  if( kfl_spher == 1 ) then

     do igaus = 1,pgaus
        gpesp(igaus) = lcuanorb_qua*(lcuanorb_qua+1)*0.5_rp/(gpcod(1,igaus)-x0_qua)**2    
     end do

  else

     do igaus = 1,pgaus
        gpesp(igaus) = 0.0_rp
     end do

  end if

  if( itask == 0 ) then
     if(kfl_dftgs_qua/=0 .or. kfl_alele_qua/=0) then

        ! aca se agrega el termino de potencial nuclear en caso DFT
        ! en caso all electron esto seria como coulomb  
        rho_pro = 0.0_rp
        do inode = 1,pnode
           rho_pro = rho_pro + elpro(inode)
        enddo
        rho_pro = rho_pro / real(pnode)


        ion_pro = 0.0_rp
        do inode = 1,pnode
           ipoin = lnods(inode)
           ion_pro = ion_pro + v_pot_ps(ipoin)
        enddo
        ion_pro = ion_pro / real(pnode)
        do igaus = 1,pgaus
           gpion(igaus) = ion_pro
        end do

        ! for test!!
        !      do igaus=1,pgaus
        !          gpion(igaus)= -3.0/sqrt( (gpcod(1,igaus)-x0_qua)**2)    
        !      enddo

        ! Hartree y XC
        Har_pro = 0.0_rp
        do inode = 1,pnode
           ipoin = lnods(inode)
           Har_pro = Har_pro + v_hartree(ipoin)
        enddo
        Har_pro = Har_pro / real(pnode)

        do igaus=1,pgaus
           gphar(igaus)=Har_pro
        end do

        if(kfl_potxc_qua==0) then 

           do igaus=1,pgaus
              gppxc(igaus)=0.0
           end do

        elseif(kfl_potxc_qua==1) then
           ! EX-CHANGE-correlation
           thrd_1 = 0.333333333333333333_rp
           aa_1   = 0.93222_rp
           bb_1   = 0.00947362_rp
           ax_1   = 0.738558766382022406_rp

           xcoef_1 = bb_1*rho_pro**thrd_1
           ycoef_1 = xcoef_1/(xcoef_1 + 1.0_rp)
           zcoef_1 = -aa_1*xcoef_1/bb_1

           do igaus=1,pgaus
              gppxc(igaus) = thrd_1 * zcoef_1 * (4.0_rp + 5.0_rp*xcoef_1*log(ycoef_1) + ycoef_1)
           end do
           ! esto es electron-gas based LDA exchange
           !         VXCALFA =-(3/PI)**0.33333333 *alfa* RHOPRO**0.33333333   !(* 0.7385587*(4./3.))
        endif

     else    

        do igaus = 1,pgaus
           gphar(igaus) = 0.0_rp
           gpion(igaus) = 0.0_rp
           gppxc(igaus) = 0.0_rp
        end do

     endif
  else

     do igaus = 1,pgaus
        gphar(igaus) = 0.0_rp
        gpion(igaus) = 0.0_rp
        gppxc(igaus) = 0.0_rp
     end do
  end if


  ! electric field x On
  do igaus=1,pgaus
     gpelf(igaus)= 0.0_rp    
  end do

  if(kfl_efieldx_qua==1) then

     if( ndime == 3 ) then
        do igaus = 1,pgaus
           gpelf(igaus) = gpelf(igaus) - efieldx_qua * ( gpcod(1,igaus) - x0_qua )    
        end do
     else if( ndime == 2 ) then
        do igaus = 1,pgaus
           gpelf(igaus) = gpelf(igaus) - efieldx_qua * ( gpcod(1,igaus) - x0_qua )
        end do
     endif

  endif


  ! electric field y On

  if(kfl_efieldy_qua==1) then
     if(ndime==3) then
        do igaus=1,pgaus
           gpelf(igaus)=gpelf(igaus) -efieldy_qua* (gpcod(2,igaus)-y0_qua)    
        end do
     elseif(ndime==2) then
        do igaus=1,pgaus
           gpelf(igaus)= gpelf(igaus)-efieldy_qua* (gpcod(2,igaus)-y0_qua)
        end do
     endif

  endif

  ! electric field z On

  if(kfl_efieldz_qua==1) then
     if(ndime==3) then
        do igaus=1,pgaus
           gpelf(igaus)=gpelf(igaus) -efieldz_qua* (gpcod(3,igaus)-z0_qua)    
        end do
     endif
  endif



  ! other potential term (por ahora es oscilador armonico): GPOSC
  !
  if(kfl_vother_qua==1) then
     if(ndime==3) then
        do igaus=1,pgaus
           gposc(igaus)= vother_qua*0.5_rp*( (gpcod(1,igaus)-x0_qua)**2 + &
                (gpcod(2,igaus)-y0_qua)**2+(gpcod(3,igaus)-z0_qua)**2 )

           ! agrego termino de acomplamiento
           gposc(igaus)=gposc(igaus)+ w_vother_qua*((gpcod(1,igaus)-gpcod(2,igaus))**2 + &
                (gpcod(2,igaus)-gpcod(3,igaus))**2 )

        end do
     elseif(ndime==2) then
        do igaus=1,pgaus
           gposc(igaus)= vother_qua*0.5_rp*( (gpcod(1,igaus)-x0_qua)**2 + &
                (gpcod(2,igaus)-y0_qua)**2 )
           ! agrego termino de acomplamiento
           gposc(igaus)=gposc(igaus)+ w_vother_qua* (gpcod(1,igaus)-gpcod(2,igaus))**2     
        end do
     endif
  else
     do igaus=1,pgaus
        gposc(igaus)=0.0_rp
     end do
  end if

end subroutine qua_elmpre

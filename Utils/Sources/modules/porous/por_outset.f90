subroutine por_outset()
  !-----------------------------------------------------------------------
  !****f* porous/por_outset
  ! NAME 
  !    por_outset
  ! DESCRIPTION
  !    Compute and write results on sets
  ! USES
  ! USED BY
  !    npor_output
  !***
  !----------------------------------------------------------------------- 
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_porous
  use mod_memory
  implicit none
  integer(ip)         :: inset,ipoin,iwell,dummi
  real(rp)            :: auxii,auxiw,auxio,pbh,qw,qo,qt,wi,Bo
  real(rp)            :: Bw,deltp,xx_o,xx_w,p0,Sw
  real(rp),   pointer :: smokr(:,:) 
 
  !----------------------------------------------------------------------
  !
  ! Element sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setse) > 0 ) then
  end if

  !----------------------------------------------------------------------
  !
  ! Boundary sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setsb) > 0 ) then
  end if

  !----------------------------------------------------------------------
  !
  ! Node sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setsn) > 0 ) then

     if( INOTMASTER ) then

        nullify(smokr)
        call memory_alloca(mem_modul(1:2,modul),'SMOKR','por_outset',smokr,2_ip,npoin) 
        call por_smokrx(smokr)

        do inset = 1,nnset

           ipoin = lnsec(inset)
           if( ipoin /= 0 ) then
              iwell = iwell_por(ipoin)
              if (iwell /= 0 ) then

                 p0    = press(ipoin,1)
                 deltp = p0 - prref_por
                 xx_o  = comoi_por * deltp
                 xx_w  = comwa_por * deltp
                 Bo    = boref_por / (1.0_rp + xx_o + (0.5_rp*xx_o*xx_o) )
                 Bw    = bwref_por / (1.0_rp + xx_w + (0.5_rp*xx_w*xx_w) )
                 Sw    = wasat(ipoin,1)

                 if( tywel_por(iwell) % itype == 3_ip .or. tywel_por(iwell) % itype == 4_ip ) then  
                    !
                    ! Pbh
                    !
                    wi    = winde_por(ipoin)
                    pbh   = xwell_por(ipoin)
                    auxii = wi * ( pbh - press(ipoin,1) )
                    qw    = auxii*denwa_por*smokr(1,ipoin)/(muwat_por*Bw)
                    qo    = auxii*denoi_por*smokr(2,ipoin)/(muoil_por*Bo)
                    qt    = qo + qw

                 else if( tywel_por(iwell) % itype == 1_ip ) then    
                    !  
                    ! Water flow rate prescribed 
                    !
                    wi    = winde_por(ipoin)
                    qw    = xwell_por(ipoin)
                    qt    = qw
                    qo    = 0.0_rp
                    pbh   = press(ipoin,1) + ( qw * muwat_por*Bw / ( wi * denwa_por*smokr(1,ipoin) ) )

                 else if( tywel_por(iwell) % itype == 2_ip ) then      
                    !
                    ! Total flow rate prescribed
                    !
                    wi    = winde_por(ipoin)
                    qt    = xwell_por(ipoin)
                    auxio = denoi_por*smokr(2,ipoin)/(muoil_por*Bo) 
                    auxiw = denwa_por*smokr(1,ipoin)/(muwat_por*Bw)
                    auxii = 1.0_rp / ( auxio + auxiw )
                    pbh   = press(ipoin,1) + qt * auxii / wi
                    qw    = qt * auxiw * auxii 
                    qo    = qt * auxio * auxii 

                 else if( tywel_por(iwell) % itype == 5_ip ) then   
                    !   
                    ! Oil flow rate prescribed at well 
                    !
                    wi    = winde_por(ipoin)
                    qo    = xwell_por(ipoin)

                    auxio = denoi_por*smokr(2,ipoin) / (muoil_por*Bo) 
                    auxiw = denwa_por*smokr(1,ipoin) / (muwat_por*Bw)
                    auxii = 1.0_rp / ( auxio + auxiw )

                    qw    = qo * auxiw / auxio
                    qt    = qw + qo
                    pbh   = press(ipoin,1) + qt * auxii / wi

                 end if

                 if( postp(1) % npp_setsn(1) /= 0 ) vnset(1,inset) = wi     ! WI
                 if( postp(1) % npp_setsn(2) /= 0 ) vnset(2,inset) = qt     ! QT
                 if( postp(1) % npp_setsn(3) /= 0 ) vnset(3,inset) = qw     ! QW
                 if( postp(1) % npp_setsn(4) /= 0 ) vnset(4,inset) = qo     ! QO
                 if( postp(1) % npp_setsn(5) /= 0 ) vnset(5,inset) = pbh    ! PBH
                 if( postp(1) % npp_setsn(6) /= 0 ) vnset(6,inset) = pbh-p0 ! DELTP
                 if( postp(1) % npp_setsn(7) /= 0 ) vnset(7,inset) = Sw     ! SW
                 if( postp(1) % npp_setsn(8) /= 0 ) vnset(8,inset) = real(iwell_por(ipoin),rp)     ! Well number

              end if
           end if
        end do
        call memory_deallo(mem_modul(1:2,modul),'SMOKR','por_outset',smokr) 

     end if
     call posdef(23_ip,dummi)

  end if

end subroutine por_outset

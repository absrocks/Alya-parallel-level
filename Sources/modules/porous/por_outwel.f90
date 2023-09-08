!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_outwel.f90
!> @date    29/08/2013
!> @author  Herbert Owen
!> @brief   Output results (flow rate and Pbh) at well  -- no longer used guillaume created por_outset
!> @details Output results (flow rate and Pbh) at well
!> @} 
!------------------------------------------------------------------------
subroutine por_outwel()
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  mem_modul,cutim,modul,press,kfl_paral,wasat
  use def_domain, only       :  npoin
  use def_porous, only       :  iwell_por,iheip_por,boref_por,bwref_por,comoi_por,comwa_por,denoi_por,denwa_por,&
                                muoil_por,muwat_por,nwell_por,prref_por,kfl_wellc_por,winde_por,dataw_por,xwell_por,tywel_por
  use mod_memory 

  implicit none
  integer(ip)         :: ipoin,iheiw,iwell
  real(rp)            :: auxii,auxiw,auxio,pbh,qw,qo,qt,wi,Bo,Bw,deltp,xx_o,xx_w,p0,Sw
  real(rp),pointer    :: smokr(:,:)   ! perhaps I could send to def_porous and usar the one I have already calculated in por_wellte

  ! for paralell we surelly can do it tidier but I guess this should work

  if (nwell_por==0) return
  nullify(smokr)
  call memory_alloca(mem_modul(1:2,modul),'smokr','por_wellte',smokr,2_ip,npoin) 
  !
  !  Smooth kro & krw 
  !
  call por_smokrx(smokr)
 
  if (cutim<0.00001_rp)write(kfl_paral+900,'(a)')' ipoin   iwell   iheiw           wi           qt           qw           qo          pbh      pbh-p0       Sw       cutim'
  do ipoin=1,npoin
     iwell = iwell_por(ipoin)
     iheiw = 0
     if (iwell /= 0 ) then
        p0 = press(ipoin,1)
        deltp = p0 - prref_por
        xx_o = comoi_por * deltp
        xx_w = comwa_por * deltp
        Bo = boref_por / (1.0_rp + xx_o + (0.5_rp*xx_o*xx_o) )
        Bw = bwref_por / (1.0_rp + xx_w + (0.5_rp*xx_w*xx_w) )

        if ( (tywel_por(iwell) % itype == 3_ip) .or. (tywel_por(iwell) % itype == 4_ip) ) then  ! Pbh

           wi = winde_por(ipoin)
           pbh = xwell_por(ipoin)
           auxii = wi * ( pbh - press(ipoin,1) )
           qw = auxii*denwa_por*smokr(1,ipoin)/(muwat_por*Bw)
           qo = auxii*denoi_por*smokr(2,ipoin)/(muoil_por*Bo)
           qt = qo + qw

        else if (tywel_por(iwell) % itype == 1_ip) then      ! Water flow rate prescribed 

           wi = winde_por(ipoin)
           qw = xwell_por(ipoin)
           qt = qw
           qo = 0.0_rp
           pbh = press(ipoin,1) + ( qw * muwat_por*Bw / ( wi * denwa_por*smokr(1,ipoin) ) )

        else if (tywel_por(iwell) % itype == 2_ip) then      ! Total flow rate prescribed

           wi = winde_por(ipoin)
           qt = xwell_por(ipoin)
           auxio = denoi_por*smokr(2,ipoin)/(muoil_por*Bo) 
           auxiw = denwa_por*smokr(1,ipoin)/(muwat_por*Bw)
           auxii = 1.0_rp / ( auxio + auxiw )
           pbh = press(ipoin,1) + qt * auxii / wi
           qw = qt * auxiw * auxii 
           qo = qt * auxio * auxii 

        else if (tywel_por(iwell) % itype == 5_ip) then      ! Oil flow rate prescribed at well 

           Sw = wasat(ipoin,1)
           wi = winde_por(ipoin)
           qo = xwell_por(ipoin)
           qw = qo * ( denwa_por * smokr(1,ipoin) * muoil_por * Bo ) / ( denoi_por * smokr(2,ipoin) * muwat_por * Bw )
           qt = qw + qo

           auxio = denoi_por*smokr(2,ipoin)/(muoil_por*Bo) 
           auxiw = denwa_por*smokr(1,ipoin)/(muwat_por*Bw)
           auxii = 1.0_rp / ( auxio + auxiw )
           pbh = press(ipoin,1) + qt * auxii / wi


        end if
        write(kfl_paral+900,'(a,3(i7,1x),8(e12.5,1x))')'A',ipoin,iwell,iheiw,wi,qt,qw,qo,pbh,pbh-p0,wasat(ipoin,1),cutim
        flush(kfl_paral+900)
!        write(965,'(a,3(i7,1x),4(e12.5,1x))')'B',ipoin,iwell,iheiw,denoi_por,smokr(2,ipoin),muoil_por,Bo
!        write(965,'(a,3(i7,1x),4(e12.5,1x))')'C',ipoin,iwell,iheiw,denwa_por,smokr(1,ipoin),muwat_por,Bw
     end if
  end do

  call memory_deallo(mem_modul(1:2,modul),'SMOKR','por_wellte',smokr )

end subroutine por_outwel

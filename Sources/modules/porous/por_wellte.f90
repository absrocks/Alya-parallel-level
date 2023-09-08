!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_wellte.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Adds Well Terms
!> @details Adds Well Terms
!! for the pressure : 
!! 
!> @} 
!------------------------------------------------------------------------
subroutine por_wellte()
  ! 
  use def_kintyp, only :  ip,rp
  use def_domain, only :  npoin,coord
  use def_porous, only :  winde_por,iwell_por,prref_por,comwa_por,&
       &                  bwref_por,muwat_por,muoil_por,kprsa_por,&
       &                  nwell_por,boref_por,comoi_por,denoi_por,&
       &                  denwa_por,kfl_wellc_por,kfl_malta_por,&
       &                  tywel_por,xwell_por,wmass_por
  use def_master, only :  ID_POROUS,mem_modul,modul,press,amatr,&
       &                  rhsid,solve_sol,cutim,wasat,IMASTER,&
       &                  npoi1,npoi2,npoi3
  use mod_memory, only :  memory_alloca,memory_deallo 
  use mod_por_interp
  implicit none
  integer(ip)          :: ipoin,iwell,izmat,ntime,nposi
  real(rp)             :: deltp,xx_w,Bwinv,wellt,xvalu
  real(rp)             :: praux,xx_o,Bo,Bw
  real(rp)             :: auxii,f_w,lambd_o,lambd_w
  real(rp),    pointer :: smokr(:,:)
  real(rp),    pointer :: table_value(:,:)
  real(rp),    pointer :: table_coord(:)
  integer(ip), pointer :: lpoin(:)

  if( nwell_por == 0 .or. IMASTER ) return
  !
  ! Allocate memory
  !
  nullify(smokr)
  nullify(lpoin)
  call memory_alloca(mem_modul(1:2,modul),'SMOKR','por_wellte',smokr,2_ip,npoin) 
  call memory_alloca(mem_modul(1:2,modul),'LPOIN','por_wellte',lpoin,npoin) 
  !
  ! Decide where to compute source term
  !
  do ipoin = 1,npoi1
     lpoin(ipoin) = iwell_por(ipoin)
  end do
  do ipoin = npoi2,npoi3
     lpoin(ipoin) = iwell_por(ipoin)
  end do
  !
  ! Smooth kro & krw 
  !
  call por_smokrx(smokr)
  !
  ! Obtain well index per node
  !
  do ipoin = 1,npoin
     iwell = lpoin(ipoin)
     if( iwell > 0 ) then
        ntime       =  tywel_por(iwell) % ntime
        nposi       =  tywel_por(iwell) % nposi
        if( tywel_por(iwell) % itype == 1_ip .or. tywel_por(iwell) % itype == 2_ip .or. tywel_por(iwell) % itype == 5_ip ) then 
           !
           ! Water/oil/Total flow rate
           !
           table_value => tywel_por(iwell) % q_table
           call por_interp(ntime,nposi,coord(3,ipoin),cutim,xvalu,table_value)
           xwell_por(ipoin) = xvalu
           xwell_por(ipoin) = xwell_por(ipoin) * wmass_por(ipoin)

        else if( tywel_por(iwell) % itype == 3_ip ) then 
           !
           ! Interpolate Pbh 
           !
           table_value => tywel_por(iwell) % pbh_table
           table_coord => tywel_por(iwell) % pbh_coord
           call por_interp(ntime,nposi,coord(3,ipoin),cutim,xvalu,table_value,table_coord)
           xwell_por(ipoin) = xvalu

        else if( tywel_por(iwell) % itype == 4_ip ) then 
           !
           ! pbh given at middle height
           !
           call runend('DE MOMENTO NADA')
        end if
     end if
  end do
  !
  ! Assemble matrix and RHS
  !
  do ipoin = 1,npoin
     iwell = lpoin(ipoin)
     if( iwell > 0 ) then

        if ( tywel_por(iwell) % itype == 3_ip .or. tywel_por(iwell) % itype == 4_ip ) then  ! Pbh

           deltp = press(ipoin,1) - prref_por
           xx_w  = comwa_por * deltp
           Bwinv = (1.0_rp + xx_w + (0.5_rp*xx_w*xx_w) )  /  bwref_por   ! 1/Bw

           if ( kprsa_por == 1_ip ) then 
              !
              ! Pressure
              !
              wellt = (smokr(1,ipoin) / muwat_por) + (smokr(2,ipoin) / muoil_por)
           else 
              !
              ! Saturation
              !
              if ( kfl_malta_por == 0 ) then 
                 wellt = (smokr(1,ipoin) / muwat_por)
              else if ( kfl_malta_por == 1 ) then
                 wellt = 0.0_rp   ! They end up cancelling out see 2phase..lyx
              end if
           end if

           wellt = wellt * winde_por(ipoin) * Bwinv
           praux = xwell_por(ipoin)

           if ( kprsa_por == 1_ip ) then  
              !
              ! Pressure equation 
              !
              rhsid (ipoin) = rhsid(ipoin) + wellt * praux
              call csrdia(ipoin,solve_sol(1)%kfl_symme,izmat)   ! Here we should have a non symmetric solver for both pr and sat??? 
              amatr(izmat) = amatr(izmat) + wellt
           else  
              ! Saturation
              !
              rhsid (ipoin) = rhsid(ipoin) + denwa_por * wellt * ( praux - press(ipoin,1) ) 
           end if


        else if (tywel_por(iwell) % itype == 1_ip) then      ! Water flow rate prescribed 
           !
           ! Beware: for the moment we are solving the Saturation equation including rho_ws
           ! Instead in the pressure equation everything is divided by density
           ! The flow rate we prescribe is  a mass flow rate and therefore it needs to be divided by denwa_por for the pressure eq.
           ! Perhaps at some time divide all the density equation by denwa_por.
           !
           if (kprsa_por == 1_ip) then  
              !
              ! Pressure equation
              !
              rhsid (ipoin) = rhsid(ipoin) + ( xwell_por(ipoin) / denwa_por )
           else
              if ( kfl_malta_por == 0 ) then 
                 rhsid (ipoin) = rhsid(ipoin) + xwell_por(ipoin)
              else if ( kfl_malta_por == 1 ) then
                 !
                 ! Obtain f_w
                 !
                 lambd_o = smokr(2,ipoin) / muoil_por
                 lambd_w = smokr(1,ipoin) / muwat_por
                 f_w = lambd_w  / (lambd_w + lambd_o)
                 rhsid (ipoin) = rhsid(ipoin) + ( 1.0_rp - f_w ) * xwell_por(ipoin)   ! They end up cancelling out see 2phase..lyx
              end if
           end if

        else if (tywel_por(iwell) % itype == 2_ip) then      ! Total flow rate prescribed

           call runend('por_wellte:The units are not ok for pressure - needs to be revised')

           deltp = press(ipoin,1) - prref_por
           xx_o = comoi_por * deltp
           xx_w = comwa_por * deltp
           Bo = boref_por / (1.0_rp + xx_o + (0.5_rp*xx_o*xx_o) )
           Bw = bwref_por / (1.0_rp + xx_w + (0.5_rp*xx_w*xx_w) )
           !
           ! auxii = 1/( (rho_w*krw/mu_w) + (rho_o*kro/mu_o))
           !
           auxii = (denwa_por*smokr(1,ipoin)/(muwat_por*Bw) ) +  (denoi_por*smokr(2,ipoin)/(muoil_por*Bo) )
           auxii = 1.0_rp / auxii

           if (kprsa_por == 1_ip) then  ! pressure equation
              rhsid (ipoin) = rhsid(ipoin) + ( ( xwell_por(ipoin) * auxii ) * &
                   ( ( smokr(1,ipoin)/muwat_por ) + ( smokr(2,ipoin)/muoil_por ) ) / Bw )
           else
              rhsid (ipoin) = rhsid(ipoin) + xwell_por(ipoin) * auxii * (denwa_por*smokr(1,ipoin) / (muwat_por*Bw) )
           end if

        else if (tywel_por(iwell) % itype == 5_ip ) then      ! Oil flow rate prescribed - only for well

           if (kprsa_por == 1_ip) then  ! pressure equation
              deltp = press(ipoin,1) - prref_por
              xx_o = comoi_por * deltp
              xx_w = comwa_por * deltp
              Bo = boref_por / (1.0_rp + xx_o + (0.5_rp*xx_o*xx_o) )
              Bw = bwref_por / (1.0_rp + xx_w + (0.5_rp*xx_w*xx_w) )

              auxii = ( Bo / ( Bw * denoi_por )  ) * &
                   ( 1.0_rp + ( ( smokr(1,ipoin) * muoil_por ) / ( smokr(2,ipoin) * muwat_por ) ) )
              rhsid (ipoin) = rhsid(ipoin) + ( xwell_por(ipoin) * auxii )
           else
              if ( kfl_malta_por /= 1 ) call runend ('por_wellte: oil flow rate prescribed only ready in kfl_malta_por == 1') 
              ! in the malta case nothing needs to be added to the RHS for saturation
           end if
        else
           call runend('POR_WELLTE: kfl_wellc_por must be 1,2,3,4 or 5')
        end if
     end if
  end do

  call memory_deallo(mem_modul(1:2,modul),'LPOIN','por_wellte',lpoin)
  call memory_deallo(mem_modul(1:2,modul),'SMOKR','por_wellte',smokr)

end subroutine por_wellte

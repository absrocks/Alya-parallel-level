subroutine chm_nsacfi()
  !-----------------------------------------------------------------------
  !****f* partis/chm_heatso
  ! NAME 
  !    chm_nsacfi
  ! DESCRIPTION
  !    CFI combustion model: before leaving chemic, send chemical 
  !    heat source to gauss points as required in nastal
  ! USES
  ! USED BY
  !   
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      : ip,rp
  use def_master, only      : inotmaster,radiative_heat,wmean,prthe
  use def_domain, only      : nelem,ltype,nnode,ngaus,mnode,lnods,elmar,mgaus
  use def_chemic, only      : tempe_chm,kfl_radia_chm,rspec_chm,radwt_chm
  implicit none
  integer(ip)               :: ipoin,ielem,inode,igaus,ispec,icoef
  integer(ip)               :: pelty,pnode,pgaus
  real(rp)                  :: gpspe(2,mgaus), gpwme(mgaus), gptem(mgaus)
  real(rp)                  :: rcoef(2,7), mwspe(2)
  real(rp)                  :: tempi, gprak, apk, gpram

  if (INOTMASTER) then

    !
    ! Simple radiation model based on CO2 and H2O mass fractions
    !
    if (kfl_radia_chm > 0) then

      mwspe(1) = 18.0153400897980e-3_rp
      mwspe(2) = 44.0099506378174e-3_rp
      rcoef(1,1) =  0.38041e1_rp
      rcoef(1,2) = -0.27808e1_rp
      rcoef(1,3) =  0.11672e1_rp
      rcoef(1,4) = -0.28491e0_rp
      rcoef(1,5) =  0.38162e-1_rp
      rcoef(1,6) = -0.26292e-2_rp
      rcoef(1,7) =  0.37774e-4_rp
      rcoef(2,1) =  0.22317e1_rp
      rcoef(2,2) = -0.15829e1_rp
      rcoef(2,3) =  0.13296e1_rp
      rcoef(2,4) = -0.50707e0_rp
      rcoef(2,5) =  0.93334e-1_rp
      rcoef(2,6) = -0.83108e-2_rp
      rcoef(2,7) =  0.28834e-3_rp
      !
      ! Loop over elements
      !
      elements: do ielem=1,nelem
        !
        ! Element properties and dimensions
        !
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        !
        ! Gather operations
        !
        gpspe = 0.0_rp
        gpwme = 0.0_rp
        gptem = 0.0_rp
        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           do igaus=1,pgaus
              gpspe(1,igaus)=gpspe(1,igaus)&
                    +elmar(pelty)%shape(inode,igaus)*rspec_chm(1,ipoin)
              gpspe(2,igaus)=gpspe(2,igaus)&
                   +elmar(pelty)%shape(inode,igaus)*rspec_chm(2,ipoin)
              gpwme(igaus)=gpwme(igaus)&
                   +elmar(pelty)%shape(inode,igaus)*wmean(ipoin,1)
              gptem(igaus)=gptem(igaus)&
                   +elmar(pelty)%shape(inode,igaus)*tempe_chm(ipoin)
           end do
        end do
        !
        ! Compute 
        !
        do igaus=1,pgaus
          tempi = gptem(igaus)/300.0_rp
          gprak = 0.0_rp
          do ispec = 1,2
            apk = rcoef(ispec,7)
            do icoef = 1,6
              apk = apk*tempi + rcoef(ispec,7-icoef)
            end do
            apk = 10.0_rp**apk
            gprak = gprak + apk * gpspe(ispec,igaus) * gpwme(igaus) / mwspe(ispec) * prthe(1)/101325.0_rp
          end do
          radiative_heat(ielem)%a(igaus) = -4.0_rp*5.6696e-8_rp*gptem(igaus)**4*gprak
        end do

        if (kfl_radia_chm == 2) then
          do igaus=1,pgaus
            tempi = radwt_chm/300.0_rp
            gprak = 0.0_rp
            do ispec = 1,2
              apk = rcoef(ispec,7)
              do icoef = 1,6
                apk = apk*tempi + rcoef(ispec,7-icoef)
              end do
              apk = 10.0_rp**apk
              gprak = gprak + apk * gpspe(ispec,igaus) * gpwme(igaus) / mwspe(ispec) * prthe(1)/101325.0_rp
            end do
            gpram = gprak * radwt_chm/gptem(igaus)
            radiative_heat(ielem)%a(igaus) = radiative_heat(ielem)%a(igaus) + 4.0_rp*5.6696e-8_rp*radwt_chm**4*gpram
          end do
        end if
    
      end do elements

    end if

  end if
  
end subroutine chm_nsacfi
